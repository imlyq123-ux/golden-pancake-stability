"""Independent multiprecision cross-check of the Golden Pancake orbits and Floquet data.

Integrator: tools/gbs_mpfr (Gragg-Bulirsch-Stoer extrapolation, MPFR), not the Taylor code.
Newton shooting, invariants and Floquet analysis: mpmath, sharing no code with
high_precision.py. The Floquet analysis here is reduction-free:

* the 8 trivial multipliers of the exact monodromy equal 1, so the nontrivial quartic
  follows from the power sums p_k = tr(Phi^k) - 8 (k = 1..4) by Newton's identities, and
  D = 2 p_2 + 8 - p_1^2 needs no symmetry reduction or eigenvector selection;
* all 12 eigenvalues and right eigenvectors are also computed with mpmath.eig;
* Krein forms use full-space eigenvectors and Omega = [[0, M], [-M, 0]], M = diag(1,1,mu,mu,1,1).
  Eigenvectors of nontrivial multipliers are Omega-orthogonal to the generalized 1-eigenspace,
  so their Krein signs agree with those of any symplectic quotient.

Seeds come from binary64 continuation (threebody_expand.py) of the originally supplied
sweep, or from an explicit JSON file; they never come from the Taylor/Decimal MP records.
"""
import argparse
import csv
import hashlib
import json
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import mpmath
from mpmath import mp, mpf, mpc

ROOT = Path(__file__).resolve().parent
GBS = ROOT / 'tools' / 'gbs_mpfr'
NAMES = ['q1x', 'q1y', 'q2x', 'q2y', 'q3x', 'q3y', 'v1x', 'v1y', 'v2x', 'v2y', 'v3x', 'v3y']
PAIRS = [(0, 1), (0, 2), (1, 2)]


def s(x, digits=None):
    """Decimal string of an mpmath number at the working precision."""
    return mpmath.nstr(x, digits or mp.dps, min_fixed=-5, max_fixed=5, strip_zeros=False)


def initial_state(mu, u, v):
    mu, u, v = mpf(mu), mpf(u), mpf(v)
    return [mpf(-1), mpf(0), mpf(0), mpf(0), mpf(1), mpf(0), u, v, -2 * u / mu, -2 * v / mu, u, v]


def acceleration(x, mu):
    m = [mpf(1), mpf(mu), mpf(1)]
    acc = [mpf(0)] * 6
    for i, j in PAIRS:
        dx, dy = x[2 * j] - x[2 * i], x[2 * j + 1] - x[2 * i + 1]
        g = (dx * dx + dy * dy) ** mpf(-1.5)
        acc[2 * i] += m[j] * dx * g; acc[2 * i + 1] += m[j] * dy * g
        acc[2 * j] -= m[i] * dx * g; acc[2 * j + 1] -= m[i] * dy * g
    return acc


def invariants(x, mu):
    m = [mpf(1), mpf(mu), mpf(1)]
    kinetic = sum(m[i] * (x[6 + 2 * i] ** 2 + x[7 + 2 * i] ** 2) for i in range(3)) / 2
    dist = [mpmath.sqrt((x[2 * j] - x[2 * i]) ** 2 + (x[2 * j + 1] - x[2 * i + 1]) ** 2) for i, j in PAIRS]
    energy = kinetic - sum(m[i] * m[j] / d for (i, j), d in zip(PAIRS, dist))
    angular = sum(m[i] * (x[2 * i] * x[7 + 2 * i] - x[2 * i + 1] * x[6 + 2 * i]) for i in range(3))
    com = [sum(m[i] * x[2 * i + k] for i in range(3)) / sum(m) for k in range(2)]
    momentum = [sum(m[i] * x[6 + 2 * i + k] for i in range(3)) for k in range(2)]
    return energy, angular, com, momentum, dist


def run_gbs(mu, p, digits, tol, stages, columns, path, stride):
    command = [str(GBS), str(digits), str(tol), str(stages), str(mu), *[str(x) for x in p], str(columns), str(path), str(stride)]
    started = time.perf_counter()
    done = subprocess.run(command, capture_output=True, text=True)
    elapsed = time.perf_counter() - started
    path.with_suffix('.log').write_text(' '.join(command) + '\n' + done.stdout + done.stderr)
    if done.returncode != 0:
        raise RuntimeError(f'GBS integration failed for mu={mu}: {done.stderr.strip()}')
    lines = path.read_text().splitlines()
    if not lines[-1].startswith('# success=1'):
        raise RuntimeError('incomplete GBS integration')
    rows = [line for line in lines if not line.startswith('#')]
    info = {'integrator': 'GBS extrapolation (tools/gbs_mpfr)', 'digits': digits, 'working_bits': lines[0],
            'tolerance': str(tol), 'stages': stages, 'columns': columns, 'stride': stride,
            'completion': lines[-1], 'runtime_seconds': elapsed,
            'raw_file': str(path.resolve().relative_to(ROOT)) if path.resolve().is_relative_to(ROOT) else str(path),
            'raw_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    return rows, info


def parse(row, count=None):
    values = row.split(',')
    return [mpf(v) for v in (values if count is None else values[:count])]


def newton(mu, seed, args, raw_dir, prefix):
    """Least-squares Newton on all 12 return components, unknowns (u, v, T)."""
    p = [mpf(x) for x in seed]
    history = []
    dd = [[mpf(0)] * 2 for _ in range(12)]
    dd[6] = [mpf(1), mpf(0)]; dd[7] = [mpf(0), mpf(1)]
    dd[8] = [-2 / mpf(mu), mpf(0)]; dd[9] = [mpf(0), -2 / mpf(mu)]
    dd[10] = [mpf(1), mpf(0)]; dd[11] = [mpf(0), mpf(1)]
    previous = None
    for iteration in range(args.max_newton):
        rows, info = run_gbs(mu, [s(x) for x in p], args.digits, args.newton_tol, args.stages, 2,
                             raw_dir / f'{prefix}_newton_{iteration}.csv', 10 ** 9)
        final = parse(rows[-1])
        x0 = initial_state(mu, p[0], p[1])
        residual = [a - b for a, b in zip(final[1:13], x0)]
        norm = max(abs(r) for r in residual)
        history.append({'iteration': iteration, 'param': [s(x) for x in p], 'max_residual': s(norm, 8),
                        'integration': info})
        print(f'{prefix} Newton {iteration}: residual {mpmath.nstr(norm, 5)}', flush=True)
        if norm < mpf(args.target):
            return p, history
        if previous is not None and norm > previous / 10:
            history[-1]['note'] = 'slow convergence'
        previous = norm
        field = final[7:13] + acceleration(final[1:13], mu)  # vector field (v, a) at X(T): d/dT column
        jac = [[final[13 + k] - dd[k][0], final[25 + k] - dd[k][1], field[k]] for k in range(12)]
        J = mpmath.matrix(jac)
        r = mpmath.matrix(residual)
        step = mpmath.lu_solve(J.T * J, -(J.T * r))
        scale = max(abs(step[0]) / mpf('0.002'), abs(step[1]) / mpf('0.002'), abs(step[2]) / mpf('0.04'))
        if scale > 1:
            raise RuntimeError(f'{prefix}: Newton step leaves the local branch (scaled size {mpmath.nstr(scale, 3)})')
        p = [p[i] + step[i] for i in range(3)]
    raise RuntimeError(f'{prefix}: Newton target {args.target} not reached')


def omega_matrix(mu):
    om = mpmath.zeros(12, 12)
    for i, m in enumerate([1, 1, mpf(mu), mpf(mu), 1, 1]):
        om[i, 6 + i] = m; om[6 + i, i] = -m
    return om


def floquet(phi, mu):
    n = 12
    om = omega_matrix(mu)
    sym = phi.T * om * phi - om
    traces, power = [], phi
    for k in range(4):
        traces.append(sum(power[i, i] for i in range(n)))
        power = power * phi
    p1, p2, p3, p4 = [t - 8 for t in traces]
    e1 = p1
    e2 = (e1 * p1 - p2) / 2
    e3 = (e2 * p1 - e1 * p2 + p3) / 3
    e4 = (e3 * p1 - e2 * p2 + e1 * p3 - p4) / 4
    disc = 2 * p2 + 8 - p1 ** 2
    roots = mpmath.polyroots([e4, -e3, e2, -e1, 1], maxsteps=400, extraprec=4 * mp.prec, asc=True)
    sroot = mpmath.sqrt(mpc(disc))
    pair_sums = [(e1 + sroot) / 2, (e1 - sroot) / 2]
    eigenvalues, right = mpmath.eig(phi, right=True)
    order = sorted(range(n), key=lambda i: abs(eigenvalues[i] - 1))
    trivial, nontrivial = order[:8], order[8:]
    modes = []
    for i in sorted(nontrivial, key=lambda i: (float(mpmath.im(eigenvalues[i])), float(mpmath.re(eigenvalues[i])))):
        lam = eigenvalues[i]
        z = right[:, i]
        z = z / mpmath.norm(z)
        residual = mpmath.norm(phi * z - lam * z)
        form = (z.H * om * z)[0]
        kappa = mpmath.re(1j * form)  # convention: kappa = real(i conj(z)^T Omega z)
        on_circle = abs(abs(lam) - 1) < mpf(10) ** (-20)
        modes.append({'lambda': [s(mpmath.re(lam)), s(mpmath.im(lam))], 'modulus': s(abs(lam)),
                      'eigen_residual': s(residual, 8), 'krein_form': s(kappa, 12),
                      'krein_signature': (1 if kappa > 0 else -1) if on_circle else None})
    quartic_vs_eig = max(min(abs(r - eigenvalues[i]) for r in roots) for i in nontrivial)
    return {
        'D_from_traces': s(disc),
        'pair_sums_s': [[s(mpmath.re(x)), s(mpmath.im(x))] for x in pair_sums],
        'quartic_coefficients': [s(e1), s(e2), s(e3), s(e4)],
        'palindromic_defects': {'abs(e3-e1)': s(abs(e3 - e1), 8), 'abs(e4-1)': s(abs(e4 - 1), 8)},
        'quartic_roots': [[s(mpmath.re(r)), s(mpmath.im(r))] for r in roots],
        'spectral_radius_nontrivial': s(max(abs(eigenvalues[i]) for i in nontrivial)),
        'max_trivial_distance_from_1': s(max(abs(eigenvalues[i] - 1) for i in trivial), 8),
        'min_nontrivial_distance_from_1': s(min(abs(eigenvalues[i] - 1) for i in nontrivial), 8),
        'quartic_roots_vs_eig_max_difference': s(quartic_vs_eig, 8),
        'full_symplectic_defect': s(max(abs(sym[i, j]) for i in range(n) for j in range(n)), 8),
        'modes': modes,
        'method': 'power sums tr(Phi^k)-8 with Newton identities; mpmath.eig on the full 12x12 monodromy; '
                  'full-space Krein forms with Omega=[[0,M],[-M,0]]',
    }


def verify(mu, p, history, args, raw_dir, prefix):
    rows, info = run_gbs(mu, [s(x) for x in p], args.digits, args.verify_tol, args.verify_stages, 12,
                         raw_dir / f'{prefix}_verification.csv', args.stride)
    x0 = initial_state(mu, p[0], p[1])
    e0, l0, _, _, _ = invariants(x0, mu)
    drift = [mpf(0)] * 4
    nearest, t_near = mpf('inf'), mpf(0)
    for row in rows:
        values = parse(row, 13)
        e, l, com, mom, dist = invariants(values[1:13], mu)
        drift = [max(drift[0], abs(e - e0)), max(drift[1], abs(l)), max(drift[2], *map(abs, com)), max(drift[3], *map(abs, mom))]
        if min(dist) < nearest:
            nearest, t_near = min(dist), values[0]
    final = parse(rows[-1])
    residual = [a - b for a, b in zip(final[1:13], x0)]
    phi = mpmath.matrix(12, 12)
    for c in range(12):
        for r in range(12):
            phi[r, c] = final[13 + 12 * c + r]
    record = {'mu': str(mu), 'u': s(p[0]), 'v': s(p[1]), 'T': s(p[2]), 'E': s(e0),
              'S': s(p[2] * abs(e0) ** mpf(1.5)),
              'return_vector': {k: s(r, 12) for k, r in zip(NAMES, residual)},
              'max_closure': s(max(abs(r) for r in residual), 8),
              'sampled_max_energy_drift': s(drift[0], 6), 'sampled_max_abs_angular_momentum': s(drift[1], 6),
              'sampled_max_com': s(drift[2], 6), 'sampled_max_momentum': s(drift[3], 6),
              'sampled_min_distance_upper_estimate': s(nearest, 15), 'sampled_t_nearest': s(t_near, 15),
              'newton_history': history, 'verification_integration': info,
              'full_monodromy': [[s(phi[r, c]) for c in range(12)] for r in range(12)]}
    started = time.perf_counter()
    record['floquet'] = floquet(phi, mu)
    record['floquet']['runtime_seconds'] = time.perf_counter() - started
    print(f"{prefix}: closure {record['max_closure']} D {mpmath.nstr(mpf(record['floquet']['D_from_traces']), 12)} "
          f"rho {mpmath.nstr(mpf(record['floquet']['spectral_radius_nontrivial']), 15)}", flush=True)
    return record


def double_seed(mu):
    """Binary64 continuation from the originally supplied sweep (independent of all MP records)."""
    import numpy as np
    from threebody_expand import continue_to
    rows = list(csv.DictReader((ROOT / 'threebody_codex_project' / 'pancake_mass_sweep.csv').open()))
    start = min(rows, key=lambda r: abs(float(r['mu']) - float(mu)))
    p = np.array([float(start[k]) for k in ['u', 'v', 'T']])
    if float(start['mu']) != float(mu):
        p, _, _ = continue_to(float(start['mu']), float(mu), p, step=.0002)
    return [repr(float(x)) for x in p], f"binary64 continuation from archived sweep point mu={start['mu']}"


def job(task):
    """One mass: seed, independent Newton, stricter full-period verification and Floquet analysis.

    Runs in its own process: the mpmath precision context is process-global and some mpmath
    routines change it temporarily, so threads must not share it.
    """
    mu, args, seeds = task
    mp.dps = args.digits
    prefix = 'gbs_' + mu.replace('.', 'p') + f'_d{args.digits}'
    if mu in seeds:
        seed, source = seeds[mu], f'seed file {args.seed_file}'
    else:
        seed, source = double_seed(mu)
    started = time.perf_counter()
    p, history = newton(mu, seed, args, args.raw_dir, prefix)
    record = verify(mu, p, history, args, args.raw_dir, prefix)
    record['seed'] = {'values': seed, 'source': source}
    record['total_runtime_seconds'] = time.perf_counter() - started
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--masses', required=True, help='comma-separated decimal masses (exact decimal strings)')
    parser.add_argument('--seed-file', type=Path, help='JSON {mu: [u, v, T]} decimal strings; otherwise binary64 continuation')
    parser.add_argument('--digits', type=int, default=70)
    parser.add_argument('--stages', type=int, default=18)
    parser.add_argument('--verify-stages', type=int, default=22)
    parser.add_argument('--newton-tol', default='1e-55')
    parser.add_argument('--verify-tol', default='1e-60')
    parser.add_argument('--target', default='1e-48')
    parser.add_argument('--max-newton', type=int, default=10)
    parser.add_argument('--stride', type=int, default=1)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--output', type=Path, default=ROOT / 'results' / 'independent_gbs_orbits.json')
    parser.add_argument('--raw-dir', type=Path, default=ROOT / 'results' / 'gbs_raw')
    args = parser.parse_args()
    mp.dps = args.digits
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    seeds = json.loads(args.seed_file.read_text()) if args.seed_file else {}
    masses = [m.strip() for m in args.masses.split(',') if m.strip()]

    jobs = [(mu, args, seeds) for mu in masses]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        records = list(pool.map(job, jobs))
    result = {'status': 'independent numerical cross-check; not interval-certified',
              'equations': 'q_i_ddot = sum_{j!=i} m_j (q_j-q_i)/|q_j-q_i|^3, G=1, masses (1, mu, 1)',
              'section': 'q=(-1,0),(0,0),(1,0); v1=v3=(u,v); v2=-2(u,v)/mu; full labeled 12-component return',
              'state_order': NAMES, 'precision_digits': args.digits,
              'settings': {k: str(v) for k, v in vars(args).items()},
              'mpmath_version': mpmath.__version__, 'mpmath_backend': mpmath.libmp.BACKEND,
              'source_sha256': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                for name in ['independent_mp.py', 'tools/gbs_mpfr.cpp']},
              'orbits': records}
    args.output.write_text(json.dumps(result, indent=1))
    print('wrote', args.output)


if __name__ == '__main__':
    main()
