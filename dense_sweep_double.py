"""Dense binary64 continuation of the branch from the single mu=0.613 seed.

Unlike make_pancake_results.py, no archived per-mass seed is used: each mass is reached by
predictor-corrector continuation from its neighbour, so the archived 23-point sweep is
re-derived independently. At every mass the nontrivial quartet is obtained from an explicit
symplectic quotient of the full monodromy (COM/momentum, dH = dL = 0, time and rotation
gauges) and compared with the legacy "4 multipliers farthest from 1" heuristic. Pair sums
s = lambda + 1/lambda locate every possible spectral change: Krein collisions (s1 = s2),
period doubling (s = -2) and tangent bifurcations (s = +2).
"""
import argparse
import csv
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from threebody_expand import correction, energy, f, xinit

ROOT = Path(__file__).resolve().parent
SEED = (0.15148154107900583396835825988003861711916800727,
        0.14030301957068408729274415253719345009218867155,
        8.846608333332378559519969249427806733522940068)


def quotient(phi, x0, mu):
    """Reduced 4x4 monodromy and symplectic form on the explicit quotient (binary64)."""
    m = np.array([1., mu, 1.])
    mass = np.diag(np.repeat(m, 2))
    zero = np.zeros((6, 6))
    omega = np.block([[zero, mass], [-mass, zero]])
    flow = f(0, x0, m)
    rot = np.concatenate([np.column_stack((-x0[1:6:2], x0[0:6:2])).ravel(),
                          np.column_stack((-x0[7:12:2], x0[6:12:2])).ravel()])
    rows = []
    for start in (0, 6):
        for k in range(2):
            r = np.zeros(12); r[start + k:start + 6:2] = m; rows.append(r)
    rows += [flow @ omega, rot @ omega, flow, rot]
    constraints = np.array(rows)
    _, sv, vt = np.linalg.svd(constraints)
    q = vt[8:].T                                   # orthonormal basis of the 4-dim quotient chart
    n = np.column_stack((flow, rot))
    proj = np.eye(12) - n @ np.linalg.solve(n.T @ n, n.T)
    a = q.T @ proj @ phi @ q
    om_r = q.T @ omega @ q
    return a, om_r, {'constraint_singular_values': sv.tolist(),
                     'reconstruction_defect': float(np.max(abs(proj @ phi @ q - q @ a))),
                     'reduced_symplectic_defect': float(np.max(abs(a.T @ om_r @ a - om_r))),
                     'full_symplectic_defect': float(np.max(abs(phi.T @ omega @ phi - omega)))}


def analyse(mu, p, phi):
    x0 = xinit(mu, *p[:2])
    a, om_r, checks = quotient(phi, x0, mu)
    lam, vec = np.linalg.eig(a)
    modes = []
    for k in np.argsort(-lam.imag):
        z = vec[:, k] / np.linalg.norm(vec[:, k])
        kappa = float(np.real(1j * np.conj(z) @ om_r @ z))
        on_circle = abs(abs(lam[k]) - 1) < 1e-7
        modes.append({'lambda': [float(lam[k].real), float(lam[k].imag)], 'modulus': float(abs(lam[k])),
                      'krein_form': kappa, 'krein_signature': (int(np.sign(kappa)) if on_circle else None)})
    tr1, tr2 = np.trace(a), np.trace(a @ a)
    disc = float(2 * tr2 + 8 - tr1 ** 2)
    root = np.sqrt(complex(disc))
    s_values = [(tr1 + root) / 2, (tr1 - root) / 2]
    full = np.linalg.eigvals(phi)
    heuristic = full[np.argsort(abs(full - 1))[::-1][:4]]
    mismatch = max(min(abs(h - lam)) for h in heuristic)
    trivial = full[np.argsort(abs(full - 1))[:8]]
    return {'mu': mu, 'u': p[0], 'v': p[1], 'T': p[2], 'E': energy(mu, x0),
            'S': p[2] * abs(energy(mu, x0)) ** 1.5, 'D': disc,
            's_values': [[z.real, z.imag] for z in s_values],
            'spectral_radius': float(max(abs(lam))), 'modes': modes,
            'heuristic_vs_quotient_max_difference': float(mismatch),
            'heuristic_spectral_radius': float(max(abs(heuristic))),
            'trivial_max_distance_from_1': float(max(abs(trivial - 1))),
            'nontrivial_min_distance_from_1': float(min(abs(lam - 1))),
            'nontrivial_min_distance_from_minus_1': float(min(abs(lam + 1))),
            'checks': checks}


def chain(task):
    start, stop, step = task
    count = int(round(abs(stop - start) / step))
    sign = 1 if stop > start else -1
    out, history = [], []
    p_prev = None
    p = np.array(SEED)
    for k in range(0, count + 1):
        mu = round(start + sign * k * step, 10)
        guess = p if p_prev is None else 2 * p - p_prev    # secant predictor in mu
        started = time.perf_counter()
        try:
            p_new, residual, sol = correction(mu, guess)
        except RuntimeError as exc:
            out.append({'mu': mu, 'failure': str(exc)})
            break
        phi = sol.y[12:, -1].reshape(12, 12)
        rec = analyse(mu, p_new, phi)
        rec.update(closure=residual, newton_iterations=len([h for h in sol.correction_info['history'] if 'iteration' in h]),
                   runtime_seconds=time.perf_counter() - started)
        out.append(rec)
        print(f"mu={mu:.4f} D={rec['D']:+.6e} rho={rec['spectral_radius']:.9f} closure={residual:.1e} "
              f"heur-diff={rec['heuristic_vs_quotient_max_difference']:.1e}", flush=True)
        p_prev, p = (p if k > 0 else None), p_new
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--low', type=float, default=0.580)
    parser.add_argument('--high', type=float, default=0.6345)
    parser.add_argument('--step', type=float, default=0.0005)
    parser.add_argument('--output', type=Path, default=ROOT / 'results' / 'dense_sweep_double.json')
    args = parser.parse_args()
    with ProcessPoolExecutor(max_workers=2) as pool:
        down, up = pool.map(chain, [(0.613, args.low, args.step), (0.613, args.high, args.step)])
    records = sorted([r for r in down + up[1:]], key=lambda r: r['mu'])
    archived = {float(r['mu']): r for r in csv.DictReader((ROOT / 'threebody_codex_project' / 'pancake_mass_sweep.csv').open())}
    comparison = []
    for r in records:
        if 'failure' in r:
            continue
        match = next((a for m, a in archived.items() if abs(m - r['mu']) < 1e-12), None)
        if match:
            comparison.append({'mu': r['mu'], **{k: r[k] - float(match[k]) for k in ['u', 'v', 'T', 'E', 'D']},
                               'rho_difference': r['spectral_radius'] - float(match['max_lambda'])})
    result = {'status': 'binary64 exploratory scan; quartet from explicit symplectic quotient',
              'seed': list(SEED), 'step': args.step, 'records': records,
              'comparison_with_archived_23_point_sweep': comparison}
    args.output.write_text(json.dumps(result, indent=1))
    print('archived masses matched:', len(comparison))
    if comparison:
        for key in ['u', 'v', 'T', 'D', 'rho_difference']:
            print(key, 'max |difference| =', max(abs(c[key]) for c in comparison))


if __name__ == '__main__':
    main()
