"""Phase E (exploratory, binary64): stability transitions for unequal outer masses (1+delta, mu, 1-delta).

Shooting set-up (derived for general masses; the symmetric Euler ansatz v1 = v3 is NOT used):
  section   collinear configuration on the x-axis at t = 0 (fixes rotation phase and time phase):
            q1 = (x1, 0), q2 = (a, 0), q3 = (x1 + 2, 0) with x1 = -(mu a + 2 m3)/(m1 + m3)
            (centre of mass at the origin, scale |q3 - q1| = 2);
  velocities v1 = (u1, w1), v3 = (u3, w3) with w3 = -m1 w1 (x1 - a) / (m3 (x3 - a)) (L = 0) and
            v2 = -(m1 v1 + m3 v3)/mu (zero momentum);
  unknowns  (a, u1, u3, w1, T, theta);
  equations X(T) - R_theta X(0) = 0, all 12 labelled components (relative-periodic form; 12
            equations of rank 6 because of the COM, momentum, energy and L identities).
For delta = 0 this reduces to the Euler ansatz (a = 0, u1 = u3 = u, w1 = v, theta = 0). Reflections
commute with the flow for any masses, so the half-period reflection symmetry X(T/2) = sigma X(0)
and with it theta = 0 are expected to persist; theta is nevertheless solved for and reported.

Floquet analysis of M = R_{-theta} Phi(T) on the explicit quotient (COM, momentum, dH = dL = 0,
time and rotation gauges) with the mass-weighted symplectic form. Transition curves mu_*(delta)
are followed by secant iteration on D(mu, delta) = 2 tr(A^2) + 8 - tr(A)^2 at each delta.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from threebody_expand import ext_numba, f

ROOT = Path(__file__).resolve().parent
PAIRS = [(0, 1), (0, 2), (1, 2)]


def masses(mu, delta):
    return np.array([1. + delta, mu, 1. - delta])


def initial_state(p, mu, delta):
    """X(0) from p = (a, u1, u3, w1); works for complex p (complex-step derivatives)."""
    a, u1, u3, w1 = p
    m1, m2, m3 = 1. + delta, mu, 1. - delta
    x1 = -(m2 * a + 2 * m3) / (m1 + m3)
    x3 = x1 + 2
    w3 = -m1 * w1 * (x1 - a) / (m3 * (x3 - a))
    v2x = -(m1 * u1 + m3 * u3) / m2
    v2y = -(m1 * w1 + m3 * w3) / m2
    zero = 0 * a
    return np.array([x1, zero, a, zero, x3, zero, u1, w1, v2x, v2y, u3, w3])


def d_initial_state(p, mu, delta, h=1e-30):
    cols = []
    for k in range(4):
        z = np.array(p, dtype=complex)
        z[k] += 1j * h
        cols.append(initial_state(z, mu, delta).imag / h)
    return np.column_stack(cols)


def rot(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.kron(np.eye(6), np.array([[c, -s], [s, c]]))


JROT = np.kron(np.eye(6), np.array([[0., -1.], [1., 0.]]))


def shoot(z, mu, delta, rtol=3e-13, atol=1e-15):
    p, T, theta = z[:4], z[4], z[5]
    m = masses(mu, delta)
    x0 = initial_state(p, mu, delta)
    sol = solve_ivp(lambda t, y: ext_numba(t, y, m), (0, T), np.concatenate((x0, np.eye(12).ravel())),
                    method='DOP853', rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    xt, phi = sol.y[:12, -1], sol.y[12:, -1].reshape(12, 12)
    r_t = rot(theta)
    dx0 = d_initial_state(p, mu, delta)
    jac = np.column_stack((phi @ dx0 - r_t @ dx0, f(0, xt, m), -(JROT @ r_t @ x0)))
    pos = sol.y[:6]
    rmin = min(float(np.min(np.linalg.norm(pos[2 * j:2 * j + 2] - pos[2 * i:2 * i + 2], axis=0))) for i, j in PAIRS)
    return xt - r_t @ x0, jac, phi, x0, rmin


def newton(z, mu, delta, tol=2e-11, max_iter=12):
    z = np.array(z, dtype=float)
    for it in range(max_iter):
        r, jac, phi, x0, rmin = shoot(z, mu, delta)
        if np.max(abs(r)) < tol:
            return z, r, phi, x0, rmin, it
        step = np.linalg.lstsq(jac, -r, rcond=None)[0]
        scale = max(1., np.max(abs(step) / np.array([.01, .003, .003, .003, .05, .01])))
        z = z + step / scale
    raise RuntimeError(f'Newton failed at mu={mu}, delta={delta}: |r|={np.max(abs(r)):.2e}')


def quotient_D(phi, theta, x0, mu, delta):
    m = masses(mu, delta)
    mm = np.diag(np.repeat(m, 2))
    zero = np.zeros((6, 6))
    omega = np.block([[zero, mm], [-mm, zero]])
    M = rot(-theta) @ phi
    flow = f(0, x0, m)
    rgen = JROT @ x0
    rows = []
    for start in (0, 6):
        for k in range(2):
            r = np.zeros(12); r[start + k:start + 6:2] = m; rows.append(r)
    rows += [flow @ omega, rgen @ omega, flow, rgen]
    q = np.linalg.svd(np.array(rows))[2][8:].T
    n = np.column_stack((flow, rgen))
    proj = np.eye(12) - n @ np.linalg.solve(n.T @ n, n.T)
    a = q.T @ proj @ M @ q
    lam = np.linalg.eigvals(a)
    tr1, tr2 = np.trace(a), np.trace(a @ a)
    return float(2 * tr2 + 8 - tr1 ** 2), lam, float(np.max(abs(proj @ M @ q - q @ a)))


def evaluate(z, mu, delta):
    z, r, phi, x0, rmin, its = newton(z, mu, delta)
    D, lam, rec = quotient_D(phi, z[5], x0, mu, delta)
    return z, D, lam, {'closure': float(np.max(abs(r))), 'rmin': rmin, 'newton_iterations': its, 'reconstruction_defect': rec}


def continue_delta(mu, uvT, deltas):
    """Sequential continuation in delta at fixed mu with a linear predictor.

    Returns (solutions by delta, failure or None). A failure stops the continuation: later deltas
    have no predecessor and are reported as not scanned, never silently skipped."""
    z = np.array([0., uvT[0], uvT[0], uvT[1], uvT[2], 0.])
    out, prev = {}, None
    for delta in deltas:
        guess = z if prev is None else 2 * z - prev
        try:
            znew, D, lam, info = evaluate(guess, mu, delta)
        except RuntimeError as exc:
            return out, {'delta': delta, 'failure': str(exc)}
        prev, z = z, znew
        out[delta] = (z.copy(), D, info)
        print(f"[delta-continuation] mu={mu} delta={delta:.5f} D={D:+.6f} theta={z[5]:+.1e} a={z[0]:+.4e}", flush=True)
    return out, None


def grid_masses(mu_center, mu_lo, mu_hi, step):
    """The mu grid of scan_mu: mu_center + k step within [mu_lo, mu_hi]."""
    out = []
    for direction in (-1, 1):
        k = 0 if direction == -1 else 1
        while True:
            mu = round(mu_center + direction * k * step, 10)
            if mu < mu_lo - 1e-12 or mu > mu_hi + 1e-12:
                break
            out.append(mu)
            k += 1
    return sorted(out)


def scan_mu(task):
    """mu-continuation at fixed delta over a grid; refine every D sign change by secant.

    Newton failures are kept (`failures`) and reported in `coverage`; a row is `complete` only if
    every grid point was computed and every crossing was refined. Incomplete rows must not be
    read as physical instability (missing or negative-D points can be solver failures)."""
    delta, z_center, mu_center, mu_lo, mu_hi, step = task
    expected = grid_masses(mu_center, mu_lo, mu_hi, step)
    points, failures = [], []
    for direction in (-1, 1):
        z, prev = z_center.copy(), None
        masses_dir = [m for m in expected if (m <= mu_center + 1e-12 if direction == -1 else m > mu_center + 1e-12)]
        for mu in (masses_dir[::-1] if direction == -1 else masses_dir):
            guess = z if prev is None else 2 * z - prev
            try:
                znew, D, lam, info = evaluate(guess, mu, delta)
            except RuntimeError as exc:
                failures.append({'stage': 'grid', 'mu': mu, 'failure': str(exc),
                                 'consequence': 'continuation in this direction stopped'})
                break
            prev, z = z, znew
            points.append({'mu': mu, 'D': D, 'rho': float(max(abs(lam))), 'theta': z[5], 'z': z.tolist(), **info,
                           'stable': bool(D > 0 and max(abs(lam)) < 1 + 1e-6)})
    points.sort(key=lambda g: g['mu'])
    index = {m: i for i, m in enumerate(expected)}
    crossings = []
    for g0, g1 in zip(points, points[1:]):
        if index[g1['mu']] - index[g0['mu']] != 1:      # never bracket across a gap in the grid
            continue
        if (g0['D'] > 0) == (g1['D'] > 0):
            continue
        pts = [(g0['mu'], g0['D'], np.array(g0['z'])), (g1['mu'], g1['D'], np.array(g1['z']))]
        record = {'grid_bracket': [g0['mu'], g1['mu']], 'stable_side': 'above' if g1['D'] > 0 else 'below'}
        try:
            for _ in range(40):
                (m0, d0, _), (m1, d1, z1) = pts[-2], pts[-1]
                mu_new = m1 - d1 * (m1 - m0) / (d1 - d0)
                znew, D, lam, info = evaluate(z1, mu_new, delta)
                pts.append((mu_new, D, znew))
                if abs(mu_new - m1) < 1e-13 or abs(D) < 1e-12:
                    break
            upper = sorted([w for w in lam if w.imag > 0], key=lambda w: w.real)
            record.update(refined=True, mu_star=pts[-1][0], D_at_mu_star=pts[-1][1], theta=float(pts[-1][2][5]),
                          collision_multiplier=[float(np.mean([w.real for w in upper])), float(np.mean([w.imag for w in upper]))])
        except RuntimeError as exc:
            record.update(refined=False, failure=str(exc))
            failures.append({'stage': 'crossing', 'grid_bracket': record['grid_bracket'], 'failure': str(exc)})
        crossings.append(record)
    computed = [g['mu'] for g in points]
    coverage = {'requested_range': [mu_lo, mu_hi], 'expected_points': len(expected), 'computed_points': len(points),
                'computed_range': [min(computed), max(computed)] if computed else None,
                'missing_points': [m for m in expected if m not in set(computed)],
                'complete': not failures and len(points) == len(expected)}
    status = 'complete' if coverage['complete'] else f"INCOMPLETE ({len(failures)} failures, {len(coverage['missing_points'])} missing)"
    max_d = f"{max(g['D'] for g in points):+.4f}" if points else 'n/a'
    print(f"delta={delta:.5f}: {len(points)}/{len(expected)} grid points ({status}), crossings at "
          f"{[round(c['mu_star'], 10) for c in crossings if c.get('refined')]}, max D={max_d}", flush=True)
    return {'delta': delta, 'grid': points, 'failures': failures, 'coverage': coverage, 'crossings': crossings}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--label', default='main_window')
    parser.add_argument('--mu-center', type=float, default=0.6125)
    parser.add_argument('--uvT', default='0.15175005652,0.13970811655,8.8457011',
                        help='Euler-section u,v,T of the symmetric orbit at mu-center (binary64 seed)')
    parser.add_argument('--mu-range', default='0.594,0.624')
    parser.add_argument('--mu-step', type=float, default=0.0005)
    parser.add_argument('--deltas', default='0:0.012:0.0005', help='start:stop:step (inclusive)')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'results' / 'asymmetric')
    args = parser.parse_args()
    from concurrent.futures import ProcessPoolExecutor
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    a, b, h = map(float, args.deltas.split(':'))
    deltas = [round(a + k * h, 10) for k in range(int(round((b - a) / h)) + 1)]
    # correct the binary64 seed at the centre (delta = 0) before continuing
    uvT = [float(x) for x in args.uvT.split(',')]
    centre, centre_failure = continue_delta(args.mu_center, uvT, deltas)
    lo, hi = map(float, args.mu_range.split(','))
    tasks = [(d, centre[d][0], args.mu_center, lo, hi, args.mu_step) for d in deltas if d in centre]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(scan_mu, tasks))
    # symmetry check: delta -> -delta must give the same transition masses
    mirror = None
    d_check = deltas[len(deltas) // 3]
    if d_check > 0 and d_check in centre:
        zc = centre[d_check][0]
        # rotate by pi, relabel 1<->3 and reverse time: (a, u1, u3, w1) -> (-a, u3, u1, w3) at -delta
        w3 = initial_state(zc[:4], args.mu_center, d_check)[11]
        zm = np.array([-zc[0], zc[2], zc[1], w3, zc[4], -zc[5]])
        try:
            mirror = scan_mu((-d_check, evaluate(zm, args.mu_center, -d_check)[0], args.mu_center, lo, hi, args.mu_step))
        except RuntimeError as exc:
            mirror = {'delta': -d_check, 'failure': str(exc), 'coverage': {'complete': False}}
    (args.output_dir / f'{args.label}.json').write_text(json.dumps(
        {'status': 'binary64 exploratory 2D scan; not multiprecision-verified', 'label': args.label,
         'settings': {k: str(v) for k, v in vars(args).items()}, 'runtime_seconds': time.perf_counter() - started,
         'centre_continuation': {str(d): {'z': v[0].tolist(), 'D': v[1], **v[2]} for d, v in centre.items()},
         'centre_continuation_failure': centre_failure,
         'deltas_requested': deltas, 'deltas_scanned': [r['delta'] for r in rows],
         'all_rows_complete': centre_failure is None and all(r['coverage']['complete'] for r in rows),
         'rows': rows, 'mirror_check': mirror}, indent=1, default=str))


if __name__ == '__main__':
    main()
