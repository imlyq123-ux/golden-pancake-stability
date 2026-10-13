"""Binary64 pseudo-arclength continuation of the Golden Pancake isola with a Floquet scan.

Unknowns z = (u, v, T, mu) of the Euler-section ansatz; equations: the full labeled
12-component return X(T) - X(0) = 0 plus the pseudo-arclength condition. The 12x4 shooting
Jacobian uses the 12x12 monodromy and an exact mass-sensitivity column dX(T)/dmu obtained
from the inhomogeneous variational equation (d/dt) S = A(t) S + dF/dmu,
S(0) = dX(0)/dmu, so the continuation can pass folds in mu (e.g. mu_max ~ 0.63460).

At every converged point the nontrivial quartet is taken from the explicit symplectic
quotient of dense_sweep_double.quotient; pair sums s = lambda + 1/lambda, D, spectral radius
and Krein signs are recorded. This is an exploratory binary64 scan: transitions it finds
must be verified with the multiprecision pipelines.
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from dense_sweep_double import analyse
from threebody_expand import f, jac, xinit

ROOT = Path(__file__).resolve().parent
W = np.array([1., 1., .1, 1.])          # scaling of (u, v, T, mu) in the arclength metric
PAIRS = [(0, 1), (0, 2), (1, 2)]


def rhs(t, y, m):
    x = y[:12]
    a = jac(x, m)
    q = x[:6].reshape(3, 2)
    g = np.zeros(12)
    d12, d32 = q[1] - q[0], q[1] - q[2]
    g[6:8] = d12 / np.linalg.norm(d12) ** 3      # d(a_1)/d(mu): body 2 has mass mu
    g[10:12] = d32 / np.linalg.norm(d32) ** 3    # d(a_3)/d(mu); a_2 does not depend on mu
    return np.concatenate((f(t, x, m), (a @ y[12:156].reshape(12, 12)).ravel(), a @ y[156:168] + g))


try:  # compiled right-hand side when Numba is installed (same equations as rhs above)
    from numba import njit

    @njit(cache=True)
    def rhs_compiled(t, y, m):
        x = y[:12]
        a = jac(x, m)
        out = np.empty(168)
        out[:12] = f(t, x, m)
        out[12:156] = (a @ np.ascontiguousarray(y[12:156]).reshape(12, 12)).ravel()
        g = np.zeros(12)
        d12 = x[2:4] - x[0:2]
        d32 = x[2:4] - x[4:6]
        g[6:8] = d12 / np.sqrt(d12 @ d12) ** 3
        g[10:12] = d32 / np.sqrt(d32 @ d32) ** 3
        out[156:168] = a @ np.ascontiguousarray(y[156:168]) + g
        return out
    RHS = rhs_compiled
except ImportError:
    RHS = rhs


def shoot(z, rtol, atol):
    u, v, T, mu = z
    m = np.array([1., mu, 1.])
    x0 = xinit(mu, u, v)
    s0 = np.zeros(12); s0[8], s0[9] = 2 * u / mu ** 2, 2 * v / mu ** 2
    sol = solve_ivp(RHS, (0, T), np.concatenate((x0, np.eye(12).ravel(), s0)), method='DOP853',
                    rtol=rtol, atol=atol, args=(m,))
    if not sol.success:
        raise RuntimeError(sol.message)
    xt, phi, st = sol.y[:12, -1], sol.y[12:156, -1].reshape(12, 12), sol.y[156:168, -1]
    du = np.zeros(12); du[[6, 8, 10]] = [1, -2 / mu, 1]
    dv = np.zeros(12); dv[[7, 9, 11]] = [1, -2 / mu, 1]
    jacobian = np.column_stack((phi @ du - du, phi @ dv - dv, f(0, xt, m), st - s0))
    pos = sol.y[:6]
    rmin = min(float(np.min(np.linalg.norm(pos[2 * j:2 * j + 2] - pos[2 * i:2 * i + 2], axis=0))) for i, j in PAIRS)
    return xt - x0, jacobian, phi, rmin, sol


def tangent(jacobian, previous=None):
    _, sv, vt = np.linalg.svd(jacobian / W)       # columns scaled: J W^-1
    tau = vt[-1]
    if previous is not None and tau @ previous < 0:
        tau = -tau
    return tau, sv


def correct(z0, tau, h, rtol, atol, tol=2e-11, max_iter=7):  # tol: max |X(T)-X(0)| accepted
    z = z0 + h * tau / W
    for it in range(max_iter):
        r, jacobian, phi, rmin, sol = shoot(z, rtol, atol)
        n = tau @ (W * (z - z0)) - h
        if np.max(abs(r)) < tol and abs(n) < 1e-12:
            return z, r, jacobian, phi, rmin, it
        system = np.vstack((jacobian, tau * W))
        step = np.linalg.lstsq(system, -np.concatenate((r, [n])), rcond=None)[0]
        z = z + step
    raise RuntimeError(f'corrector did not converge (|r|={np.max(abs(r)):.2e})')


def run(seed, direction, args):
    z = np.array(seed, dtype=float)
    # converge the seed at fixed mu first (Gauss-Newton on (u, v, T))
    for _ in range(8):
        r, jacobian, phi, rmin, sol = shoot(z, args.rtol, args.atol)
        if np.max(abs(r)) < args.newton_tol:
            break
        z[:3] += np.linalg.lstsq(jacobian[:, :3], -r, rcond=None)[0]
    tau, _ = tangent(jacobian)
    if np.sign(tau[3]) != np.sign(direction):
        tau = -tau
    h, s_arc, out = args.h0, 0.0, []
    started = time.perf_counter()
    for k in range(args.max_steps):
        rec = analyse(float(z[3]), z[:3], phi)
        rec.update(arclength=s_arc, sampled_rmin=rmin, closure=float(np.max(abs(r))), step_h=h,
                   tangent=tau.tolist(), elapsed=time.perf_counter() - started)
        out.append(rec)
        print(f"[{args.label}] s={s_arc:.5f} mu={z[3]:.7f} u={z[0]:.6f} T={z[2]:.5f} rmin={rmin:.2e} "
              f"D={rec['D']:+.5e} rho={rec['spectral_radius']:.6f} s=({rec['s_values'][0][0]:+.4f},{rec['s_values'][1][0]:+.4f})",
              flush=True)
        if rmin < args.rmin_stop or not (args.mu_min < z[3] < args.mu_max):
            break
        while True:
            try:
                znew, r, jacobian, phi, rmin, its = correct(z, tau, h, args.rtol, args.atol, args.newton_tol)
                if abs(znew[3] - z[3]) > args.max_dmu:
                    raise RuntimeError('mass step too large')
                break
            except RuntimeError as exc:
                h *= .5
                if h < 1e-7:
                    out.append({'stop': f'step size underflow: {exc}', 'mu': float(z[3])})
                    return out
        s_arc += h
        tau, _ = tangent(jacobian, tau)
        z = znew
        h = min(args.h_max, h * (1.4 if its <= 2 else (.7 if its >= 4 else 1.0)))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--seed', required=True, help='u,v,T,mu (binary64 is enough; the seed is corrected first)')
    parser.add_argument('--direction', type=float, required=True, help='+1: initially increasing mu, -1: decreasing')
    parser.add_argument('--label', default='scan')
    parser.add_argument('--h0', type=float, default=.002)
    parser.add_argument('--h-max', type=float, default=.004)
    parser.add_argument('--max-dmu', type=float, default=.0006)
    parser.add_argument('--max-steps', type=int, default=600)
    parser.add_argument('--rmin-stop', type=float, default=3e-4)
    parser.add_argument('--mu-min', type=float, default=.40)
    parser.add_argument('--mu-max', type=float, default=.70)
    parser.add_argument('--newton-tol', type=float, default=2e-11)
    parser.add_argument('--rtol', type=float, default=1e-12)
    parser.add_argument('--atol', type=float, default=1e-14)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    seed = [float(x) for x in args.seed.split(',')]
    records = run(seed, args.direction, args)
    args.output.write_text(json.dumps({'status': 'binary64 exploratory pseudo-arclength scan', 'seed': seed,
                                       'direction': args.direction, 'settings': {k: str(v) for k, v in vars(args).items()},
                                       'records': records}, indent=1))


if __name__ == '__main__':
    main()
