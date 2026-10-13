"""Classify the two spectral transitions from orbits at mu* +/- 10^-k (GBS + mpmath pipeline).

For a non-semisimple (Jordan-block) collision of two unit-circle multipliers with opposite
Krein signatures -- the linear part of a Hamiltonian-Hopf bifurcation -- one expects
  stable side:   |lambda_1 - lambda_2| ~ delta^(1/2), Krein forms of normalised eigenvectors -> 0
                 like delta^(1/2), and 1 - |<z_1, z_2>| ~ delta (eigenvectors coalesce);
  unstable side: rho - 1 ~ delta^(1/2) (complex quartet leaves the circle);
  at mu*:        sigma_min(Phi - lambda_bar I) ~ delta while sigma_2 = O(1) (geometric multiplicity 1).
A semisimple collision would instead give linear splitting and sigma_min ~ sigma_2.
"""
import argparse
import json
import time
from argparse import Namespace
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import mpmath
from mpmath import mp, mpf

import independent_mp as im

ROOT = Path(__file__).resolve().parent
EXPONENTS = [4, 6, 8, 10, 12, 15, 20, 25, 30, 35, 40]


def seed_from(evaluations, mu):
    near = sorted(evaluations, key=lambda r: abs(mpf(r['mu']) - mpf(mu)))[:2]
    m0, m1 = mpf(near[0]['mu']), mpf(near[1]['mu'])
    a = (mpf(mu) - m0) / (m1 - m0)
    return [im.s((1 - a) * mpf(near[0][k]) + a * mpf(near[1][k])) for k in ['u', 'v', 'T']]


def colliding_pair(phi, lam_c):
    eigenvalues, right = mpmath.eig(phi, right=True)
    idx = sorted(range(12), key=lambda i: abs(eigenvalues[i] - lam_c))[:2]
    return [(eigenvalues[i], right[:, i] / mpmath.norm(right[:, i])) for i in idx], eigenvalues


def job(task):
    name, mu, side, delta_exp, seed, lam_c, args = task
    mp.dps = args.digits
    gargs = Namespace(digits=args.digits, stages=20, verify_stages=24, newton_tol='1e-62', verify_tol='1e-66',
                      target='1e-55', max_newton=10, stride=100)
    prefix = f"ta_{name}_{side}{delta_exp}"
    p, history = im.newton(mu, seed, gargs, args.raw_dir, prefix)
    rec = im.verify(mu, p, history, gargs, args.raw_dir, prefix)
    phi = mpmath.matrix(12, 12)
    for r in range(12):
        for c in range(12):
            phi[r, c] = mpf(rec['full_monodromy'][r][c])
    lam_c = mpmath.mpc(*map(mpf, lam_c))
    pair, eigenvalues = colliding_pair(phi, lam_c)
    om = im.omega_matrix(mu)
    (l1, z1), (l2, z2) = pair
    kappas = [mpmath.re(1j * (z.H * om * z)[0]) for z in (z1, z2)]
    overlap = abs((z1.H * z2)[0])
    out = {'boundary': name, 'mu': mu, 'side': side, 'delta_exponent': delta_exp,
           'D': rec['floquet']['D_from_traces'], 'closure': rec['max_closure'],
           'lambda_1': [im.s(l1.real, 40), im.s(l1.imag, 40)], 'lambda_2': [im.s(l2.real, 40), im.s(l2.imag, 40)],
           'split_abs': im.s(abs(l1 - l2), 12), 'moduli_minus_1': [im.s(abs(l1) - 1, 12), im.s(abs(l2) - 1, 12)],
           'krein_forms_normalised': [im.s(k, 12) for k in kappas], 'eigenvector_overlap': im.s(overlap, 30),
           'one_minus_overlap': im.s(1 - overlap, 12)}
    if delta_exp is None:  # at mu*: Jordan-structure test
        lam_bar = (l1 + l2) / 2
        sv = mpmath.svd_c(phi - lam_bar * mpmath.eye(12), compute_uv=False)
        sv = sorted(mpf(abs(sv[i])) for i in range(sv.rows))
        out['singular_values_Phi_minus_lambda_bar'] = [im.s(x, 10) for x in sv]
        out['lambda_bar'] = [im.s(lam_bar.real, 40), im.s(lam_bar.imag, 40)]
    print(f"{prefix}: split={mpmath.nstr(abs(l1 - l2), 5)} rho-1={mpmath.nstr(max(abs(l1), abs(l2)) - 1, 5)} "
          f"kappa={[mpmath.nstr(k, 4) for k in kappas]}", flush=True)
    return out


def slope(points):
    xs = [mpmath.log10(mpf(10) ** (-k)) for k, _ in points]
    ys = [mpmath.log10(v) for _, v in points]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    return float(sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--refinement-dir', type=Path, default=ROOT / 'results' / 'boundary_refinement')
    parser.add_argument('--pipeline', default='gbs80')
    parser.add_argument('--digits', type=int, default=80)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--output', type=Path, default=ROOT / 'results' / 'transition_analysis.json')
    parser.add_argument('--raw-dir', type=Path, default=ROOT / 'results' / 'gbs_raw_transition')
    args = parser.parse_args()
    mp.dps = args.digits
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    tasks, meta = [], {}
    for name in ['lower', 'upper']:
        ref = json.loads((args.refinement_dir / f'{name}_{args.pipeline}.json').read_text())
        star = mpf(ref['best_estimate'])
        evals = ref['evaluations']
        at_star = next(r for r in evals if r['mu'] == ref['best_estimate'])
        upper_half = [m for m in at_star['multipliers'] if mpf(m[1]) > 0]
        lam_c = [im.s(sum(mpf(m[0]) for m in upper_half) / 2), im.s(sum(mpf(m[1]) for m in upper_half) / 2)]
        # stable side is where D > 0
        d_lo, d_hi = map(mpf, ref['bracket_D'])
        stable_sign = 1 if d_hi > 0 else -1
        meta[name] = {'mu_star': ref['best_estimate'], 'lambda_c_estimate': lam_c, 'stable_direction': stable_sign}
        tasks.append((name, ref['best_estimate'], 'star', None, [at_star[k] for k in ['u', 'v', 'T']], lam_c, args))
        for k in EXPONENTS:
            for side, sign in [('stable', stable_sign), ('unstable', -stable_sign)]:
                mu = mpmath.nstr(star + sign * mpf(10) ** (-k), args.digits - 10, strip_zeros=False)
                tasks.append((name, mu, side, k, seed_from(evals, mu), lam_c, args))
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(job, tasks))
    summary = {}
    for name in ['lower', 'upper']:
        st = [(r['delta_exponent'], mpf(r['split_abs'])) for r in rows if r['boundary'] == name and r['side'] == 'stable']
        un = [(r['delta_exponent'], max(map(mpf, r['moduli_minus_1']))) for r in rows if r['boundary'] == name and r['side'] == 'unstable']
        kap = [(r['delta_exponent'], abs(mpf(r['krein_forms_normalised'][0]))) for r in rows if r['boundary'] == name and r['side'] == 'stable']
        ov = [(r['delta_exponent'], mpf(r['one_minus_overlap'])) for r in rows if r['boundary'] == name and r['side'] == 'stable']
        star_row = next(r for r in rows if r['boundary'] == name and r['side'] == 'star')
        signs = [[1 if mpf(k) > 0 else -1 for k in r['krein_forms_normalised']] for r in rows if r['boundary'] == name and r['side'] == 'stable']
        summary[name] = {**meta[name],
                         'loglog_slope_stable_split': slope(st), 'loglog_slope_unstable_rho_minus_1': slope(un),
                         'loglog_slope_krein_form': slope(kap), 'loglog_slope_one_minus_overlap': slope(ov),
                         'stable_side_krein_sign_pairs': signs,
                         'singular_values_at_mu_star': star_row['singular_values_Phi_minus_lambda_bar'],
                         'lambda_bar_at_mu_star': star_row['lambda_bar']}
    out = {'description': __doc__, 'pipeline': 'GBS + mpmath (independent_mp.py)', 'digits': args.digits,
           'runtime_seconds': time.perf_counter() - started, 'summary': summary,
           'rows': sorted(rows, key=lambda r: (r['boundary'], r['side'], r['delta_exponent'] or 0))}
    args.output.write_text(json.dumps(out, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
