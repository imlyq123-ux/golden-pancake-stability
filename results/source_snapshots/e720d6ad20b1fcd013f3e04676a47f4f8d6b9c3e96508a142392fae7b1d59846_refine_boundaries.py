"""Refine the two candidate spectral transitions by secant iteration on D(mu).

Two pipelines that share no integration, shooting or Floquet code:
  taylor  tools/taylor_mpfr + high_precision.correct/verify_orbit (Decimal quotient reduction)
  gbs     tools/gbs_mpfr + independent_mp.newton/verify (mpmath, reduction-free trace formula)

Every D(mu) evaluation is a fresh Newton-corrected periodic orbit (all 12 return components)
followed by a separate stricter full-period monodromy integration. After the secant
iterates converge, D is evaluated at mu* -/+ half_width to verify a sign bracket.
The bracket is an empirically verified numerical bracket, not an interval-arithmetic proof.
"""
import argparse
import json
import time
from argparse import Namespace
from decimal import Decimal, getcontext
from pathlib import Path

import mpmath
from mpmath import mp, mpf

import high_precision as hp
import independent_mp as im

ROOT = Path(__file__).resolve().parent
BRACKETS = {'lower': ('0.59904378', '0.59904379'), 'upper': ('0.61784426', '0.61784428')}
SETTINGS = {  # pipeline label -> numerical settings
    'taylor50': dict(backend='taylor', digits=50, order=32, tol='1e-40', vorder=36, vtol='1e-44', target='1e-36'),
    'taylor70': dict(backend='taylor', digits=70, order=44, tol='1e-55', vorder=48, vtol='1e-60', target='1e-50'),
    'taylor90': dict(backend='taylor', digits=90, order=56, tol='1e-70', vorder=60, vtol='1e-75', target='1e-64'),
    'gbs80': dict(backend='gbs', digits=80, stages=20, tol='1e-62', vstages=24, vtol='1e-66', target='1e-55'),
}


class Evaluator:
    def __init__(self, label, raw_dir):
        self.label, self.cfg, self.raw_dir = label, SETTINGS[label], raw_dir
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.orbits = {}  # mu string -> record
        mp.dps = self.cfg['digits']
        getcontext().prec = self.cfg['digits']

    def seed(self, mu):
        known = sorted(self.orbits.values(), key=lambda r: mpf(r['mu']))
        if len(known) >= 2:  # linear interpolation/extrapolation from the two nearest MP orbits
            near = sorted(known, key=lambda r: abs(mpf(r['mu']) - mpf(mu)))[:2]
            m0, m1 = mpf(near[0]['mu']), mpf(near[1]['mu'])
            a = (mpf(mu) - m0) / (m1 - m0)
            return [im.s((1 - a) * mpf(near[0][k]) + a * mpf(near[1][k])) for k in ['u', 'v', 'T']], 'MP interpolation'
        if known:
            return [known[0][k] for k in ['u', 'v', 'T']], 'nearest MP orbit'
        seed, source = im.double_seed(mu)
        return seed, source

    def __call__(self, mu):
        if mu in self.orbits:
            return self.orbits[mu]
        cfg = self.cfg
        seed, source = self.seed(mu)
        prefix = f"{self.label}_{mu.replace('.', 'p')}"
        started = time.perf_counter()
        if cfg['backend'] == 'taylor':
            p, history = hp.correct(mu, seed, cfg['digits'], cfg['order'], cfg['tol'], self.raw_dir, prefix,
                                    target=cfg['target'], stride=10 ** 9)
            dat = hp.verify_orbit(mu, p, history, cfg['digits'], cfg['vorder'], cfg['vtol'], self.raw_dir, prefix, stride=50)
            rec = {'mu': mu, 'u': dat['u'], 'v': dat['v'], 'T': dat['T'], 'D': dat['floquet']['D'],
                   'max_closure': dat['max_closure'], 'spectral_radius': dat['floquet']['spectral_radius'],
                   'multipliers': [m['lambda'] for m in dat['floquet']['multipliers']],
                   'krein': [m['krein_signature'] for m in dat['floquet']['multipliers']],
                   'reciprocal_defect': dat['floquet']['reciprocal_defect'],
                   'raw_verification': dat['verification_integration']['raw_file']}
        else:
            args = Namespace(digits=cfg['digits'], stages=cfg['stages'], verify_stages=cfg['vstages'],
                             newton_tol=cfg['tol'], verify_tol=cfg['vtol'], target=cfg['target'], max_newton=10, stride=50)
            p, history = im.newton(mu, seed, args, self.raw_dir, prefix)
            dat = im.verify(mu, p, history, args, self.raw_dir, prefix)
            fl = dat['floquet']
            rec = {'mu': mu, 'u': dat['u'], 'v': dat['v'], 'T': dat['T'], 'D': fl['D_from_traces'],
                   'max_closure': dat['max_closure'], 'spectral_radius': fl['spectral_radius_nontrivial'],
                   'multipliers': [m['lambda'] for m in fl['modes']], 'krein': [m['krein_signature'] for m in fl['modes']],
                   'palindromic_defects': fl['palindromic_defects'], 'quartic_coefficients': fl['quartic_coefficients'],
                   'raw_verification': dat['verification_integration']['raw_file']}
        rec.update(seed_source=source, newton_iterations=len(history), runtime_seconds=time.perf_counter() - started)
        self.orbits[mu] = rec
        print(f"[{self.label}] mu={mu} D={mpmath.nstr(mpf(rec['D']), 15)} closure={rec['max_closure'][:10]}", flush=True)
        return rec


def value(rec, kind):
    """Transition function: D (Krein collision), q(-1) (period doubling) or q(+1) (tangent), where
    q(x) = prod_i (x - lambda_i) over the nontrivial quartet; q(-1) = (s1+2)(s2+2), q(1) = (2-s1)(2-s2)."""
    if kind == 'krein':
        return mpf(rec['D'])
    x = -1 if kind == 'pd' else 1
    if 'quartic_coefficients' in rec:  # GBS pipeline: coefficients from power sums of the full monodromy
        e1, e2, e3, e4 = map(mpf, rec['quartic_coefficients'])
        return x ** 4 - e1 * x ** 3 + e2 * x ** 2 - e3 * x + e4
    prod = mpmath.mpc(1)
    for re_, im_ in rec['multipliers']:
        prod *= (x - mpmath.mpc(mpf(re_), mpf(im_)))
    return mpmath.re(prod)


def refine(name, label, half_width, max_iter, out_dir, bracket=None, kind='krein'):
    ev = Evaluator(label, out_dir / 'raw' / f'{name}_{label}')
    a, b = bracket or BRACKETS[name]
    history = []
    F = lambda mu: value(ev(mu), kind)
    pts = [(a, F(a)), (b, F(b))]
    if pts[0][1] * pts[1][1] >= 0:
        raise RuntimeError(f'{label} {name}: initial bracket has no sign change')
    stop = mpf(10) ** (-(ev.cfg['digits'] // 2 + 5))
    for it in range(max_iter):
        (m0, d0), (m1, d1) = pts[-2], pts[-1]
        m_new = mpf(m1) - d1 * (mpf(m1) - mpf(m0)) / (d1 - d0)
        mu_new = mpmath.nstr(m_new, ev.cfg['digits'] - 10, strip_zeros=False)
        step = abs(mpf(mu_new) - mpf(m1))
        history.append({'iteration': it, 'mu': mu_new, 'secant_step': mpmath.nstr(step, 5)})
        pts.append((mu_new, F(mu_new)))
        if step < stop:
            break
    best = min(pts, key=lambda q: abs(q[1]))[0]
    hw = mpf(half_width)
    lo, hi = mpmath.nstr(mpf(best) - hw, ev.cfg['digits'] - 10, strip_zeros=False), mpmath.nstr(mpf(best) + hw, ev.cfg['digits'] - 10, strip_zeros=False)
    d_lo, d_hi = F(lo), F(hi)
    return {'boundary': name, 'kind': kind, 'pipeline': label, 'settings': SETTINGS[label], 'initial_bracket': [a, b],
            'secant_history': history, 'best_estimate': best, 'best_estimate_D': ev.orbits[best]['D'],
            'best_estimate_function_value': im.s(F(best), 12),
            'verified_bracket': [lo, hi], 'bracket_width': mpmath.nstr(2 * hw, 3),
            'bracket_D': [ev.orbits[lo]['D'], ev.orbits[hi]['D']],
            'bracket_function_values': [im.s(d_lo, 20), im.s(d_hi, 20)], 'sign_change': bool(d_lo * d_hi < 0),
            'evaluations': sorted(ev.orbits.values(), key=lambda r: mpf(r['mu']))}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--boundary', required=True, help='lower, upper, or a label for a custom --bracket')
    parser.add_argument('--bracket', help='a,b decimal masses (required for labels other than lower/upper)')
    parser.add_argument('--kind', choices=['krein', 'pd', 'tangent'], default='krein')
    parser.add_argument('--pipeline', choices=list(SETTINGS), required=True)
    parser.add_argument('--half-width', default='1e-30')
    parser.add_argument('--max-iter', type=int, default=8)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'results' / 'boundary_refinement')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    bracket = tuple(args.bracket.split(',')) if args.bracket else None
    if bracket is None and args.boundary not in BRACKETS:
        parser.error('--bracket is required for custom boundary labels')
    result = refine(args.boundary, args.pipeline, args.half_width, args.max_iter, args.output_dir, bracket, args.kind)
    result['runtime_seconds'] = time.perf_counter() - started
    result['status'] = 'empirically verified numerical sign bracket; not an interval-arithmetic enclosure'
    out = args.output_dir / f'{args.boundary}_{args.pipeline}.json'
    out.write_text(json.dumps(result, indent=1))
    print(f"{args.boundary} {args.pipeline} ({args.kind}): mu* = {result['best_estimate']}  "
          f"bracket values: {result['bracket_function_values']}  -> {out}")


if __name__ == '__main__':
    main()
