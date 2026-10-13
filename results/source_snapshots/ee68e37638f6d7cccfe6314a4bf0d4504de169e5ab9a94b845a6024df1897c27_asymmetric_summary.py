"""Summarise the Phase E (1+delta, mu, 1-delta) binary64 scan: transition curves and stability region."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'results' / 'asymmetric' / 'main_window.json')
    parser.add_argument('--figure', type=Path, default=ROOT / 'figures' / 'asymmetric_stability_region.png')
    parser.add_argument('--mu0', default='0.59904378420889985122,0.61784426820906317873',
                        help='symmetric (delta = 0) transition masses from the multiprecision refinement')
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    rows = sorted(data['rows'], key=lambda r: r['delta'])
    mu0 = [float(x) for x in args.mu0.split(',')]
    curves, theta_max, closure_max = [], 0.0, 0.0
    for r in rows:
        for g in r['grid']:
            theta_max = max(theta_max, abs(g['theta']))
            closure_max = max(closure_max, g['closure'])
        lower = [c for c in r['crossings'] if c['stable_side'] == 'above']
        upper = [c for c in r['crossings'] if c['stable_side'] == 'below']
        curves.append({'delta': r['delta'], 'mu_lower': lower[0]['mu_star'] if lower else None,
                       'mu_upper': upper[-1]['mu_star'] if upper else None,
                       'all_crossings': [(c['mu_star'], c['stable_side']) for c in r['crossings']],
                       'max_D': max(g['D'] for g in r['grid']), 'stable_points': sum(g['stable'] for g in r['grid']),
                       'mu_of_max_D': max(r['grid'], key=lambda g: g['D'])['mu']})
    open_rows = [c for c in curves if c['mu_lower'] is not None and c['mu_upper'] is not None]
    fits = {}
    for key, m0 in zip(['mu_lower', 'mu_upper'], mu0):
        d = np.array([c['delta'] for c in open_rows if c['delta'] <= 0.006])
        m = np.array([c[key] for c in open_rows if c['delta'] <= 0.006])
        # mu*(delta) - mu*(0) = c2 delta^2 + c4 delta^4 (even in delta by the 1<->3 relabelling symmetry)
        A = np.column_stack((d ** 2, d ** 4))
        coef, *_ = np.linalg.lstsq(A, m - m0, rcond=None)
        resid = float(np.max(abs(A @ coef - (m - m0))))
        small = [c for c in open_rows if 0 < c['delta'] <= 0.0010001]
        q = [(c[key] - m0) / c['delta'] ** 2 for c in small[:2]]
        c2_richardson = (4 * q[0] - q[1]) / 3 if len(q) == 2 else None  # q(delta) = c2 + c4 delta^2, deltas h and 2h
        fits[key] = {'c2_small_delta_richardson': c2_richardson, 'difference_quotients': q,
                     'c2': float(coef[0]), 'c4': float(coef[1]), 'max_fit_residual': resid, 'fit_range': 'delta <= 0.006',
                     'symmetric_value': m0,
                     'delta0_scan_value': curves[0][key], 'delta0_scan_error': (curves[0][key] - m0) if curves[0][key] else None}
    # closing of the window: last delta with stable points and first without (max D changes sign)
    closing = None
    for c0, c1 in zip(curves, curves[1:]):
        if c0['max_D'] > 0 >= c1['max_D']:
            w = c0['max_D'] / (c0['max_D'] - c1['max_D'])
            closing = {'delta_between': [c0['delta'], c1['delta']], 'delta_linear_estimate': c0['delta'] + w * (c1['delta'] - c0['delta']),
                       'mu_near': c0['mu_of_max_D']}
    mirror = data.get('mirror_check')
    mirror_cmp = None
    if mirror and 'crossings' in mirror:
        ref = next(c for c in curves if abs(c['delta'] + mirror['delta']) < 1e-12)
        mirror_cmp = {'delta': mirror['delta'], 'mirror_crossings': [c['mu_star'] for c in mirror['crossings']],
                      'positive_delta_crossings': [x[0] for x in ref['all_crossings']],
                      'max_difference': max(abs(a - b) for a, b in zip([c['mu_star'] for c in mirror['crossings']],
                                                                           [x[0] for x in ref['all_crossings']]))}
    summary = {'status': data['status'], 'max_abs_theta': theta_max, 'max_closure': closure_max,
               'curves': curves, 'quadratic_fits_mu_star_vs_delta2': fits, 'window_closing': closing,
               'mirror_symmetry_check': mirror_cmp}
    out = args.input.with_name(args.input.stem + '_summary.json')
    out.write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != 'curves'}, indent=1))
    for c in curves:
        print(f"delta={c['delta']:.4f} mu_L={c['mu_lower']} mu_R={c['mu_upper']} maxD={c['max_D']:+.4f} stable_pts={c['stable_points']}")
    # figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    mus = sorted({g['mu'] for r in rows for g in r['grid']})
    deltas = [r['delta'] for r in rows]
    Z = np.full((len(deltas), len(mus)), np.nan)
    for i, r in enumerate(rows):
        for g in r['grid']:
            Z[i, mus.index(g['mu'])] = g['D']
    im = axes[0].pcolormesh(mus, deltas, np.clip(Z, -1, 1), cmap='RdBu', vmin=-1, vmax=1, shading='nearest')
    axes[0].contour(mus, deltas, Z, levels=[0], colors='k', linewidths=1)
    for sign in (1, -1):
        axes[0].plot([c['mu_lower'] for c in open_rows], [sign * c['delta'] for c in open_rows], 'g.', ms=4)
        axes[0].plot([c['mu_upper'] for c in open_rows], [sign * c['delta'] for c in open_rows], 'm.', ms=4)
    axes[0].set(xlabel='middle mass mu', ylabel='outer-mass asymmetry delta  (masses 1+delta, mu, 1-delta)',
                title='D(mu, delta) (clipped); black: D = 0; points: refined transitions (mirrored to delta<0)')
    fig.colorbar(im, ax=axes[0], label='D (clipped to [-1, 1])')
    axes[1].plot([c['delta'] for c in open_rows], [c['mu_lower'] for c in open_rows], 'g.-', label='lower transition')
    axes[1].plot([c['delta'] for c in open_rows], [c['mu_upper'] for c in open_rows], 'm.-', label='upper transition')
    axes[1].set(xlabel='delta', ylabel='transition mass mu*(delta)', title='Main-window boundaries vs asymmetry')
    axes[1].legend(); axes[1].grid(alpha=.2)
    fig.tight_layout(); fig.savefig(args.figure, dpi=170); plt.close(fig)


if __name__ == '__main__':
    main()
