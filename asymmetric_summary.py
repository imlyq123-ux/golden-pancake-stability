"""Summarise the Phase E (1+delta, mu, 1-delta) binary64 scan: transition curves and window closing.

Only rows whose scan is *complete* (every grid point computed and every crossing refined; see
asymmetric_double.scan_mu) enter the edge fits, the window-closing estimate and the mirror check:
a missing or failed grid point can masquerade as instability (negative D) and is therefore
never interpreted. Incomplete rows are listed with their failures. Figures: make_figures.py.
"""
import argparse
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def row_complete(row):
    cov = row.get('coverage') or {}
    return bool(cov.get('complete')) and all(c.get('refined', False) for c in row.get('crossings', []))


def summarize(data, mu0):
    rows = sorted(data['rows'], key=lambda r: r['delta'])
    requested = data.get('deltas_requested') or [r['delta'] for r in rows]
    curves, theta_max, closure_max, incomplete = [], 0.0, 0.0, []
    for r in rows:
        ok = row_complete(r)
        for g in r['grid']:
            theta_max = max(theta_max, abs(g['theta']))
            closure_max = max(closure_max, g['closure'])
        if not ok:
            incomplete.append({'delta': r['delta'], 'failures': r.get('failures', 'no coverage record'),
                               'missing_points': (r.get('coverage') or {}).get('missing_points')})
        refined = [c for c in r['crossings'] if c.get('refined')]
        lower = [c for c in refined if c['stable_side'] == 'above']
        upper = [c for c in refined if c['stable_side'] == 'below']
        curves.append({'delta': r['delta'], 'complete': ok,
                       'mu_lower': lower[0]['mu_star'] if lower else None,
                       'mu_upper': upper[-1]['mu_star'] if upper else None,
                       'all_crossings': [(c['mu_star'], c['stable_side']) for c in refined],
                       'max_D': max((g['D'] for g in r['grid']), default=None),
                       'stable_points': sum(g['stable'] for g in r['grid']),
                       'mu_of_max_D': max(r['grid'], key=lambda g: g['D'])['mu'] if r['grid'] else None})
    usable = [c for c in curves if c['complete'] and c['mu_lower'] is not None and c['mu_upper'] is not None]
    fits = {}
    for key, m0 in zip(['mu_lower', 'mu_upper'], mu0):
        sel = [c for c in usable if c['delta'] <= 0.006]
        d = np.array([c['delta'] for c in sel]); m = np.array([c[key] for c in sel])
        entry = {'symmetric_value': m0, 'rows_used': len(sel)}
        if len(sel) >= 3:
            # mu*(delta) - mu*(0) = c2 delta^2 + c4 delta^4 (even in delta by the 1<->3 relabelling symmetry)
            A = np.column_stack((d ** 2, d ** 4))
            coef, *_ = np.linalg.lstsq(A, m - m0, rcond=None)
            entry.update(c2=float(coef[0]), c4=float(coef[1]), max_fit_residual=float(np.max(abs(A @ coef - (m - m0)))),
                         fit_range='delta <= 0.006')
        small = [c for c in usable if 0 < c['delta'] <= 0.0010001]
        q = [(c[key] - m0) / c['delta'] ** 2 for c in small[:2]]
        entry['difference_quotients'] = q
        entry['c2_small_delta_richardson'] = (4 * q[0] - q[1]) / 3 if len(q) == 2 else None
        zero = next((c for c in curves if c['delta'] == 0 and c['complete']), None)
        entry['delta0_scan_value'] = zero[key] if zero else None
        entry['delta0_scan_error'] = (zero[key] - m0) if zero and zero[key] is not None else None
        fits[key] = entry
    # window closing: max D changes sign between two *consecutive requested* deltas whose rows are complete
    by_delta = {c['delta']: c for c in curves}
    closing, blocked = None, []
    for d0, d1 in zip(requested, requested[1:]):
        c0, c1 = by_delta.get(d0), by_delta.get(d1)
        if c0 is None or c1 is None or not (c0['complete'] and c1['complete']):
            if c0 and c1 and c0['max_D'] is not None and c1['max_D'] is not None and c0['max_D'] > 0 >= c1['max_D']:
                blocked.append([d0, d1])
            continue
        if c0['max_D'] > 0 >= c1['max_D']:
            w = c0['max_D'] / (c0['max_D'] - c1['max_D'])
            closing = {'delta_between': [d0, d1], 'delta_linear_estimate': d0 + w * (d1 - d0), 'mu_near': c0['mu_of_max_D']}
    mirror = data.get('mirror_check')
    mirror_cmp = None
    if mirror and row_complete(mirror):
        ref = next((c for c in curves if abs(c['delta'] + mirror['delta']) < 1e-12 and c['complete']), None)
        if ref is not None:
            mine = [c['mu_star'] for c in mirror['crossings'] if c.get('refined')]
            theirs = [x[0] for x in ref['all_crossings']]
            mirror_cmp = {'delta': mirror['delta'], 'mirror_crossings': mine, 'positive_delta_crossings': theirs,
                          'max_difference': max(abs(a - b) for a, b in zip(mine, theirs)) if len(mine) == len(theirs) else None}
    return {'status': data['status'], 'all_rows_complete': not incomplete and not data.get('centre_continuation_failure'),
            'centre_continuation_failure': data.get('centre_continuation_failure'),
            'incomplete_rows': incomplete, 'closing_sign_changes_not_interpreted_due_to_incomplete_rows': blocked,
            'max_abs_theta': theta_max, 'max_closure': closure_max, 'curves': curves,
            'quadratic_fits_mu_star_vs_delta2': fits, 'window_closing': closing, 'mirror_symmetry_check': mirror_cmp}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--input', type=Path, default=ROOT / 'results' / 'asymmetric' / 'main_window.json')
    parser.add_argument('--mu0', default='0.59904378420889985122,0.61784426820906317873',
                        help='symmetric (delta = 0) transition masses from the multiprecision refinement')
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    summary = summarize(data, [float(x) for x in args.mu0.split(',')])
    out = args.input.with_name(args.input.stem + '_summary.json')
    out.write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != 'curves'}, indent=1))
    for c in summary['curves']:
        print(f"delta={c['delta']:.4f} complete={c['complete']} mu_L={c['mu_lower']} mu_R={c['mu_upper']} "
              f"maxD={c['max_D']} stable_pts={c['stable_points']}")


if __name__ == '__main__':
    main()
