"""Summarise the binary64 isola scans: stable stretches, transition types, figures.

A point is a spectral-stability candidate when both pair sums s = lambda + 1/lambda of the
nontrivial quartet are real and lie in (-2, 2) (all four multipliers on the unit circle).
Transitions between consecutive scan points are classified as
  krein          D changes sign (pair sums collide; complex quartet on one side)
  period-doubling a real pair sum crosses -2 (multiplier pair through -1)
  tangent         a real pair sum crosses +2 (multiplier pair through +1, e.g. a fold in mu)
Locations are linear interpolations between binary64 scan points (exploratory only).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def load(paths):
    segments = {}
    for path in paths:
        data = json.loads(path.read_text())
        recs = [r for r in data['records'] if 'D' in r]
        segments[path.stem] = {'records': recs, 'stop': [r for r in data['records'] if 'stop' in r]}
    return segments


def merged(segments, a, b):
    """Join the 'down' and 'up' scans started at the same seed into one ordered curve."""
    down = segments[a]['records'][::-1]
    up = segments[b]['records'][1:]
    return down + up


def stable(r, tol=1e-6):
    s = r['s_values']
    return r['D'] > 0 and all(abs(x[1]) < 1e-9 and -2 < x[0] < 2 for x in s) and r['spectral_radius'] < 1 + tol


def transitions(curve):
    found = []
    for a, b in zip(curve, curve[1:]):
        t = lambda fa, fb: fa / (fa - fb)
        if (a['D'] > 0) != (b['D'] > 0):
            w = t(a['D'], b['D'])
            found.append({'type': 'krein', 'mu': a['mu'] + w * (b['mu'] - a['mu']),
                          'between': [a['mu'], b['mu']], 'D': [a['D'], b['D']],
                          's_at_collision': (a['s_values'][0][0] + a['s_values'][1][0]) / 2})
        if a['D'] > 0 and b['D'] > 0:
            sa = sorted(x[0] for x in a['s_values']); sb = sorted(x[0] for x in b['s_values'])
            for k in range(2):
                for level, kind in [(-2, 'period-doubling'), (2, 'tangent')]:
                    if (sa[k] - level) * (sb[k] - level) < 0:
                        w = t(sa[k] - level, sb[k] - level)
                        found.append({'type': kind, 'mu': a['mu'] + w * (b['mu'] - a['mu']),
                                      'between': [a['mu'], b['mu']], 's': [sa[k], sb[k]]})
    return found


def intervals(curve):
    out, start = [], None
    for prev, r in zip([None] + curve[:-1], curve):
        if stable(r) and start is None:
            start = (prev, r)
        if not stable(r) and start is not None:
            out.append({'first_stable_point': start[1]['mu'], 'last_stable_point': prev['mu'],
                        'outer_unstable_points': [start[0]['mu'] if start[0] else None, r['mu']],
                        'points': None})
            start = None
    if start is not None:
        out.append({'first_stable_point': start[1]['mu'], 'last_stable_point': curve[-1]['mu'],
                    'outer_unstable_points': [start[0]['mu'] if start[0] else None, None]})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--scan-dir', type=Path, default=ROOT / 'results' / 'isola_scan')
    parser.add_argument('--figure-dir', type=Path, default=ROOT / 'figures')
    parser.add_argument('--output', type=Path, default=ROOT / 'results' / 'isola_scan_summary.json')
    args = parser.parse_args()
    segments = load(sorted(args.scan_dir.glob('*.json')))
    curves = {}
    if {'W_down', 'W_up'} <= set(segments):
        curves['AABAbbaabaBB (W) segment'] = merged(segments, 'W_down', 'W_up')
    if {'C_down', 'C_up'} <= set(segments):
        curves['(ABab)^3 segment'] = merged(segments, 'C_down', 'C_up')
    summary = {}
    for name, curve in curves.items():
        mus = [r['mu'] for r in curve]
        summary[name] = {'points': len(curve), 'mu_range': [min(mus), max(mus)],
                         'end_rmin': [curve[0]['sampled_rmin'], curve[-1]['sampled_rmin']],
                         'end_mu': [curve[0]['mu'], curve[-1]['mu']],
                         'max_closure': max(r['closure'] for r in curve),
                         'max_heuristic_vs_quotient': max(r['heuristic_vs_quotient_max_difference'] for r in curve),
                         'stable_points': sum(stable(r) for r in curve),
                         'stable_intervals': intervals(curve), 'transitions': transitions(curve),
                         'stops': {k: segments[k]['stop'] for k in segments if k[0] == name[0] or name.startswith('(') and k[0] == 'C'}}
    args.output.write_text(json.dumps(summary, indent=1))
    for name, s in summary.items():
        print(f"== {name}: {s['points']} points, mu in [{s['mu_range'][0]:.6f}, {s['mu_range'][1]:.6f}], "
              f"ends at rmin {s['end_rmin'][0]:.1e} / {s['end_rmin'][1]:.1e}")
        for it in s['stable_intervals']:
            print('   stable stretch: points', it['first_stable_point'], '..', it['last_stable_point'],
                  '(bounded by unstable points', it['outer_unstable_points'], ')')
        for t in s['transitions']:
            print(f"   {t['type']:>15s} near mu = {t['mu']:.7f} (between {t['between'][0]:.7f} and {t['between'][1]:.7f})")
    # Figures
    args.figure_dir.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    for name, curve in curves.items():
        u = np.array([r['u'] for r in curve]); mu = np.array([r['mu'] for r in curve])
        st = np.array([stable(r) for r in curve])
        axes[0].plot(u, mu, '-', color='0.6', lw=1)
        axes[0].plot(u[st], mu[st], 'o', ms=2.5, color='tab:green', label='spectrally stable' if name.startswith('AAB') else None)
        axes[0].plot(u[~st], mu[~st], 'o', ms=1.5, color='tab:red', label='unstable' if name.startswith('AAB') else None)
        rho = np.array([r['spectral_radius'] for r in curve])
        axes[1].semilogy(np.array([r['mu'] for r in curve]), rho, '.', ms=3, label=name)
    axes[0].set(xlabel='u (Euler-section velocity)', ylabel='middle mass mu', title='Isola coloured by binary64 Floquet scan')
    axes[0].legend(); axes[0].grid(alpha=.2)
    axes[1].set(xlabel='middle mass mu', ylabel='largest nontrivial |lambda|', title='Spectral radius along the isola')
    axes[1].legend(fontsize=8); axes[1].grid(alpha=.2)
    fig.tight_layout(); fig.savefig(args.figure_dir / 'isola_stability_map.png', dpi=180); plt.close(fig)
    if 'AABAbbaabaBB (W) segment' in curves:
        curve = curves['AABAbbaabaBB (W) segment']
        fig, ax = plt.subplots(figsize=(9, 5))
        mu = np.array([r['mu'] for r in curve])
        s = np.array([[x[0] for x in r['s_values']] for r in curve]); si = np.array([[x[1] for x in r['s_values']] for r in curve])
        real = np.all(abs(si) < 1e-9, axis=1)
        for k in range(2):
            ax.plot(mu[real], s[real, k], '.', ms=3, color='tab:blue')
        ax.plot(mu[~real], s[~real, 0], '.', ms=3, color='tab:red', label='Re s of complex pair (quartet off the circle)')
        ax.axhline(2, color='k', lw=.8, ls='--'); ax.axhline(-2, color='k', lw=.8, ls='--')
        ax.set(xlabel='middle mass mu', ylabel='pair sums s = lambda + 1/lambda',
               title='W segment: stability iff both s real in (-2, 2)')
        ax.set_ylim(-6, 6); ax.legend(fontsize=8); ax.grid(alpha=.2)
        fig.tight_layout(); fig.savefig(args.figure_dir / 'isola_pair_sums_W.png', dpi=180); plt.close(fig)


if __name__ == '__main__':
    main()
