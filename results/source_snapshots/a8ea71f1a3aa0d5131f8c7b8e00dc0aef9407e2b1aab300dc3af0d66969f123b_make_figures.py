"""Figures for the second-session results (transition scaling, isola stability map, Phase E region).

Palette: reference data-viz palette (categorical slots 1-3, status good/critical with text labels,
blue<->red diverging with a gray midpoint). Inputs are the JSON results; nothing is recomputed.
"""
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpmath import mpf

ROOT = Path(__file__).resolve().parent
RES, FIG = ROOT / 'results', ROOT / 'figures'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a']          # categorical slots 1-3
GOOD, CRITICAL = '#0ca30c', '#d03b3b'               # status: stable / unstable (always labelled)
INK, INK2, MUTED, GRID, BASE, SURFACE = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7', '#fcfcfb'
DIVERGING = LinearSegmentedColormap.from_list('blue_gray_red', ['#104281', '#3987e5', '#f0efec', '#e34948', '#8c2424'])

plt.rcParams.update({'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'axes.edgecolor': BASE,
                     'axes.labelcolor': INK2, 'xtick.color': MUTED, 'ytick.color': MUTED, 'text.color': INK,
                     'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': .6, 'grid.linestyle': '-',
                     'lines.linewidth': 1.5, 'font.size': 9.5, 'legend.frameon': False, 'axes.titlesize': 10.5})


def transition_scaling():
    data = json.loads((RES / 'transition_analysis.json').read_text())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for ax, name, title in zip(axes, ['lower', 'upper'], ['K1 = mu_L (lower edge of main window)', 'K2 = mu_R (upper edge)']):
        rows = [r for r in data['rows'] if r['boundary'] == name and r['side'] != 'star']
        st = sorted((r for r in rows if r['side'] == 'stable'), key=lambda r: r['delta_exponent'])
        un = sorted((r for r in rows if r['side'] == 'unstable'), key=lambda r: r['delta_exponent'])
        d_st = [10.0 ** -r['delta_exponent'] for r in st]
        d_un = [10.0 ** -r['delta_exponent'] for r in un]
        split = [float(mpf(r['split_abs'])) for r in st]
        rho = [float(max(map(mpf, r['moduli_minus_1']))) for r in un]
        kap = [float(sum(abs(mpf(k)) for k in r['krein_forms_normalised']) / 2) for r in st]
        ax.loglog(d_st, split, 'o-', color=SERIES[0], ms=4.5, label='|lambda_1 - lambda_2|  (stable side)')
        ax.loglog(d_un, rho, 's-', color=SERIES[1], ms=4.5, label='rho - 1  (unstable side)')
        ax.loglog(d_st, kap, '^-', color=SERIES[2], ms=4.5, label='|Krein form|  (stable side)')
        ref = np.array([1e-40, 1e-4])
        ax.loglog(ref, 0.5 * np.sqrt(ref), color=MUTED, lw=1)
        ax.text(1e-22, 0.5 * np.sqrt(1e-22) / 12, 'reference slope 1/2', color=INK2, rotation=16, ha='center')
        s = data['summary'][name]
        ax.set(xlabel='delta = |mu - mu*|', title=title)
        ax.text(0.03, 0.04, f"fitted slopes: {s['loglog_slope_stable_split']:.4f}, {s['loglog_slope_unstable_rho_minus_1']:.4f}, "
                f"{s['loglog_slope_krein_form']:.4f}\nsigma_min / sigma_2 at mu*: "
                f"{float(mpf(s['singular_values_at_mu_star'][0])):.1e} / {float(mpf(s['singular_values_at_mu_star'][1])):.3f}",
                transform=ax.transAxes, color=INK2, fontsize=8.5)
    axes[0].set_ylabel('magnitude')
    axes[0].legend(loc='upper left', fontsize=8.5)
    fig.suptitle('Square-root unfolding at both main transitions: non-semisimple Krein collisions (80-digit GBS)', color=INK)
    fig.tight_layout(); fig.savefig(FIG / 'transition_scaling.png', dpi=180); plt.close(fig)


def stable(r):
    return r['D'] > 0 and all(abs(x[1]) < 1e-9 and -2 < x[0] < 2 for x in r['s_values']) and r['spectral_radius'] < 1 + 1e-6


def isola_overview():
    scans = {k: [r for r in json.loads((RES / 'isola_scan' / f'{k}.json').read_text())['records'] if 'D' in r]
             for k in ['W_down', 'W_up', 'C_down', 'C_up']}
    W = scans['W_down'][::-1] + scans['W_up'][1:]
    C = scans['C_down'][::-1] + scans['C_up'][1:]
    trans = json.loads((RES / 'stability_transitions.json').read_text()) if (RES / 'stability_transitions.json').exists() else {}
    fig = plt.figure(figsize=(15.5, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.3, 0.62])
    ax = fig.add_subplot(gs[0])
    for curve, name, xy in [(W, 'AABAbbaabaBB segment', (0.1545, 0.592)), (C, '(ABab)^3 segment', (0.2385, 0.566))]:
        u = np.array([r['u'] for r in curve]); mu = np.array([r['mu'] for r in curve])
        st = np.array([stable(r) for r in curve])
        ax.plot(u, mu, color=BASE, lw=1.2, zorder=1)
        ax.plot(u, np.where(st, mu, np.nan), color=GOOD, lw=3.2, solid_capstyle='round', zorder=2)
        ax.plot([u[0], u[-1]], [mu[0], mu[-1]], 'o', mfc=SURFACE, mec=INK2, ms=6, zorder=3)
        ax.text(*xy, name, color=INK2, fontsize=8.5, rotation=0)
    ax.plot([], [], color=GOOD, lw=3.2, label='spectrally stable (all four nontrivial multipliers on |lambda| = 1)')
    ax.plot([], [], color=BASE, lw=1.2, label='unstable')
    ax.plot([], [], 'o', mfc=SURFACE, mec=INK2, ms=6, label='scan end (binary64 limit; collision limits not reached)')
    ax.plot(0.151837333844166, 0.612356088827995, marker='*', color=INK, ms=10, ls='none',
            label='atlas representative (Golden Pancake Drama)')
    ax.set(xlabel='u  (Euler-section velocity component)', ylabel='middle mass mu',
           title='Isola in the (u, mu) plane (binary64 pseudo-arclength scan)')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=7.6)
    # right: W-segment pair sums versus mu, with an inset for the fold region
    ax = fig.add_subplot(gs[1])
    near = [r for r in W if r['arclength'] < 0.1060 or r is W[0] or r['mu'] < 0.634]
    far = [r for r in scans['W_up'] if r['arclength'] >= 0.1060]
    def draw(axis, rows, marker_real='o', marker_cplx='o', size=2.2, label=True):
        mu = np.array([r['mu'] for r in rows])
        sv = np.array([[x[0] for x in r['s_values']] for r in rows]); si = np.array([[x[1] for x in r['s_values']] for r in rows])
        real = np.all(abs(si) < 1e-9, axis=1)
        axis.plot(np.repeat(mu[real], 2), sv[real].ravel(), marker_real, ms=size, color=SERIES[0], ls='none',
                  label='pair sums s = lambda + 1/lambda (both real)' if label else None)
        axis.plot(mu[~real], sv[~real, 0], marker_cplx, ms=size, color=SERIES[1], ls='none',
                  label='Re s of a complex quartet (unstable)' if label else None)
    ax.axhspan(-2, 2, color='#cde2fb', alpha=.35, lw=0, zorder=0)
    draw(ax, near)
    for y in (-2, 2):
        ax.axhline(y, color=BASE, lw=.9)
    marks = {'K0_lowwindow_top': 'K0', 'lower': 'K1', 'upper': 'K2'}
    for key, lab in marks.items():
        if key in trans:
            m = float(mpf(next(iter(trans[key]['best_estimates'].values()))))
            ax.axvline(m, color=MUTED, lw=.8)
            ax.text(m, 3.55, lab, color=INK2, ha='center', fontsize=8.5)
    ax.text(0.6339, 3.55, 'K3, fold:\nsee right', color=INK2, ha='center', fontsize=8, va='center')
    ax.annotate('touches of s = -2 (N-multipliers at +-i)', (0.5695, -2.0), xytext=(0.548, -2.75), color=INK2, fontsize=8,
                arrowprops=dict(arrowstyle='-', color=MUTED, lw=.8))
    ax.set(xlabel='middle mass mu', ylabel='s', ylim=(-3.2, 3.9), xlim=(0.536, 0.6349),
           title='AABAbbaabaBB segment: stable iff both s real and in (-2, 2) (shaded band)')
    ax.legend(loc='upper left', fontsize=7.8, bbox_to_anchor=(0.0, 0.93))
    ins = fig.add_subplot(gs[2], sharey=ax)
    ins.axhspan(-2, 2, color='#cde2fb', alpha=.35, lw=0)
    draw(ins, [r for r in near if r['mu'] > 0.6336], size=2.6, label=False)
    draw(ins, [r for r in far if r['mu'] > 0.6336], marker_real='x', marker_cplx='x', size=3.2, label=False)
    for y in (-2, 2):
        ins.axhline(y, color=BASE, lw=.8)
    if 'K3_nearfold' in trans:
        m = float(mpf(next(iter(trans['K3_nearfold']['best_estimates'].values()))))
        ins.axvline(m, color=MUTED, lw=.8); ins.text(m, 3.55, 'K3', color=INK2, ha='center', fontsize=8.5)
    ins.axvline(0.63460425, color=MUTED, lw=.8); ins.text(0.63460425, 3.55, 'fold', color=INK2, ha='center', fontsize=8.5)
    ins.set_xlim(0.6338, 0.63466)
    ins.plot([], [], 'x', color=SERIES[0], ms=4, label='past the fold (real s)')
    ins.legend(loc='lower left', fontsize=7.5)
    ins.tick_params(axis='x', labelsize=7, labelrotation=30)
    ins.set(xlabel='middle mass mu', title='Fold region (zoom)')
    fig.tight_layout(); fig.savefig(FIG / 'isola_stability_overview.png', dpi=180, bbox_inches='tight'); plt.close(fig)


def asymmetric_region():
    path = RES / 'asymmetric' / 'main_window.json'
    if not path.exists():
        return
    data = json.loads(path.read_text())
    rows = sorted(data['rows'], key=lambda r: r['delta'])
    mus = sorted({g['mu'] for r in rows for g in r['grid']})
    deltas = np.array([r['delta'] for r in rows])
    Z = np.full((len(deltas), len(mus)), np.nan)
    for i, r in enumerate(rows):
        for g in r['grid']:
            Z[i, mus.index(g['mu'])] = g['D']
    full_d = np.concatenate((-deltas[::-1][:-1], deltas))
    full_Z = np.vstack((Z[::-1][:-1], Z))       # mirrored: masses (1+delta, mu, 1-delta) ~ (1-delta, mu, 1+delta)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={'width_ratios': [1.25, 1]})
    ax = axes[0]
    im = ax.pcolormesh(mus, full_d, np.clip(full_Z, -1, 1), cmap=DIVERGING.reversed(), vmin=-1, vmax=1,
                       shading='nearest', rasterized=True)  # red: D < 0 (complex quartet), blue: D > 0
    ax.contour(mus, full_d, full_Z, levels=[0], colors=INK, linewidths=1.1)
    from asymmetric_summary import row_complete
    complete = [r for r in rows if row_complete(r)]
    for r in complete:
        for c in r['crossings']:
            for sgn in (1, -1):
                ax.plot(c['mu_star'], sgn * r['delta'], 'o', ms=2.6, color=INK, zorder=3)
    ax.set(xlabel='middle mass mu', ylabel='delta   (masses 1+delta, mu, 1-delta)',
           title='Discriminant D(mu, delta); black curve and dots: D = 0')
    cb = fig.colorbar(im, ax=ax, label='D (clipped to [-1, 1]);  D < 0: complex quartet (unstable)')
    cb.outline.set_edgecolor(BASE)
    ax.text(0.6085, 0.0, 'stable\nwindow', color=INK, ha='center', va='center', fontsize=9)
    summ = path.with_name('main_window_summary.json')
    if summ.exists():
        sj = json.loads(summ.read_text())
        dc = sj['window_closing']['delta_linear_estimate']
        c2 = [sj['quadratic_fits_mu_star_vs_delta2'][k]['c2_small_delta_richardson'] for k in ('mu_lower', 'mu_upper')]
        closing_text = f'edges meet (window closes)\nnear |delta| = {dc:.4f}, mu = 0.611'
    ax = axes[1]
    lower = [(r['delta'], c['mu_star']) for r in complete for c in r['crossings'] if c['stable_side'] == 'above']
    upper = [(r['delta'], c['mu_star']) for r in complete for c in r['crossings'] if c['stable_side'] == 'below']
    if len(complete) < len(rows):
        ax.set_title(f'{len(rows) - len(complete)} incomplete scan rows omitted from edges and dots', fontsize=8, loc='left')
    ax.plot([d for d, _ in lower], [m for _, m in lower], 'o-', color=SERIES[0], ms=3.5, label='lower transition')
    ax.plot([d for d, _ in upper], [m for _, m in upper], 's-', color=SERIES[1], ms=3.5, label='upper transition')
    ax.set(xlabel='delta', ylabel='transition mass mu*(delta)', title='Window edges versus outer-mass asymmetry')
    if summ.exists():
        ax.annotate(closing_text, (dc, 0.6108), xytext=(0.0047, 0.6105), color=INK2, fontsize=8.5,
                    arrowprops=dict(arrowstyle='-', color=MUTED, lw=.8))
        ax.text(0.0002, 0.6035, f'small delta: mu_L + {c2[0]:.1f} delta^2\n                 mu_R {c2[1]:+.1f} delta^2',
                color=INK2, fontsize=8.5)
    ax.legend(loc='upper left', fontsize=8.5)
    fig.tight_layout(); fig.savefig(FIG / 'asymmetric_stability_region.png', dpi=180); plt.close(fig)


if __name__ == '__main__':
    transition_scaling()
    isola_overview()
    asymmetric_region()
    print('figures written')
