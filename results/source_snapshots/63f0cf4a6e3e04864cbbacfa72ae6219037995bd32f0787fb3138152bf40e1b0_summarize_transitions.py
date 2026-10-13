"""Consolidate multiprecision transition refinements (results/boundary_refinement/*.json).

For every transition label the best estimates of all pipelines are compared digit by digit,
the verified sign brackets are listed with their endpoint function values, and the local
slope of the transition function is reported (nonzero slope = transversal crossing).
Outputs results/stability_transitions.json and results/stability_transitions.md.
"""
import json
from collections import defaultdict
from pathlib import Path

import mpmath
from mpmath import mp, mpf

ROOT = Path(__file__).resolve().parent
mp.dps = 100
ORDER = ['taylor50', 'taylor70', 'taylor90', 'gbs80']
LABELS = {'lower': 'K1 = mu_L (lower edge of main window)', 'upper': 'K2 = mu_R (upper edge of main window)',
          'K0_lowwindow_top': 'K0 (upper edge of lower stable stretch)',
          'K3_nearfold': 'K3 (lower edge of near-fold stable stretch)',
          'T1_abab_postfold': 'T1 ((ABab)^3 segment, s = +2 crossing past mu_min)'}


def agreeing_digits(a, b):
    diff = abs(mpf(a) - mpf(b))
    return 999 if diff == 0 else int(mpmath.floor(-mpmath.log10(diff / abs(mpf(a)))))


def main():
    groups = defaultdict(dict)
    for path in sorted((ROOT / 'results' / 'boundary_refinement').glob('*.json')):
        data = json.loads(path.read_text())
        groups[data['boundary']][data['pipeline']] = data
    out = {}
    lines = ['# Multiprecision spectral transitions', '',
             'Each row: one transition refined independently by each pipeline (secant on the transition',
             'function, every evaluation a fresh Newton-corrected orbit plus separate monodromy integration).',
             'Brackets are empirically verified sign changes, not interval-arithmetic enclosures.', '',
             '| transition | kind | pipelines | estimate (digits confirmed by two independent pipelines) | best / worst pairwise agreement (digits) | bracket half-width | endpoint values | slope |',
             '|---|---|---|---|---|---|---|---|']
    for name, runs in sorted(groups.items(), key=lambda kv: mpf(next(iter(kv[1].values()))['best_estimate'])):
        pipes = [p for p in ORDER if p in runs]
        ref = runs[pipes[-1]]
        best = {p: runs[p]['best_estimate'] for p in pipes}
        pair_digits = {f'{a}~{b}': agreeing_digits(best[a], best[b]) for i, a in enumerate(pipes) for b in pipes[i + 1:]}
        min_digits = min(pair_digits.values()) if pair_digits else None
        max_digits = max(pair_digits.values()) if pair_digits else None
        best_pair = max(pair_digits, key=pair_digits.get) if pair_digits else None
        # digits shown: the most precise estimate, truncated to the best agreement between two pipelines
        precise = runs['taylor90'] if 'taylor90' in runs else ref
        text = precise['best_estimate']          # decimal string 0.xxxx...; truncate (never round) to the agreed digits
        common = text[:2 + min(max_digits or 30, 70)] + '...'
        brackets = {p: {'bracket': runs[p]['verified_bracket'],
                        'values': runs[p].get('bracket_function_values', runs[p]['bracket_D']),
                        'sign_change': runs[p]['sign_change']} for p in pipes}
        a, b = ref['verified_bracket']
        va, vb = map(mpf, ref.get('bracket_function_values', ref['bracket_D']))
        slope = (vb - va) / (mpf(b) - mpf(a))
        out[name] = {'description': LABELS.get(name, name), 'kind': ref.get('kind', 'krein'), 'pipelines': pipes,
                     'best_estimates': best, 'pairwise_agreeing_digits': pair_digits,
                     'best_agreeing_pair': best_pair, 'best_agreement_digits': max_digits,
                     'common_digits_estimate': common, 'verified_brackets': brackets,
                     'all_brackets_change_sign': all(v['sign_change'] for v in brackets.values()),
                     'transition_function_slope': mpmath.nstr(slope, 10)}
        lines.append(f"| {LABELS.get(name, name)} | {out[name]['kind']} | {', '.join(pipes)} | {common} | "
                     f"{max_digits} ({best_pair}) / {min_digits} | {mpmath.nstr((mpf(b) - mpf(a)) / 2, 2)} | "
                     f"{mpmath.nstr(va, 4)}, {mpmath.nstr(vb, 4)} | {mpmath.nstr(slope, 8)} |")
    (ROOT / 'results' / 'stability_transitions.json').write_text(json.dumps(out, indent=1))
    (ROOT / 'results' / 'stability_transitions.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
