"""Consolidate multiprecision transition refinements (results/boundary_refinement/*.json).

Independent confirmation is counted only between *different implementations*: a Taylor/Decimal
pipeline (taylor50/70/90) against the GBS/mpmath pipeline (gbs80). Agreement between two Taylor
runs measures precision consistency of one implementation and is reported separately; it never
raises the number of confirmed digits.

Confirmed digits = length of the common prefix of the two decimal strings (significant digits),
capped one below each string's printed length because the last printed digit may be rounded.
The displayed estimate is the more precise member of the confirming pair, truncated (never
rounded) to the confirmed digits. Outputs stability_transitions.json/.md in --output-dir.
"""
import argparse
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


def implementation(pipeline):
    return 'taylor' if pipeline.startswith('taylor') else 'gbs' if pipeline.startswith('gbs') else pipeline


def significant(text):
    """Significant-digit string of a positive decimal written without exponent."""
    return text.replace('.', '').lstrip('0')


def common_digits(a, b):
    """Confirmed significant digits shared by decimal strings a and b (last printed digit excluded)."""
    if mpmath.floor(mpmath.log10(mpf(a))) != mpmath.floor(mpmath.log10(mpf(b))):
        return 0
    da, db = significant(a), significant(b)
    n = 0
    for x, y in zip(da, db):
        if x != y:
            break
        n += 1
    return max(0, min(n, len(da) - 1, len(db) - 1))


def truncate(text, n):
    """Truncate a positive decimal string to n significant digits (no rounding)."""
    out, count, started = [], 0, False
    for ch in text:
        if ch == '.':
            out.append(ch)
            continue
        if count >= n:
            break
        out.append(ch)
        if ch != '0' or started:
            started = True
            count += 1
    return ''.join(out).rstrip('.') + '...'


def confirmation(runs):
    """runs: pipeline -> best-estimate string. Independent (cross-implementation) and same-implementation agreement."""
    pipes = [p for p in ORDER if p in runs] + sorted(p for p in runs if p not in ORDER)
    cross, same = {}, {}
    for i, a in enumerate(pipes):
        for b in pipes[i + 1:]:
            target = cross if implementation(a) != implementation(b) else same
            target[f'{a}~{b}'] = common_digits(runs[a], runs[b])
    if cross:
        best_pair = max(cross, key=cross.get)
        digits = cross[best_pair]
        a, b = best_pair.split('~')
        precise = a if len(significant(runs[a])) >= len(significant(runs[b])) else b
        display = truncate(runs[precise], digits)
    else:
        best_pair, digits, display = None, None, None
    return {'pipelines': pipes, 'independent_pairs': cross, 'same_implementation_pairs': same,
            'confirmed_digits': digits, 'confirming_pair': best_pair, 'confirmed_estimate': display}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--refinement-dir', type=Path, default=ROOT / 'results' / 'boundary_refinement')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'results')
    args = parser.parse_args()
    groups = defaultdict(dict)
    for path in sorted(args.refinement_dir.glob('*.json')):
        data = json.loads(path.read_text())
        groups[data['boundary']][data['pipeline']] = data
    out = {}
    lines = ['# Multiprecision spectral transitions', '',
             'Each row: one transition refined by several pipelines (secant on the transition function; every',
             'evaluation is a fresh Newton-corrected orbit plus a separate monodromy integration). Confirmed digits',
             'count only Taylor-vs-GBS agreement (independent implementations); Taylor-vs-Taylor agreement is a',
             'precision-consistency check. Brackets are empirically verified sign changes, not interval enclosures.', '',
             '| transition | kind | pipelines | confirmed estimate (truncated) | confirmed digits (pair) | '
             'Taylor-Taylor agreement | bracket half-width | endpoint values | slope |',
             '|---|---|---|---|---|---|---|---|---|']
    for name, runs in sorted(groups.items(), key=lambda kv: mpf(next(iter(kv[1].values()))['best_estimate'])):
        conf = confirmation({p: d['best_estimate'] for p, d in runs.items()})
        ref = runs[conf['pipelines'][-1]]
        brackets = {p: {'bracket': runs[p]['verified_bracket'],
                        'values': runs[p].get('bracket_function_values', runs[p]['bracket_D']),
                        'sign_change': runs[p]['sign_change']} for p in conf['pipelines']}
        a, b = ref['verified_bracket']
        va, vb = map(mpf, ref.get('bracket_function_values', ref['bracket_D']))
        slope = (vb - va) / (mpf(b) - mpf(a))
        out[name] = {'description': LABELS.get(name, name), 'kind': ref.get('kind', 'krein'),
                     'best_estimates': {p: runs[p]['best_estimate'] for p in conf['pipelines']}, **conf,
                     'verified_brackets': brackets,
                     'all_brackets_change_sign': all(v['sign_change'] for v in brackets.values()),
                     'transition_function_slope': mpmath.nstr(slope, 10)}
        same = conf['same_implementation_pairs']
        lines.append(f"| {LABELS.get(name, name)} | {out[name]['kind']} | {', '.join(conf['pipelines'])} | "
                     f"{conf['confirmed_estimate'] or 'single implementation only'} | "
                     f"{conf['confirmed_digits']} ({conf['confirming_pair']}) | "
                     f"{max(same.values()) if same else '-'} | {mpmath.nstr((mpf(b) - mpf(a)) / 2, 2)} | "
                     f"{mpmath.nstr(va, 4)}, {mpmath.nstr(vb, 4)} | {mpmath.nstr(slope, 8)} |")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'stability_transitions.json').write_text(json.dumps(out, indent=1))
    (args.output_dir / 'stability_transitions.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
