"""Compare the independent GBS/mpmath records with the Taylor/Decimal records and the atlas data.

Inputs:  results/independent_gbs_orbits.json, results/{high_precision_orbits,high_precision_orbits_85,
         boundary_orbits_65,boundary_orbits_85}.json, results/literature/atlas_reference_data.json
Outputs: results/independent_crosscheck.json and results/independent_crosscheck.md
"""
import json
from pathlib import Path

import mpmath
from mpmath import mp, mpf, mpc

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
mp.dps = 90


def load(name):
    return json.loads((RESULTS / name).read_text())


def e(x, n=3):
    return mpmath.nstr(mpf(x), n)


def taylor_records():
    out = {}
    for name in ['high_precision_orbits.json', 'boundary_orbits_65.json', 'high_precision_orbits_85.json', 'boundary_orbits_85.json']:
        data = load(name)
        for r in data['orbits']:
            out.setdefault(mpf(r['mu']), {})[data['precision_digits']] = (name, r)
    return out


def compare_taylor(g, name, t):
    fl = g['floquet']
    tf = t['floquet']
    gmodes = [(mpc(*map(mpf, m['lambda'])), m['krein_signature']) for m in fl['modes']]
    tmodes = [(mpc(*map(mpf, m['lambda'])), m['krein_signature']) for m in tf['multipliers']]
    dist, krein_ok = mpf(0), True
    for lam, sig in gmodes:
        other = min(tmodes, key=lambda z: abs(z[0] - lam))
        dist = max(dist, abs(other[0] - lam))
        krein_ok &= (sig == other[1])
    phi_g = [[mpf(x) for x in row] for row in g['full_monodromy']]
    phi_t = [[mpf(x) for x in row] for row in tf['full_monodromy']]
    return {'taylor_file': name, 'taylor_closure': t['max_closure'][:12], 'gbs_closure': g['max_closure'],
            'abs_du': e(abs(mpf(g['u']) - mpf(t['u']))), 'abs_dv': e(abs(mpf(g['v']) - mpf(t['v']))),
            'abs_dT': e(abs(mpf(g['T']) - mpf(t['T']))),
            'abs_dD': e(abs(mpf(fl['D_from_traces']) - mpf(tf['D']))),
            'D_gbs': mpmath.nstr(mpf(fl['D_from_traces']), 20),
            'quartet_max_difference': e(dist), 'krein_signs_agree': bool(krein_ok),
            'monodromy_max_difference': e(max(abs(a - b) for ra, rb in zip(phi_g, phi_t) for a, b in zip(ra, rb))),
            'monodromy_max_entry': e(max(abs(a) for row in phi_g for a in row), 6)}


def main():
    gbs = {mpf(r['mu']): r for r in load('independent_gbs_orbits.json')['orbits']}
    tay = taylor_records()
    rows = []
    for mu, g in sorted(gbs.items()):
        for digits, (name, t) in sorted(tay.get(mu, {}).items()):
            rows.append({'mu': g['mu'], 'taylor_digits': digits, **compare_taylor(g, name, t)})
    ref = load('literature/atlas_reference_data.json')
    draft = ref['draft_note']
    gp = gbs[mpf(draft['table1']['mu'])]
    modes = [mpc(*map(mpf, m['lambda'])) for m in gp['floquet']['modes']]
    draft_modes = [mpc(*map(mpf, z)) for z in draft['multipliers']]
    draft_cmp = {'mu': draft['table1']['mu'],
                 **{f'abs_d{k}': e(abs(mpf(gp[k]) - mpf(draft['table1'][k]))) for k in ['u', 'v', 'T']},
                 'our_u': gp['u'][:50], 'our_v': gp['v'][:50], 'our_T': gp['T'][:50],
                 'abs_dE': e(abs(mpf(gp['E']) - mpf(draft['E']))), 'our_E': gp['E'][:25],
                 'our_T_star': mpmath.nstr(mpf(gp['S']), 20), 'abs_dT_star': e(abs(mpf(gp['S']) - mpf(draft['T_star']))),
                 'multiplier_max_difference': e(max(min(abs(z - w) for w in modes) for z in draft_modes)),
                 'our_multipliers': [[mpmath.nstr(z.real, 22), mpmath.nstr(z.imag, 22)] for z in modes],
                 'our_D': mpmath.nstr(mpf(gp['floquet']['D_from_traces']), 20),
                 'our_spectral_radius_minus_1': e(mpf(gp['floquet']['spectral_radius_nontrivial']) - 1),
                 'closure': gp['max_closure']}
    cat = []
    decoded = ref['orbit_files']['decoded']
    for entry in ref['catalogue']['pancake_entries']:
        mu = mpf(repr(entry['masses'][1]))
        g = next((r for m, r in gbs.items() if abs(m - mu) < mpf('1e-15')), None)
        if g is None:
            continue
        ic = decoded.get(entry['key'])
        rec = {'key': entry['key'], 'atlas_mu': repr(entry['masses'][1]), 'our_mu': g['mu'],
               'atlas_stab': entry['stab'], 'our_rho': mpmath.nstr(mpf(g['floquet']['spectral_radius_nontrivial']), 8),
               'atlas_E': entry['E'], 'our_E': mpmath.nstr(mpf(g['E']), 10),
               'atlas_rmin': entry['rmin'], 'our_sampled_rmin': mpmath.nstr(mpf(g['sampled_min_distance_upper_estimate']), 6),
               'atlas_T0': entry['T0'], 'our_T': mpmath.nstr(mpf(g['T']), 18)}
        if ic:
            pub = [mpf(x) for x in ic['ic_published']]
            rec.update({'abs_du_vs_published_f64': e(abs(mpf(g['u']) - pub[6])),
                        'abs_dv_vs_published_f64': e(abs(mpf(g['v']) - pub[7])),
                        'abs_dT_vs_file_f64': e(abs(mpf(g['T']) - mpf(ic['T'])))})
        cat.append(rec)
    out = {'description': 'Independent GBS/mpmath pipeline versus Taylor/Decimal records and atlas data',
           'taylor_comparisons': rows, 'atlas_draft_comparison': draft_cmp, 'atlas_catalogue_comparison': cat}
    (RESULTS / 'independent_crosscheck.json').write_text(json.dumps(out, indent=1))
    lines = ['# Independent cross-check (GBS + mpmath vs Taylor + Decimal vs atlas)', '',
             '| mu | Taylor digits | GBS closure | |du| | |dv| | |dT| | |dD| | quartet diff | Krein agree | max |dPhi| |',
             '|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        lines.append(f"| {r['mu']} | {r['taylor_digits']} | {r['gbs_closure']} | {r['abs_du']} | {r['abs_dv']} | {r['abs_dT']} | "
                     f"{r['abs_dD']} | {r['quartet_max_difference']} | {r['krein_signs_agree']} | {r['monodromy_max_difference']} |")
    lines += ['', '## Atlas draft (42-digit Table 1, 19-digit multipliers)', '',
              '```', json.dumps(draft_cmp, indent=1), '```', '', '## Atlas catalogue entries on this branch', '',
              '```', json.dumps(cat, indent=1), '```']
    (RESULTS / 'independent_crosscheck.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
