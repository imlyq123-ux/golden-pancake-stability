"""Acceptance tests for the second-session result files; each assertion backs a documented claim."""
import json
from pathlib import Path

import pytest
from mpmath import mp, mpf

mp.dps = 100
R = Path('results')


def load(name):
    path = R / name
    if not path.exists():
        pytest.skip(f'{name} not generated yet')
    return json.loads(path.read_text())


def test_independent_gbs_orbits_close_and_match_taylor_and_atlas():
    data = load('independent_gbs_orbits.json')
    assert len(data['orbits']) == 15
    for o in data['orbits']:
        assert mpf(o['max_closure']) < mpf('1e-55')
        assert mpf(o['floquet']['full_symplectic_defect']) < mpf('1e-50')
        assert mpf(o['sampled_max_energy_drift']) < mpf('1e-55')
    cross = load('independent_crosscheck.json')
    for row in cross['taylor_comparisons']:
        assert row['krein_signs_agree']
        limit = mpf('1e-55') if row['taylor_digits'] == 85 and row['mu'] in ('0.59904378', '0.59904379', '0.61784426', '0.61784428') else mpf('1e-36')
        assert mpf(row['abs_dD']) < limit
    draft = cross['atlas_draft_comparison']
    assert mpf(draft['abs_du']) < mpf('1e-37') and mpf(draft['abs_dv']) < mpf('1e-37') and mpf(draft['abs_dT']) < mpf('1e-35')
    assert mpf(draft['multiplier_max_difference']) < mpf('1e-19')
    for entry in cross['atlas_catalogue_comparison']:
        assert abs(float(entry['our_rho']) - entry['atlas_stab']) < 6e-5


def test_refined_transitions_have_verified_sign_brackets_and_agree():
    data = load('stability_transitions.json')
    for name in ['lower', 'upper', 'K0_lowwindow_top', 'K3_nearfold']:
        t = data[name]
        assert t['all_brackets_change_sign']
        assert len(t['pipelines']) >= 2
        # confirmation must come from the independent implementations (Taylor vs GBS)
        assert t['confirming_pair'] in t['independent_pairs'] and t['confirmed_digits'] >= 39
        assert abs(mpf(t['transition_function_slope'])) > 1
    for name in ['lower', 'upper']:
        assert data[name]['confirmed_digits'] >= 59
    assert abs(mpf(data['lower']['best_estimates']['gbs80']) - mpf('0.5990437842088998512207288604620684192754')) < mpf('1e-40')
    assert abs(mpf(data['upper']['best_estimates']['gbs80']) - mpf('0.6178442682090631787273667106004179967021')) < mpf('1e-40')


def test_main_transitions_are_non_semisimple_krein_collisions():
    data = load('transition_analysis.json')
    for name in ['lower', 'upper']:
        s = data['summary'][name]
        for key in ['loglog_slope_stable_split', 'loglog_slope_unstable_rho_minus_1', 'loglog_slope_krein_form']:
            assert abs(s[key] - 0.5) < 0.005
        assert abs(s['loglog_slope_one_minus_overlap'] - 1.0) < 0.005
        assert all(sorted(pair) == [-1, 1] for pair in s['stable_side_krein_sign_pairs'])
        sv = [mpf(x) for x in s['singular_values_at_mu_star']]
        assert sv[0] < mpf('1e-55') and sv[1] > mpf('1e-3')


def test_half_period_reflection_symmetry():
    data = load('half_period_symmetry.json')
    assert len(data['results']) == 15
    for r in data['results']:
        assert mpf(r['symmetry_defect_max_abs']) < mpf('1e-60')
        assert mpf(r['N_squared_minus_M_max_abs']) < mpf('1e-55')
        assert mpf(r['pair_sum_identity_defect']) < mpf('1e-50')
    d_n = {r['mu']: mpf(r['D_N']) for r in data['results']}
    assert d_n['0.59904378'] < 0 < d_n['0.59904379'] and d_n['0.61784428'] < 0 < d_n['0.61784426']


def test_isola_scan_transitions_match_multiprecision_values():
    data = load('isola_scan_summary.json')
    w = data['AABAbbaabaBB (W) segment']
    krein = sorted(t['mu'] for t in w['transitions'] if t['type'] == 'krein')
    expected = [0.5800944142884690, 0.5990437842088999, 0.6178442682090632, 0.6341970669305546]
    assert len(krein) == 4
    assert all(abs(a - b) < 5e-5 for a, b in zip(krein, expected))   # linear interpolation between scan points
    assert w['max_heuristic_vs_quotient'] < 1e-6
    tangent = [t['mu'] for t in w['transitions'] if t['type'] == 'tangent']
    assert len(tangent) == 1 and abs(tangent[0] - 0.63460425) < 2e-5
    assert len(w['stable_intervals']) == 3


def test_dense_sweep_reproduces_archived_samples():
    data = load('dense_sweep_double.json')
    comp = data['comparison_with_archived_23_point_sweep']
    assert len(comp) == 23
    assert max(abs(c['u']) for c in comp) < 1e-11 and max(abs(c['T']) for c in comp) < 1e-9
    assert max(abs(c['D']) for c in comp) < 1e-6
    assert max(r['heuristic_vs_quotient_max_difference'] for r in data['records'] if 'D' in r) < 1e-6


def test_asymmetric_mass_scan():
    data = load('asymmetric/main_window_summary.json')
    assert data['all_rows_complete'] and not data['incomplete_rows']
    assert data['max_abs_theta'] < 1e-9
    assert 0.0085 < data['window_closing']['delta_linear_estimate'] < 0.009
    fits = data['quadratic_fits_mu_star_vs_delta2']
    assert abs(fits['mu_lower']['delta0_scan_error']) < 1e-9 and abs(fits['mu_upper']['delta0_scan_error']) < 1e-9
    assert 120 < fits['mu_lower']['c2_small_delta_richardson'] < 130
    assert -42 < fits['mu_upper']['c2_small_delta_richardson'] < -36
    mirror = data['mirror_symmetry_check']
    assert mirror is not None and mirror['max_difference'] < 1e-9
