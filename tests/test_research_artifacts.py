"""Acceptance tests run after finalize_results.py; do not mistake partial sweeps for completion."""
import hashlib
import json
from decimal import Decimal as D, localcontext
from pathlib import Path
import csv


def test_original_checksums_preserved():
    root=Path('threebody_codex_project')
    for line in (root/'SHA256SUMS.txt').read_text().splitlines():
        expected,name=line.split(maxsplit=1)
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==expected


def test_complete_double_sweep_and_boundary_checks():
    rows=list(csv.DictReader(open('results/pancake_mass_sweep.csv')))
    assert len(rows)==23
    differences=list(csv.DictReader(open('results/baseline_numerical_differences.csv')))
    assert len(differences)==23
    for row in differences:
        assert abs(float(row['D']))<1e-6
        assert abs(float(row['T']))<1e-8
    convergence=json.load(open('results/boundary_convergence_double.json'))
    assert len(convergence)==24
    assert all(r['integration_success'] and r['closure']<1e-8 for r in convergence)
    for mu in [.5990,.5991,.6178,.6179]:
        assert len({r['D']>0 for r in convergence if r['mu']==mu})==1


def test_complete_high_precision_sweeps_and_actual_closure():
    for filename,count in [('high_precision_orbits.json',8),('high_precision_orbits_85.json',8),
                            ('boundary_orbits_65.json',4),('boundary_orbits_85.json',4)]:
        data=json.load(open('results/'+filename));assert len(data['orbits'])==count
        with localcontext() as ctx:
            ctx.prec=data['precision_digits']
            for orbit in data['orbits']:
                assert D(orbit['max_closure'])<D('1e-28')
                assert D(orbit['sampled_max_energy_drift'])<D('1e-28')
                assert D(orbit['floquet']['full_symplectic_defect'])<D('1e-28')
                assert D(orbit['floquet']['reduced_symplectic_defect'])<D('1e-28')
                assert D(orbit['floquet']['reciprocal_defect'])<D('1e-28')
                assert D(orbit['floquet']['determinant_defect'])<D('1e-28')
                assert all(D(m['eigen_residual'])<D('1e-28') for m in orbit['floquet']['multipliers'])


def test_precision_comparison_and_opposite_krein_brackets():
    convergence=json.load(open('results/high_precision_convergence.json'))['comparisons']
    assert len(convergence)==12
    assert all(D(r['D_difference'])<D('1e-28') and D(r['quartet_max_difference'])<D('1e-28') for r in convergence)
    report=json.load(open('results/stability_boundaries.json'))
    assert len(report['boundaries'])==2
    for b in report['boundaries']:
        assert D(b['endpoint_D'][0])*D(b['endpoint_D'][1])<0
        assert D(b['width'])<=D('2e-8')
        assert set(b['positive_imaginary_unit_modes_krein_signatures'])=={-1,1}


def test_full_period_time_reversal_identity():
    data=json.load(open('results/time_reversal_verification.json'))
    assert D(data['monodromy_max_identity_defect'])<D('1e-28')
    assert D(data['reversed_orbit']['max_closure'])<D('1e-28')
