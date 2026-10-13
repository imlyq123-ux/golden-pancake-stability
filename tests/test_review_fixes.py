"""Regression tests for the review findings (failure propagation, failure bookkeeping, independence)."""
import os
import stat
import subprocess
from pathlib import Path

import numpy as np

import asymmetric_double as ad
from asymmetric_summary import row_complete, summarize
from summarize_transitions import common_digits, confirmation

ROOT = Path(__file__).resolve().parent.parent


# ---- 2. parallel worker failures must fail the Make target -------------------------------------

def test_run_parallel_waits_on_every_pid():
    ok = subprocess.run(['bash', 'tools/run_parallel.sh', 'exit 0', 'sleep 0.2'], cwd=ROOT, capture_output=True, text=True)
    assert ok.returncode == 0
    bad = subprocess.run(['bash', 'tools/run_parallel.sh', 'sleep 0.2', 'exit 3', 'exit 0'], cwd=ROOT,
                         capture_output=True, text=True)
    assert bad.returncode != 0 and 'FAILED (exit 3)' in bad.stderr


def fake_python(tmp_path, failing_pattern):
    calls = tmp_path / 'calls.log'
    script = tmp_path / 'fakepython'
    script.write_text('#!/usr/bin/env bash\n'
                      f'echo "$*" >> {calls}\n'
                      f'case "$*" in {failing_pattern}) exit 3;; esac\n'
                      'exit 0\n')
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    return script, calls


def run_make_refine(tmp_path, failing_pattern):
    script, calls = fake_python(tmp_path, failing_pattern)
    env = dict(os.environ)
    done = subprocess.run(['make', '-s', 'refine', f'PYTHON={script}', f'REFINE_DIR={tmp_path}/refine',
                           f'SUMMARY_DIR={tmp_path}/summary', f'SCAN_DIR={tmp_path}/scan'],
                          cwd=ROOT, capture_output=True, text=True, env=env)
    lines = calls.read_text().splitlines() if calls.exists() else []
    return done, lines


def test_make_refine_fails_and_skips_summary_when_one_worker_fails(tmp_path):
    done, calls = run_make_refine(tmp_path, '*refine_boundaries.py*upper*gbs80*')
    assert done.returncode != 0
    assert not any('summarize_transitions.py' in c for c in calls)
    assert not any('K0_lowwindow_top' in c for c in calls)          # second batch never started
    assert sum('refine_boundaries.py' in c for c in calls) == 8      # all first-batch workers ran


def test_make_refine_fails_when_every_worker_fails(tmp_path):
    done, calls = run_make_refine(tmp_path, '*refine_boundaries.py*')
    assert done.returncode != 0 and not any('summarize_transitions.py' in c for c in calls)


def test_make_refine_succeeds_and_summarizes_when_all_workers_succeed(tmp_path):
    done, calls = run_make_refine(tmp_path, 'NEVER_MATCHES_ANYTHING')
    assert done.returncode == 0, done.stderr
    summaries = [c for c in calls if 'summarize_transitions.py' in c]
    assert len(summaries) == 1 and f'{tmp_path}/refine' in summaries[0]


# ---- 3. the asymmetric scan must keep failures and coverage ------------------------------------

def synthetic_evaluate(fails):
    """D(mu) = 0.5 - 8000 (mu - 0.608)^2: zeros near 0.6001 and 0.6159; `fails(mu)` -> Newton failure."""
    def evaluate(z, mu, delta):
        if fails(mu):
            raise RuntimeError('synthetic Newton failure')
        D = 0.5 - 8000.0 * (mu - 0.608) ** 2
        lam = np.exp(1j * np.array([2.0, -2.0, 1.2, -1.2])) * (1.0 if D > 0 else 1.05)
        info = {'closure': 1e-12, 'rmin': 0.01, 'newton_iterations': 1, 'reconstruction_defect': 0.0}
        return np.array(z, dtype=float), D, lam, info
    return evaluate


def test_scan_mu_records_failures_and_marks_row_incomplete(monkeypatch):
    monkeypatch.setattr(ad, 'evaluate', synthetic_evaluate(lambda mu: mu > 0.6125))
    row = ad.scan_mu((0.003, np.zeros(6), 0.608, 0.600, 0.620, 0.001))
    assert row['failures'] and row['failures'][0]['stage'] == 'grid' and abs(row['failures'][0]['mu'] - 0.613) < 1e-12
    assert row['coverage']['complete'] is False and 0.62 in row['coverage']['missing_points']
    assert not row_complete(row)
    assert [c['stable_side'] for c in row['crossings']] == ['above']   # upper crossing lies in the failed range
    monkeypatch.setattr(ad, 'evaluate', synthetic_evaluate(lambda mu: False))
    full = ad.scan_mu((0.003, np.zeros(6), 0.608, 0.600, 0.620, 0.001))
    assert full['coverage']['complete'] and row_complete(full)
    assert [c['stable_side'] for c in full['crossings']] == ['above', 'below']


def test_summary_never_reads_incomplete_rows_as_instability(monkeypatch):
    monkeypatch.setattr(ad, 'evaluate', synthetic_evaluate(lambda mu: False))
    complete = ad.scan_mu((0.0, np.zeros(6), 0.608, 0.600, 0.620, 0.001))
    # scan centred in the unstable part; the positive-D region fails, leaving only D < 0 points
    monkeypatch.setattr(ad, 'evaluate', synthetic_evaluate(lambda mu: mu < 0.6165))
    broken = ad.scan_mu((0.0005, np.zeros(6), 0.618, 0.600, 0.620, 0.001))
    assert broken['grid'] and max(g['D'] for g in broken['grid']) < 0
    data = {'status': 'synthetic', 'rows': [complete, broken], 'deltas_requested': [0.0, 0.0005]}
    out = summarize(data, [0.6001, 0.6159])
    assert out['window_closing'] is None                       # not misread as the window closing
    assert out['closing_sign_changes_not_interpreted_due_to_incomplete_rows'] == [[0.0, 0.0005]]
    assert out['incomplete_rows'][0]['delta'] == 0.0005 and out['all_rows_complete'] is False


# ---- 4. confirmed digits must come from independent implementations ----------------------------

def test_confirmation_ignores_taylor_taylor_agreement():
    runs = {'taylor70': '0.59904378420889985122072886046206841927541372566',
            'taylor90': '0.59904378420889985122072886046206841927541372566178',
            'gbs80': '0.5981234'}
    out = confirmation(runs)
    assert out['confirmed_digits'] == 2 and out['confirming_pair'].endswith('gbs80')
    assert out['same_implementation_pairs']['taylor70~taylor90'] >= 40


def test_common_digits_excludes_possibly_rounded_last_digit():
    assert common_digits('0.5800944142884689965488971809648796402745', '0.58009441428846899654889718096487964027444868') == 39
    assert common_digits('0.123', '0.124') == 2
    assert common_digits('0.12', '1.12') == 0
