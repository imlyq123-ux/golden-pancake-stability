"""Tests for the independent GBS/mpmath pipeline and the binary64 isola scan."""
import json
import subprocess
from pathlib import Path

import mpmath
import numpy as np
from mpmath import mp, mpf

import independent_mp as im
from threebody_expand import integrate

SEED = ['0.15148154107900583396835825988003861711916800727',
        '0.14030301957068408729274415253719345009218867155',
        '8.846608333332378559519969249427806733522940068']


def last_row(path, columns=None):
    rows = [line for line in Path(path).read_text().splitlines() if not line.startswith('#')]
    return rows[-1].split(',')


def test_gbs_short_segment_against_dop853(tmp_path):
    out = tmp_path / 'gbs.csv'
    subprocess.run(['tools/gbs_mpfr', '40', '1e-30', '10', '.613', '.151481541079005834', '.140303019570684087',
                    '.1', '12', str(out), '1000'], check=True, capture_output=True)
    last = np.array(list(map(float, last_row(out))))
    sol = integrate(.613, .151481541079005834, .140303019570684087, .1, variational=True, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(last[1:13], sol.y[:12, -1], atol=2e-13, rtol=2e-13)
    np.testing.assert_allclose(last[13:].reshape(12, 12).T, sol.y[12:, -1].reshape(12, 12), atol=2e-12, rtol=2e-12)


def test_gbs_and_taylor_full_period_return_vectors_agree(tmp_path):
    """Two different integration algorithms reproduce the seed's 5.14e-40 return vector."""
    mp.dps = 60
    gbs, tay = tmp_path / 'g.csv', tmp_path / 't.csv'
    subprocess.run(['tools/gbs_mpfr', '60', '1e-52', '16', '0.613', *SEED, '0', str(gbs), '1000000'], check=True, capture_output=True)
    subprocess.run(['tools/taylor_mpfr', '60', '40', '1e-52', '0.613', *SEED, '0', str(tay), '1000000'], check=True, capture_output=True)
    x0 = im.initial_state('0.613', SEED[0], SEED[1])
    rg = [mpf(a) - b for a, b in zip(last_row(gbs)[1:13], x0)]
    rt = [mpf(a) - b for a, b in zip(last_row(tay)[1:13], x0)]
    assert max(abs(a - b) for a, b in zip(rg, rt)) < mpf('1e-46')
    assert mpf('5.1e-40') < max(map(abs, rg)) < mpf('5.2e-40')


def test_trace_floquet_and_full_space_krein_on_synthetic_matrix():
    """Omega-symplectic conjugate of two rotations of opposite Krein sign plus a unipotent block."""
    mp.dps = 50
    mu = mpf('0.613')
    n = 6
    J = mpmath.zeros(12, 12)
    for i in range(n):
        J[i, n + i], J[n + i, i] = 1, -1
    th1, th2 = mpf('2.1'), mpf('1.3')
    C = mpmath.eye(12)
    for k, th, sign in [(0, th1, 1), (1, th2, -1)]:  # H = +-(q^2+p^2)/2 in the (q_k, p_k) plane
        c, s = mpmath.cos(th), sign * mpmath.sin(th)
        C[k, k], C[k, n + k], C[n + k, k], C[n + k, n + k] = c, s, -s, c
    for k in range(2, 6):  # unipotent (Jordan) trivial block: q_k += p_k
        C[k, n + k] = 1
    rng = np.random.default_rng(1)
    H = mpmath.matrix(rng.normal(size=(12, 12)) * .15)
    H = (H + H.T) / 2
    R = mpmath.expm(J * H)
    A = R * C * mpmath.inverse(R)                      # J-symplectic
    B = mpmath.diag([1] * 6 + [1, 1, mu, mu, 1, 1])
    phi = mpmath.inverse(B) * A * B                     # Omega-symplectic, Omega = B^T J B
    out = im.floquet(phi, mu)
    expected_D = (2 * mpmath.cos(th1) - 2 * mpmath.cos(th2)) ** 2
    assert abs(mpf(out['D_from_traces']) - expected_D) < mpf('1e-35')
    assert mpf(out['full_symplectic_defect']) < mpf('1e-40')
    upper = [m for m in out['modes'] if mpf(m['lambda'][1]) > 0]
    assert len(upper) == 2 and sorted(m['krein_signature'] for m in upper) == [-1, 1]
    for m in upper:  # the conjugation preserves the sign attached to each rotation angle
        angle = mpmath.atan2(mpf(m['lambda'][1]), mpf(m['lambda'][0]))
        assert min(abs(angle - th1), abs(angle - th2)) < mpf('1e-30')


def test_mass_sensitivity_column_matches_finite_differences():
    from isola_scan_double import shoot
    z = np.array([.151481541079005834, .140303019570684087, 8.846608333332379, .613])
    _, jacobian, _, _, _ = shoot(z, 1e-13, 1e-15)
    h = 1e-6
    e = np.array([0, 0, 0, h])
    fd = (shoot(z + e, 1e-13, 1e-15)[0] - shoot(z - e, 1e-13, 1e-15)[0]) / (2 * h)
    assert np.max(abs(jacobian[:, 3] - fd)) < 1e-6 * np.max(abs(fd))
    sv = np.linalg.svd(jacobian / np.array([1, 1, .1, 1]), compute_uv=False)
    assert sv[-1] < 1e-7 * sv[-2]          # one-dimensional branch tangent
