import numpy as np
import pytest
from threebody_expand import f, jac, xinit, energy, integrate, shooting_system

SEED = np.array([.151481541079005834, .140303019570684087, 8.846608333332379])


def test_force_balance_and_pair_antisymmetry():
    m = np.array([1., .613, 1.])
    x = xinit(.613, *SEED[:2]); x[:6] += np.array([.1, .2, .2, -.3, .3, .1])
    acc = f(0, x, m)[6:].reshape(3, 2)
    assert np.max(abs(m @ acc)) < 1e-14
    q = x[:6].reshape(3, 2)
    for i in range(3):
        for j in range(i+1, 3):
            d = q[j]-q[i]
            fij = m[i]*m[j]*d/np.linalg.norm(d)**3
            fji = m[j]*m[i]*(-d)/np.linalg.norm(-d)**3
            np.testing.assert_allclose(fij+fji, 0, atol=1e-15)
            isolated=m.copy();isolated[3-i-j]=0
            measured=f(0,x,isolated)[6:].reshape(3,2)*isolated[:,None]
            np.testing.assert_allclose(measured[i],fij,atol=1e-14)
            np.testing.assert_allclose(measured[j],fji,atol=1e-14)


def test_analytic_jacobian_and_mass_symplectic_generator():
    x = xinit(.613, *SEED[:2]); x[1::2][:3] += [.2, -.1, .3]
    m = np.array([1., .613, 1.]); h = 1e-5
    fd = np.column_stack([(f(0, x+h*np.eye(12)[k], m)-f(0, x-h*np.eye(12)[k], m))/(2*h) for k in range(12)])
    a = jac(x, m)
    np.testing.assert_allclose(a, fd, atol=3e-9, rtol=3e-9)
    mass = np.diag(np.repeat(m, 2)); zero = np.zeros((6, 6))
    omega = np.block([[zero, mass], [-mass, zero]])
    assert np.max(abs(a.T@omega+omega@a)) < 1e-14


def test_shooting_jacobian_full_period_finite_difference():
    rr, jj, sol = shooting_system(.613, SEED)
    assert sol.success
    h = 1e-7
    fd = []
    for k in range(3):
        d = np.eye(3)[k]*h
        plus = integrate(.613, *(SEED+d)); minus = integrate(.613, *(SEED-d))
        assert plus.success and minus.success
        rp = plus.y[:12, -1]-xinit(.613, *(SEED+d)[:2])
        rm = minus.y[:12, -1]-xinit(.613, *(SEED-d)[:2])
        fd.append((rp-rm)/(2*h))
    np.testing.assert_allclose(jj, np.array(fd).T, atol=5e-4, rtol=5e-4)


def test_section_constraints_and_full_labeled_closure():
    x = xinit(.613, *SEED[:2]); m = np.array([1., .613, 1.])
    q, v = x[:6].reshape(3, 2), x[6:].reshape(3, 2)
    assert np.max(abs(m@q)) < 1e-15
    assert np.max(abs(m@v)) < 1e-15
    assert abs(np.sum(m*(q[:, 0]*v[:, 1]-q[:, 1]*v[:, 0]))) < 1e-15
    sol = integrate(.613, *SEED, dense=True)
    assert sol.success
    assert np.max(abs(sol.y[:, -1]-x)) < 5e-9
    states = sol.sol(np.linspace(0, SEED[2], 501)).T
    assert max(abs(energy(.613, s)-energy(.613, x)) for s in states) < 1e-9
    for s in states:
        q, v = s[:6].reshape(3, 2), s[6:].reshape(3, 2)
        assert abs(np.sum(m*(q[:, 0]*v[:, 1]-q[:, 1]*v[:, 0]))) < 1e-10
        assert np.max(abs(m@q)) < 1e-10
        assert np.max(abs(m@v)) < 1e-10
