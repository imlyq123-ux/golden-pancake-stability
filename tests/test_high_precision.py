from decimal import Decimal as D, localcontext
from pathlib import Path
import json
import subprocess
import numpy as np
from high_precision import C, nullspace, mm, transpose, inverse, spectrum, reduce_floquet, state
from threebody_expand import integrate


def test_decimal_complex_square_roots():
    with localcontext() as ctx:
        ctx.prec=65
        for value in [C(-4),C(4),C(3,4),C(-3,-4),C('1e20','1e-40')]:
            root=value.sqrt()
            assert abs(root*root-value)<D('1e-60')*(1+abs(value))


def test_nullspace_and_pivoted_inverse():
    with localcontext() as ctx:
        ctx.prec=65
        a=[[D(1),D(2),D(3)],[D(2),D(1),D(1)]]
        q,_,free=nullspace(a)
        assert len(free)==1
        assert max(abs(x) for row in mm(a,q) for x in row)<D('1e-60')
        b=[[D(2),D(3)],[D(1),D(4)]]
        product=mm(b,inverse(b))
        assert max(abs(product[i][j]-D(i==j)) for i in range(2) for j in range(2))<D('1e-60')


def test_symplectic_quartic_on_known_rotations():
    with localcontext() as ctx:
        ctx.prec=65
        a=[[D('.6'),D('.8'),D(0),D(0)],[-D('.8'),D('.6'),D(0),D(0)],
           [D(0),D(0),D('-.6'),D('.8')],[D(0),D(0),-D('.8'),D('-.6')]]
        omega=[[D(0),D(1),D(0),D(0)],[-D(1),D(0),D(0),D(0)],
               [D(0),D(0),D(0),D(1)],[D(0),D(0),-D(1),D(0)]]
        dat=spectrum(a,omega)
        assert abs(D(dat['D'])-D('5.76'))<D('1e-60')
        assert all(D(m['eigen_residual'])<D('1e-60') for m in dat['multipliers'])
        assert all(abs(D(m['modulus'])-1)<D('1e-60') for m in dat['multipliers'])


def test_reduced_quartet_matches_full_double_spectrum():
    p=['.151481541079005834','.140303019570684087','8.846608333332379']
    sol=integrate(.613,*map(float,p),variational=True)
    assert sol.success
    phi=sol.y[12:,-1].reshape(12,12)
    with localcontext() as ctx:
        ctx.prec=65
        reduced,omega,checks,_=reduce_floquet([[D(str(v)) for v in row] for row in phi],state('.613',p),'.613')
        data=spectrum(reduced,omega)
    full=np.linalg.eigvals(phi);full=full[np.argsort(abs(full-1))[-4:]]
    small=np.array([complex(float(m['lambda'][0]),float(m['lambda'][1])) for m in data['multipliers']])
    assert max(min(abs(v-full)) for v in small)<1e-6
    assert float(checks['quotient_reconstruction_defect'])<1e-7


def test_compiled_taylor_short_segment_against_dop853(tmp_path):
    output=tmp_path/'short.csv'
    command=['tools/taylor_mpfr','65','32','1e-40','.613','.151481541079005834','.140303019570684087','.1','12',str(output),'1']
    subprocess.run(command,check=True,capture_output=True)
    lines=[s for s in output.read_text().splitlines() if not s.startswith('#')]
    last=np.array(list(map(float,lines[-1].split(','))))
    sol=integrate(.613,.151481541079005834,.140303019570684087,.1,variational=True)
    # MP output stores columns contiguously, SciPy stores row-major.
    np.testing.assert_allclose(last[1:13],sol.y[:12,-1],atol=2e-13,rtol=2e-13)
    mp_phi=last[13:].reshape(12,12).T
    np.testing.assert_allclose(mp_phi,sol.y[12:,-1].reshape(12,12),atol=2e-12,rtol=2e-12)


def test_saved_high_precision_actual_closure_and_reduction():
    path=Path('results/high_precision_orbits.json')
    assert path.exists(), 'Run high_precision.py before the numerical artifact tests'
    data=json.loads(path.read_text())
    with localcontext() as ctx:
        ctx.prec=data['precision_digits']
        for orbit in data['orbits']:
            assert len(orbit['return_vector'])==12
            residual=[D(b)-D(a) for a,b in zip(orbit['initial_state'],orbit['final_state'])]
            assert residual==list(map(D,orbit['return_vector']))
            assert max(map(abs,residual))<D('1e-28')
            assert D(orbit['floquet']['full_symplectic_defect'])<D('1e-28')
            assert D(orbit['floquet']['quotient_reconstruction_defect'])<D('1e-28')
            assert max(D(m['eigen_residual']) for m in orbit['floquet']['multipliers'])<D('1e-28')
