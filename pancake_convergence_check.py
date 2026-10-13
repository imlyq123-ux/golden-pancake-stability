"""Reproduce the original failure, then compare both sides of both transitions."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from threebody_expand import correction, legacy_correction, integrate, xinit, integration_metadata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-dir', type=Path, default=Path('results'))
    args = parser.parse_args(); args.output_dir.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parent/'threebody_codex_project'
    baseline = sum((json.loads((root/name).read_text()) for name in ['pancake_continuation_C.json', 'pancake_continuation_B.json']), [])
    closest = min(baseline, key=lambda r: abs(r['mu']-.599))
    guess = [closest[k] for k in ['u', 'v', 'T']]
    diagnosis = {'historical_seed_mass': closest['mu'], 'original_guess': guess,
                 'historical_exception': 'Required step size is less than spacing between numbers.'}
    print('Original unbounded corrector, mu=.599, nearest seed:', closest['mu'], flush=True)
    try:
        p, residual, sol = legacy_correction(.599, guess, rtol=9e-14, atol=1e-15)
        diagnosis.update(reproduced_exception=False, returned_param=p.tolist(), residual=residual)
    except RuntimeError as exc:
        diagnosis.update(reproduced_exception=True, exception=str(exc))
        print('Reproduced:', exc, flush=True)
    (args.output_dir/'convergence_failure_diagnosis.json').write_text(json.dumps(diagnosis, indent=2))
    rows = []
    for mu in [.5990, .5991, .6178, .6179]:
        seed = min(baseline, key=lambda r: abs(r['mu']-mu))
        # Linear interpolation of two independent archived neighbors is a predictor,
        # never accepted as a periodic solution without correction and reintegration.
        ordered = sorted(baseline, key=lambda r: r['mu'])
        lo = max((r for r in ordered if r['mu'] <= mu), key=lambda r: r['mu'])
        hi = min((r for r in ordered if r['mu'] >= mu), key=lambda r: r['mu'])
        alpha = (mu-lo['mu'])/(hi['mu']-lo['mu'])
        guess = np.array([(1-alpha)*lo[k]+alpha*hi[k] for k in ['u','v','T']])
        p, _, correction_sol = correction(mu, guess, rtol=9e-14, atol=1e-15)
        print('mass', mu, 'corrected', p, flush=True)
        for rtol in [2e-12, 3e-13, 9e-14]:
            for max_step in [np.inf, .01]:
                started = time.perf_counter()
                sol = integrate(mu, *p, variational=True, rtol=rtol, atol=1e-15, max_step=max_step)
                if not sol.success: raise RuntimeError(sol.message)
                matrix = sol.y[12:, -1].reshape(12,12)
                z = np.linalg.eigvals(matrix); z = z[np.argsort(abs(z-1))[::-1][:4]]
                discr = (2*(z*z).sum()+8-z.sum()**2).real
                masses = np.diag(np.repeat([1,mu,1],2)); zero=np.zeros((6,6))
                omega=np.block([[zero,masses],[-masses,zero]])
                row = {'mu':mu, 'param':p.tolist(), 'rtol':rtol, 'atol':1e-15,
                    'max_step':None if np.isinf(max_step) else max_step,
                    'D':float(discr), 'spectral_radius':float(max(abs(z))),
                    'closure':float(max(abs(sol.y[:12,-1]-xinit(mu,*p[:2])))),
                    'symplectic_defect':float(np.max(abs(matrix.T@omega@matrix-omega))),
                    'runtime_seconds':time.perf_counter()-started,
                    'correction':correction_sol.correction_info, **integration_metadata(sol)}
                rows.append(row)
                print(f"rtol={rtol:g} step={max_step:g} D={discr:+.10e} closure={row['closure']:.3e}",flush=True)
                (args.output_dir/'boundary_convergence_double.json').write_text(json.dumps(rows,indent=2))
    assert len(rows)==24 and all(r['integration_success'] for r in rows)
    for mu in [.5990,.5991,.6178,.6179]:
        signs={np.sign(r['D']) for r in rows if r['mu']==mu}
        assert len(signs)==1 and 0 not in signs
    print('PASS: 24 successful integrations; consistent D signs at four endpoints.', flush=True)


if __name__=='__main__': main()
