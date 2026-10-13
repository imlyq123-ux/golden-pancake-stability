"""Independently correct and reintegrate the supplied 23 samples; preserve originals."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from threebody_expand import correction, stats, integrate, energy, xinit


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer=csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=Path('results'))
    parser.add_argument('--figure-dir',type=Path,default=Path('figures'))
    args=parser.parse_args()
    args.output_dir.mkdir(parents=True,exist_ok=True); args.figure_dir.mkdir(parents=True,exist_ok=True)
    baseline=list(csv.DictReader((Path(__file__).parent/'threebody_codex_project/pancake_mass_sweep.csv').open()))
    rows=[]; raw=[]; diffs=[]
    for source in baseline:
        mu=float(source['mu']); seed=[float(source[k]) for k in ['u','v','T']]
        p,res,sol=correction(mu,seed)
        dat=stats(mu,p,hi_variational=True)
        dat.update(correction_residual=res,correction=sol.correction_info)
        z=np.array(dat['mult_real'])+1j*np.array(dat['mult_imag'])
        dat.update(D=float((2*(z*z).sum()+8-z.sum()**2).real),max_lambda=float(max(abs(z))))
        row={k:dat[k] for k in ['mu','u','v','T','E','S','min_distance','t_nearest','pair_nearest','closure','D','max_lambda']}
        row['classification']='spectrally stable candidate' if dat['D']>0 and dat['max_lambda']<1+1e-5 else 'unstable candidate'
        row['nontrivial_multipliers']='; '.join(f'{w.real:.12g}{w.imag:+.12g}j' for w in z)
        raw.append(dat); rows.append(row)
        diffs.append({'mu':mu,**{k:float(dat[k])-float(source[k]) for k in ['u','v','T','E','S','min_distance','closure','D','max_lambda']}})
        (args.output_dir/'mass_sweep_raw_double.json').write_text(json.dumps(raw,indent=2))
        print(f"mu={mu:.6f} D={dat['D']:+.8e} rho={dat['max_lambda']:.9f} closure={dat['closure']:.2e}",flush=True)
    write_csv(args.output_dir/'pancake_mass_sweep.csv',rows)
    write_csv(args.output_dir/'baseline_numerical_differences.csv',diffs)
    for key,name,ylabel in [('max_lambda','pancake_stability_sweep','Largest nontrivial multiplier modulus'),('D','pancake_discriminant_sweep','Quartet discriminant D'),('min_distance','pancake_nearest_approach','Minimum pair separation')]:
        fig,ax=plt.subplots(figsize=(8,4.5))
        ax.plot([r['mu'] for r in rows],[r[key] for r in rows],'o-',ms=4)
        ax.set(xlabel='Middle mass mu',ylabel=ylabel,title='Recomputed binary64 data — provisional')
        if key!='min_distance': ax.axhline(1 if key=='max_lambda' else 0,color='gray',ls='--')
        # Measured endpoint sign brackets, rather than unsupported hard-coded roots.
        for lo,hi in [(.5990,.5991),(.6178,.6179)]: ax.axvspan(lo,hi,color='orange',alpha=.3)
        ax.grid(alpha=.2);fig.tight_layout();fig.savefig(args.figure_dir/(name+'.png'),dpi=180);plt.close(fig)
    conservation=[]
    for mu in [.597,.605,.613,.618,.625]:
        r=next(r for r in rows if r['mu']==mu)
        sol=integrate(mu,r['u'],r['v'],r['T'],dense=True,rtol=9e-14,atol=1e-15)
        if not sol.success: raise RuntimeError(sol.message)
        states=sol.sol(np.linspace(0,r['T'],501)).T
        e0=energy(mu,states[0]);m=np.array([1,mu,1])
        angular=[sum(m*(s[:6].reshape(3,2)[:,0]*s[6:].reshape(3,2)[:,1]-s[:6].reshape(3,2)[:,1]*s[6:].reshape(3,2)[:,0])) for s in states]
        conservation.append({'mu':mu,'sampled_max_energy_drift':max(abs(energy(mu,s)-e0) for s in states),
            'sampled_max_abs_angular_momentum':max(abs(np.array(angular))),
            'reintegrated_closure':max(abs(states[-1]-xinit(mu,r['u'],r['v'])))})
    write_csv(args.output_dir/'pancake_conservation_tests.csv',conservation)


if __name__=='__main__': main()
