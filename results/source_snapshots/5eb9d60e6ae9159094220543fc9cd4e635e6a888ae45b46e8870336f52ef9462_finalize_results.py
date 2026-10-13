"""Reanalyze saved monodromies, compare runs, and generate reports/figures."""
import csv
import hashlib
import itertools
import json
import xml.etree.ElementTree as ET
from decimal import Decimal as D, localcontext
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from high_precision import C, mm, spectrum, reduce_floquet, state, invariants

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'; FIGURES=ROOT/'figures'


def save_json(path,data):path.write_text(json.dumps(data,indent=2))


def reanalyze(path):
    data=json.loads(path.read_text())
    if 'analysis_revision' not in data:
        save_json(path.with_name(path.stem+'_initial_analysis.json'),data)
    with localcontext() as ctx:
        ctx.prec=data['precision_digits']
        for orbit in data['orbits']:
            phi=[list(map(D,row)) for row in orbit['floquet']['full_monodromy']]
            initial=list(map(D,orbit['initial_state']))
            reduced,omega,checks,basis=reduce_floquet(phi,initial,orbit['mu'])
            orbit['floquet'].update(spectrum(reduced,omega));orbit['floquet'].update(checks)
    data['analysis_revision']='Actual characteristic polynomial roots; reciprocal/unit-circle symmetry not imposed.'
    data['analysis_source_sha256']=hashlib.sha256((ROOT/'high_precision.py').read_bytes()).hexdigest()
    save_json(path,data);return data


def compare(low,high):
    lower={D(r['mu']):r for r in low['orbits']}; rows=[]
    with localcontext() as ctx:
        ctx.prec=high['precision_digits']
        for r in high['orbits']:
            if D(r['mu']) not in lower:continue
            base=lower[D(r['mu'])]
            row={'mu':r['mu'],'lower_digits':low['precision_digits'],'higher_digits':high['precision_digits'],
                'parameter_differences':{k:str(abs(D(r[k])-D(base[k]))) for k in ['u','v','T']},
                'D_difference':str(abs(D(r['floquet']['D'])-D(base['floquet']['D']))),
                'rho_difference':str(abs(D(r['floquet']['spectral_radius'])-D(base['floquet']['spectral_radius']))),
                'closure_lower':base['max_closure'],'closure_higher':r['max_closure']}
            a=[[D(v) for v in line] for line in base['floquet']['full_monodromy']]
            b=[[D(v) for v in line] for line in r['floquet']['full_monodromy']]
            row['monodromy_max_difference']=str(max(abs(x-y) for ra,rb in zip(a,b) for x,y in zip(ra,rb)))
            z=[C(*m['lambda']) for m in base['floquet']['multipliers']]
            w=[C(*m['lambda']) for m in r['floquet']['multipliers']]
            row['quartet_max_difference']=str(min(max(abs(z[i]-w[p[i]]) for i in range(4)) for p in itertools.permutations(range(4))))
            rows.append(row)
    return rows


def track(orbits):
    tracks=[];previous=None
    with localcontext() as ctx:
        ctx.prec=85
        for orbit in sorted(orbits,key=lambda r:D(r['mu'])):
            mult=orbit['floquet']['multipliers'];values=[C(*m['lambda']) for m in mult]
            if previous is None:perm=tuple(range(4));cost=D(0)
            else:
                perm=min(itertools.permutations(range(4)),key=lambda p:sum(abs(previous[i]-values[p[i]])**2 for i in range(4)))
                cost=sum(abs(previous[i]-values[perm[i]])**2 for i in range(4))
            previous=[values[i] for i in perm]
            tracks.append({'mu':orbit['mu'],'tracked_quartet':[values[i].json() for i in perm],
                'minimum_assignment_cost':str(cost),'krein_signatures':[mult[i]['krein_signature'] for i in perm]})
    return {'method':'Decimal eigenvalue-distance minimum over all 24 permutations; near exact collisions individual labels are ambiguous.',
            'samples':tracks}


def boundary_report(low,high,comparisons):
    records=[]
    hi={D(r['mu']):r for r in high['orbits']};lo={D(r['mu']):r for r in low['orbits']}
    differences={D(r['mu']):r for r in comparisons}
    with localcontext() as ctx:
        ctx.prec=85
        for name,a,b in [('lower','0.59904378','0.59904379'),('upper','0.61784426','0.61784428')]:
            ra,rb=hi[D(a)],hi[D(b)];da,db=D(ra['floquet']['D']),D(rb['floquet']['D'])
            assert da*db<0
            for value in [D(a),D(b)]:
                assert D(lo[value]['floquet']['D'])*D(hi[value]['floquet']['D'])>0
                assert D(hi[value]['max_closure'])<D('1e-28')
                assert abs(D(hi[value]['floquet']['D']))>D(differences[value]['D_difference'])*1000
            stable=ra if da>0 else rb
            modes=[m for m in stable['floquet']['multipliers'] if D(m['lambda'][1])>0]
            assert len(modes)==2 and {m['krein_signature'] for m in modes}=={-1,1}
            records.append({'name':name,'bracket':[a,b],'width':str(D(b)-D(a)),
                'endpoint_D':[str(da),str(db)],'endpoint_spectral_radius':[ra['floquet']['spectral_radius'],rb['floquet']['spectral_radius']],
                'precision_runs':[65,85],
                'endpoint_D_changes':[differences[D(a)]['D_difference'],differences[D(b)]['D_difference']],
                'positive_imaginary_unit_modes_krein_signatures':[m['krein_signature'] for m in modes],
                'stable_side_mass':stable['mu'],
                'classification':'candidate spectral transition; opposite Krein signs and off-circle quartet observed',
                'exact_collision_multiplicity':'unresolved; no integration at the exact zero of D'})
    return {'status':'empirically converged numerical sign brackets; not rigorous interval enclosures',
        'boundaries':records,'limitations':['No independent multiprecision implementation or interval-arithmetic proof.',
        'No exact-transition algebraic/geometric multiplicity proof.','No nonlinear stability or novelty claim.']}


def figures(orbits,boundaries):
    ordered=sorted(orbits,key=lambda r:D(r['mu']))
    for key,name,label in [('D','floquet_discriminant_high_precision','Reduced-quartet discriminant D'),
                           ('spectral_radius','floquet_radius_high_precision','Spectral radius')]:
        fig,ax=plt.subplots(figsize=(8,4.5))
        ax.plot([float(r['mu']) for r in ordered],[float(r['floquet'][key]) for r in ordered],'o-',ms=4)
        ax.axhline(0 if key=='D' else 1,color='gray',ls='--')
        for r in boundaries['boundaries']:
            a,b=map(float,r['bracket']);ax.axvspan(a,b,color='orange',alpha=.6);ax.axvline((a+b)/2,color='orange',ls=':')
        ax.set(xlabel='Middle mass mu',ylabel=label,title='MPFR / Decimal Floquet verification — provisional')
        ax.grid(alpha=.2);fig.tight_layout();fig.savefig(FIGURES/(name+'.png'),dpi=180);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10,8))
    for column,b in enumerate(boundaries['boundaries']):
        ax=axes[0,column]; zoom=axes[1,column]; local_points=[]
        for value in b['bracket']:
            r=next(r for r in ordered if D(r['mu'])==D(value))
            ax.scatter([float(m['lambda'][0]) for m in r['floquet']['multipliers']],
                       [float(m['lambda'][1]) for m in r['floquet']['multipliers']],label=f'mu={value}')
            positive=[m for m in r['floquet']['multipliers'] if D(m['lambda'][1])>0]
            points=np.array([[float(m['lambda'][0]),float(m['lambda'][1])] for m in positive])
            local_points.extend(points.tolist());zoom.scatter(points[:,0],points[:,1],label=f'mu={value}')
        theta=np.linspace(0,2*np.pi,500);ax.plot(np.cos(theta),np.sin(theta),color='gray',lw=.8)
        ax.set(title=b['name']+' candidate transition',xlabel='Re(lambda)',ylabel='Im(lambda)',aspect='equal')
        ax.legend(fontsize=8);ax.grid(alpha=.2)
        points=np.array(local_points);center=points.mean(axis=0);angle=np.arctan2(center[1],center[0])
        theta=np.linspace(angle-.003,angle+.003,500);zoom.plot(np.cos(theta),np.sin(theta),color='gray',lw=.8)
        margin=max(np.ptp(points[:,0]),np.ptp(points[:,1]))*.8
        zoom.set(xlim=(center[0]-margin,center[0]+margin),ylim=(center[1]-margin,center[1]+margin),
                 title='Positive-imaginary modes, enlarged',xlabel='Re(lambda)',ylabel='Im(lambda)',aspect='equal')
        zoom.ticklabel_format(axis='both',style='plain',useOffset=False);zoom.tick_params(axis='x',labelrotation=25);zoom.grid(alpha=.2)
    fig.tight_layout();fig.savefig(FIGURES/'floquet_complex_plane_high_precision.png',dpi=200);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for ax,b in zip(axes,boundaries['boundaries']):
        a,c=map(D,b['bracket']);center=(a+c)/2
        ax.plot([float((D(v)-center)/D('1e-8')) for v in b['bracket']],list(map(float,b['endpoint_D'])),'o-')
        ax.axhline(0,color='gray',ls='--');ax.axvspan(float((a-center)/D('1e-8')),float((c-center)/D('1e-8')),alpha=.15,color='orange')
        ax.set(xlabel=f'(mu - {center}) / 1e-8',ylabel='D',title=b['name']+' measured sign bracket')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(FIGURES/'stability_boundary_uncertainty.png',dpi=200);plt.close(fig)


def orbit_figure(data):
    orbit=next(r for r in data['orbits'] if D(r['mu'])==D('.613'))
    source=ROOT/orbit['verification_integration']['raw_file']
    with source.open() as stream:
        samples=[line.split(',')[:13] for line in stream if not line.startswith('#')]
    target=RESULTS/'golden_pancake_verified_trajectory.csv'
    with target.open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['t']+data['state_order']);writer.writerows(samples)
    a=np.array(samples,dtype=float)
    fig,ax=plt.subplots(figsize=(7,6))
    for i in range(3):ax.plot(a[:,1+2*i],a[:,2+2*i],lw=1,label=f'body {i+1}')
    ax.set(aspect='equal',xlabel='x',ylabel='y',title='Verified full-state orbit, mu=0.613');ax.legend();ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(FIGURES/'golden_pancake_verified_orbit.png',dpi=180);plt.close(fig)
    with localcontext() as ctx:
        ctx.prec=data['precision_digits'];initial=invariants(list(map(D,orbit['initial_state'])),orbit['mu']);rows=[]
        for raw in samples:
            inv=invariants(list(map(D,raw[1:])),orbit['mu'])
            rows.append([float(D(raw[0])),float(abs(inv[0]-initial[0])),float(abs(inv[1])),float(inv[4])])
    a=np.array(rows)
    fig,axes=plt.subplots(2,1,figsize=(8,6),sharex=True)
    axes[0].semilogy(a[:,0],np.maximum(a[:,1],1e-85),label='Absolute energy drift')
    axes[0].semilogy(a[:,0],np.maximum(a[:,2],1e-85),label='Absolute angular momentum');axes[0].legend()
    axes[1].plot(a[:,0],a[:,3]);axes[1].set(xlabel='Time',ylabel='Nearest pair separation')
    axes[0].set(title='Diagnostics at MPFR accepted-step nodes',ylabel='Absolute error')
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(FIGURES/'high_precision_conservation.png',dpi=180);plt.close(fig)


def write_status():
    datasets=[json.loads((RESULTS/name).read_text()) for name in ['high_precision_orbits.json','high_precision_orbits_85.json','boundary_orbits_65.json','boundary_orbits_85.json']]
    orbits=[r for data in datasets for r in data['orbits']]
    convergence=json.loads((RESULTS/'high_precision_convergence.json').read_text())['comparisons']
    boundaries=json.loads((RESULTS/'stability_boundaries.json').read_text())['boundaries']
    reversal=json.loads((RESULTS/'time_reversal_verification.json').read_text())
    worst_closure=max(D(r['max_closure']) for r in orbits)
    worst_energy=max(D(r['sampled_max_energy_drift']) for r in orbits)
    worst_symplectic=max(D(r['floquet']['full_symplectic_defect']) for r in orbits)
    worst_d=max(D(r['D_difference']) for r in convergence)
    worst_quartet=max(D(r['quartet_max_difference']) for r in convergence)
    table='\n'.join(f"| {b['name']} | [{b['bracket'][0]}, {b['bracket'][1]}] | {b['width']} | {D(b['endpoint_D'][0]):+.5E}, {D(b['endpoint_D'][1]):+.5E} |" for b in boundaries)
    tests='Final acceptance tests have not yet run; see the existing test logs.'
    if (RESULTS/'tests_all.xml').exists():
        tree=ET.parse(RESULTS/'tests_all.xml');suites=tree.getroot().findall('testsuite')
        count=sum(int(s.attrib['tests']) for s in suites);failed=sum(int(s.attrib['failures'])+int(s.attrib['errors']) for s in suites)
        tests=f'{count} final tests run; {failed} failures/errors. See `tests_all.xml` and `tests_all.log`.'
    text=f'''# Research status — numerical verification, not publication-ready

Phase A reproduction and the planned symmetric-mass high-precision orbit /
Floquet checks are complete for the saved samples. Phases D/E remain incomplete
or deferred. No existence proof, nonlinear stability result, exact bifurcation
classification, or novelty claim is made.

## Phase A

- The original supplied files match their SHA256 manifest and remain unchanged.
- The six-point smoke test completed successfully. The 23-point binary64 sweep
  reproduced all compared numerical columns exactly in this environment, using
  the original archived seeds and fresh integrations. This is not an independent
  seed search; the already small residuals did not require parameter changes.
- All three original plot types were regenerated from recalculated data.
  PNG hashes intentionally differ because annotations now describe provisional
  candidates and measured sign brackets instead of unsupported precise roots.
- The original working note's claim of completed tolerance checks was unsupported
  by its archived exception. The replacement check completed **24** successful
  integrations: four endpoint masses, three tolerances, and two maximum steps,
  with consistent signs at each endpoint. Detailed solver status/nfev are saved.
- The exact archived nearest-seed exception did not recur. A larger 0.613→0.599
  jump reproduced it after an off-branch near-collision trial reached separation
  ≈2.12e-7. Both successful nearest seeds and failed trial history are preserved.
  Damping, trial rejection and small mass steps address this demonstrated failure
  mechanism. The exact historical trigger remains uncertain.
- A clean isolated Python environment installed successfully from cached wheels.
  Exact Python package versions, compiler and MPFR version are recorded. Numba
  and pandas were unavailable offline; NumPy and stdlib CSV are documented
  substitutes, and MPFR/Decimal provide the arbitrary-precision backend.

## Phases B/C

Eight mass samples and four narrow-bracket endpoints were computed at each
of **65 and 85 decimal digits** (24 corrected/verified orbit records). Correction
and verification use different Taylor orders and tolerances; every verification
is a separate full-period integration. The seed at mu=0.613 retains the supplied
decimal parameters because its measured residual already meets the target.

The largest actual full-state componentwise return across these runs is
**{worst_closure:.5E}**, below 1e-28. All twelve components and complete terminal
states are printed in the machine-readable JSON records. Maximum sampled
energy drift is {worst_energy:.5E}. Maximum full mass-weighted symplectic defect
is {worst_symplectic:.5E}. Sampled nearest distances are upper estimates of
the true minimum and do not certify absence of collisions between nodes.

Across 12 matched mass samples, the largest change in D between precision/order
runs is {worst_d:.5E}; the largest optimally matched quartet change is
{worst_quartet:.5E}. These are empirical convergence measurements, not rigorous
error enclosures or evidence that every printed digit is significant.

The full monodromy is reduced explicitly by COM/momentum constraints and a
time/energy, rotation/angular-momentum quotient. Actual reduced polynomial
roots are computed in Decimal complex arithmetic, with eigenvector residuals,
reciprocal/conjugate symmetry, determinant and symplectic checks. Saved mode
tracking matches nearby eigenvalues; individual labels at exact collisions
remain ambiguous.

| Candidate transition | Numerical mass bracket | Width | Endpoint D values |
|---|---|---|---|
{table}

The positive-imaginary unit-circle modes on the stable side have **opposite
Krein signatures** at both brackets; off-circle complex quartets appear on
the other side. These support candidate spectral transitions. Exact-transition
algebraic/geometric multiplicities are unresolved. Numerical endpoint signs
do not constitute interval-arithmetic enclosures or certify all intermediate
masses as spectrally stable.

An independently initialized, time-reversed seed also closes the full state.
Its full monodromy agrees with R*inverse(Phi)*R to
{D(reversal['monodromy_max_identity_defect']):.5E}, with R=diag(I6,-I6). This is
an additional reversibility check, not an independent MP implementation.

## Tests and saved artifacts

{tests}

See `high_precision_orbits.json`, `high_precision_orbits_85.json`,
`boundary_orbits_65.json`, `boundary_orbits_85.json`,
`floquet_sweep_high_precision.csv`, `stability_boundaries.json`,
`high_precision_convergence.json`, and `time_reversal_verification.json`.
Raw accepted-step states, velocities and variational columns are in `mp_raw/`.
The verified full-state trajectory and orbit/conservation/Floquet/boundary
figures are saved. Source snapshots and file hashes preserve provenance.
The initial reciprocal-quartic analysis is retained separately; finalized
spectra refine roots on the actual polynomial instead of imposing symmetry.

## Remaining work

1. Cross-check with an independent arbitrary-precision implementation and/or
   interval arithmetic, particularly the monodromy and transition enclosures.
   Two Taylor settings share an implementation; binary64 DOP853 is the existing
   independent method but cannot verify 28-digit accuracy by itself.
2. Resolve the exact-collision multiplicities and refine Krein/Jordan analysis
   before upgrading the candidate bifurcation label.
3. Obtain the complete original atlas draft and the 2023 article/supplements.
   The draft PDF could not be fetched; only the atlas pages and publisher
   preview were inspected. The literature/novelty audit is incomplete.
4. Only after that audit consider a manuscript or asymmetric-mass continuation
   with a newly derived full shooting/phase setup. No manuscript is prepared.

Reproduce with the root `README.md`; numerical conventions and solver limits
are documented in `methodology.md`. Literature attribution/access limits are
in `literature_audit.md`.
'''
    (RESULTS/'research_status_first_session.md').write_text(text)  # current report: research_status.md


def main():
    files=['high_precision_orbits.json','high_precision_orbits_85.json','boundary_orbits_65.json','boundary_orbits_85.json']
    sets=[reanalyze(RESULTS/name) for name in files]
    assert [len(s['orbits']) for s in sets]==[8,8,4,4], 'Finish the sweeps before finalizing'
    comparisons=compare(sets[0],sets[1])+compare(sets[2],sets[3])
    save_json(RESULTS/'high_precision_convergence.json',{'method':'Independent shooting/reintegration at 65/85 digits and different Taylor orders/tolerances. Same MP implementation.', 'comparisons':comparisons})
    boundaries=boundary_report(sets[2],sets[3],comparisons)
    save_json(RESULTS/'stability_boundaries.json',boundaries)
    all_orbits=sets[1]['orbits']+sets[3]['orbits']
    save_json(RESULTS/'floquet_mode_tracking.json',track(all_orbits))
    rows=[]
    for dataset in sets:
        for orbit in dataset['orbits']:
            f=orbit['floquet']; row={k:orbit[k] for k in ['mu','u','v','T','E','S','max_closure','sampled_min_distance_upper_bound']}
            row.update(digits=dataset['precision_digits'],order=orbit['verification_integration']['order'],
                tolerance=orbit['verification_integration']['tolerance'],D=f['D'],spectral_radius=f['spectral_radius'],
                full_symplectic_defect=f['full_symplectic_defect'],reduced_symplectic_defect=f['reduced_symplectic_defect'],
                multipliers=json.dumps([m['lambda'] for m in f['multipliers']]),
                krein_signatures=json.dumps([m['krein_signature'] for m in f['multipliers']]))
            rows.append(row)
    with (RESULTS/'floquet_sweep_high_precision.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    figures(all_orbits,boundaries);orbit_figure(sets[1])
    write_status()
    for name in files:
        print(name,len(json.loads((RESULTS/name).read_text())['orbits']))
    for b in boundaries['boundaries']:print(b['name'],b['bracket'],b['endpoint_D'])
    print('Maximum observed D change:',max(D(r['D_difference']) for r in comparisons))


if __name__=='__main__':main()
