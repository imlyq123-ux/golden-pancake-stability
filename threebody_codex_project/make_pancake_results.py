import json, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from threebody_expand import integrate, energy, xinit
root=Path('/mnt/data')
raw=[]
for fname in ['pancake_continuation_A.json','pancake_continuation_B.json','pancake_continuation_C.json','pancake_transition_coarse.json']:
    raw.extend(json.load(open(root/fname)))
rows={}
for r in raw:
    mu=r['mu']
    ev=np.array([complex(*z) for z in r['multipliers_all']])
    nt=ev[np.argsort(abs(ev-1))[::-1][:4]]
    s=nt.sum(); p=(nt*nt).sum()
    rr={k:r[k] for k in ['mu','u','v','T','E','S','min_distance','t_nearest','pair_nearest','closure']}
    rr.update({'D':float((2*p+8-s*s).real),'max_lambda':max(abs(nt)),
               'classification':'stable candidate' if ((2*p+8-s*s).real>0 and max(abs(nt))<1+1e-5) else 'unstable',
               'nontrivial_multipliers':'; '.join(f'{z.real:.9g}{z.imag:+.9g}j' for z in nt)})
    rows[mu]=rr

df=pd.DataFrame([v for k,v in sorted(rows.items())]);df.to_csv(root/'pancake_mass_sweep.csv',index=False)
low=0.5990437843;high=0.6178442682

fig,ax=plt.subplots(figsize=(9,5.4))
ax.axvspan(low,high,alpha=.12,label='Numerically stable window')
ax.plot(df.mu,df.max_lambda,'o-',lw=1.7,ms=5,label='Largest nontrivial |λ|')
ax.axhline(1,ls='--',lw=.9)
ax.axvline(low,ls=':',lw=1.1); ax.axvline(high,ls=':',lw=1.1)
ax.set(xlabel=r'Middle mass $\mu$  ($m_1=m_3=1$)',ylabel=r'Maximum nontrivial Floquet modulus',title='Floquet stability along the Golden Pancake branch')
ax.grid(alpha=.20);ax.legend(loc='upper left');fig.tight_layout()
fig.savefig(root/'pancake_stability_sweep.png',dpi=230)
plt.close(fig)

fig,ax=plt.subplots(figsize=(9,5.4))
ax.axvspan(low,high,alpha=.12,label='Numerically stable window')
ax.plot(df.mu,df.D,'o-',lw=1.7,ms=5,label='Floquet discriminant')
ax.axhline(0,ls='--',lw=.9)
ax.axvline(low,ls=':',lw=1.1);ax.axvline(high,ls=':',lw=1.1)
ax.set(xlabel=r'Middle mass $\mu$',ylabel=r'$D = 2\sum_{j=1}^4\lambda_j^2+8-(\sum_{j=1}^4\lambda_j)^2$',title='Two candidate Hamiltonian-Hopf transitions')
ax.grid(alpha=.20);ax.legend(loc='lower center');fig.tight_layout()
fig.savefig(root/'pancake_discriminant_sweep.png',dpi=230)
plt.close(fig)

fig,ax=plt.subplots(figsize=(9,5.4))
ax.axvspan(low,high,alpha=.12,label='Numerically stable window')
ax.plot(df.mu,df.min_distance,'o-',lw=1.7,ms=5,label='Minimum pair separation')
ax.axvline(low,ls=':',lw=1.1);ax.axvline(high,ls=':',lw=1.1)
ax.set(xlabel=r'Middle mass $\mu$',ylabel='Minimum pair separation (over one period)',title='Close approach along the orbit branch')
ax.grid(alpha=.20);ax.legend(loc='upper left');fig.tight_layout()
fig.savefig(root/'pancake_nearest_approach.png',dpi=230)
plt.close(fig)

conservation=[]
for mu in [0.597,0.605,0.613,0.618,0.625]:
  r=rows[mu];sol=integrate(mu,r['u'],r['v'],r['T'],dense=True,rtol=9e-14,atol=1e-15)
  ts=np.linspace(0,r['T'],501); states=sol.sol(ts).T
  eng=np.array([energy(mu,x) for x in states]);mom=[]
  for x in states:
    q=x[:6].reshape(3,2);v=x[6:].reshape(3,2)
    mom.append(sum(m*(q[i,0]*v[i,1]-q[i,1]*v[i,0]) for i,m in enumerate([1,mu,1])))
  c={'mu':mu,'sampled_max_energy_drift':float(max(abs(eng-eng[0]))), 'sampled_max_abs_angular_momentum':float(max(abs(np.array(mom)))), 'reintegrated_closure':float(max(abs(states[-1]-xinit(mu,r['u'],r['v']))))}
  conservation.append(c)
pd.DataFrame(conservation).to_csv(root/'pancake_conservation_tests.csv',index=False)
print('Wrote',len(df),'mass points and',len(conservation),'conservation cases')
print(df[['mu','T','min_distance','D','max_lambda','closure']].to_string(index=False))
print(pd.DataFrame(conservation).to_string(index=False))
