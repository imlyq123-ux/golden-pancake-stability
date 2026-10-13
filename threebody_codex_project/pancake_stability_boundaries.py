"""Refine Floquet discriminant zeros for the specific Golden Pancake isola branch."""
from threebody_expand import correction, integrate, xinit, stats
from scipy.optimize import brentq
import numpy as np,json,time
p0=np.array([.15148154107900583396835825988003861711916800727,.14030301957068408729274415253719345009218867155,8.846608333332378559519969249427806733522940068])
cache={}; t0=time.time()
for path in ["/mnt/data/pancake_continuation_A.json","/mnt/data/pancake_continuation_B.json","/mnt/data/pancake_continuation_C.json"]:
    for item in json.load(open(path)):
        cache[item["mu"]]={"mu":item["mu"],"param":np.array([item["u"],item["v"],item["T"]]),"D":None,"growth":None,"residual":item["closure"],"multipliers":None}


def get(mu,rtol=2e-13,atol=2e-15):
    # use neighboring corrected seed, speed up convergence
    seed=cache[min(cache,key=lambda m:abs(m-mu))]['param'] if cache else p0
    p,res,sol=correction(mu,seed,rtol=rtol,atol=atol)
    f=integrate(mu,*p,variational=True,rtol=rtol,atol=atol)
    x=f.y[12:,-1].reshape(12,12)
    zz=np.linalg.eigvals(x)
    nt=zz[np.argsort(abs(zz-1))[::-1][:4]]
    S=sum(nt);Q=sum(nt**2)
    discr=float((2*Q+8-S*S).real)
    growth=float(max(abs(nt)))
    o={'mu':mu,'param':p,'D':discr,'growth':growth,'residual':float(np.max(abs(f.y[:12,-1]-xinit(mu,*p[:2])))),
       'multipliers':[[float(z.real),float(z.imag)] for z in nt]}
    cache[mu]=o
    print(f'mu={mu:.12f} D={discr:+.12e} growth={growth:.9f} residual={o["residual"]:.2e} took {time.time()-t0:.1f}s',flush=True)
    return discr

if __name__=='__main__':
   for a,b in [(.597,.601),(.6175,.618)]:
       fa=get(a); fb=get(b)
       root=brentq(get,a,b,xtol=3e-11,rtol=1e-12)
       print('boundary',root)
   with open('/mnt/data/pancake_boundaries.json','w') as ff:
       json.dump([{**v,'param':v['param'].tolist()} for v in sorted(cache.values(),key=lambda v:v['mu'])],ff,indent=2)
