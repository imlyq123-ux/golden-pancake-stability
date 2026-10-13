from threebody_expand import correction,integrate,xinit,energy
import json
rows=sum((json.load(open("/mnt/data/"+f)) for f in ["pancake_continuation_C.json","pancake_continuation_B.json"]),[])
def closest(mu):
  r=min(rows,key=lambda r:abs(mu-r["mu"]));return [r["u"],r["v"],r["T"]]
import numpy as np
for mu,guess in [(mu,closest(mu)) for mu in [0.5990,0.5991,0.6178,0.6179]]:
  print('mass',mu,flush=True)
  p,_,_=correction(mu,guess,rtol=9e-14,atol=1e-15)
  print('p',p,flush=True)
  for rtol in [2e-12,3e-13,9e-14]:
    sol=integrate(mu,*p,variational=True,rtol=rtol,atol=1e-15)
    M=sol.y[12:,-1].reshape(12,12)
    z=np.linalg.eigvals(M);z=z[np.argsort(abs(z-1))[::-1][:4]]
    s=z.sum();d=(2*(z*z).sum()+8-s*s).real
    cl=max(abs(sol.y[:12,-1]-xinit(mu,*p[:2])))
    print(' ',rtol,'D',f'{d:.10e}','growth',f'{max(abs(z)):.10f}','closure',f'{cl:.3e}',flush=True)
