from pancake_stability_boundaries import get,cache
from scipy.optimize import brentq
import json,numpy as np
for k in [0.6177,0.61785]:print('endpoint',k,get(k),flush=True)
root=brentq(get,0.6177,0.61785,xtol=1e-8)
print('UPPER_APPROX',repr(root))
with open('/mnt/data/pancake_upper_boundary.json','w') as f:json.dump([{**v,'param':v['param'].tolist()} for v in sorted(cache.values(),key=lambda v:v['mu']) if v['D'] is not None],f,indent=2)
