"""Portable double-precision bracket exploration; not a precision certification."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from threebody_expand import correction, stats
from high_precision import seed_for


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output-dir',type=Path,default=Path('results'))
    parser.add_argument('--upper-only',action='store_true')
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True);rows=[]
    def get(mu):
        p,res,sol=correction(mu,list(map(float,seed_for(str(mu)))))
        dat=stats(mu,p,hi_variational=True)
        z=np.array(dat['mult_real'])+1j*np.array(dat['mult_imag'])
        dat['D']=float((2*(z*z).sum()+8-z.sum()**2).real)
        dat['correction']=sol.correction_info;rows.append(dat)
        print(mu,dat['D'],flush=True)
        return dat['D']
    roots=[]
    for lo,hi in ([(.6178,.6179)] if args.upper_only else [(.5990,.5991),(.6178,.6179)]):
        roots.append(float(brentq(get,lo,hi,xtol=2e-9)))
    (args.output_dir/'pancake_boundaries.json').write_text(json.dumps({'status':'binary64 estimates; provisional','roots':roots,'evaluations':rows},indent=2))


if __name__=='__main__':main()
