"""MPFR Taylor shooting and Decimal symmetry-reduced Floquet analysis.

All inputs and saved MP quantities are decimal strings. Decimal precision is
never replaced by NumPy eigensolvers. The C++ executable is the integrator.
"""
import argparse
import csv
import hashlib
import json
import subprocess
import time
from decimal import Decimal as D, getcontext
from pathlib import Path

ROOT=Path(__file__).resolve().parent
ZERO=D(0); ONE=D(1)
SEED=['0.15148154107900583396835825988003861711916800727',
      '0.14030301957068408729274415253719345009218867155',
      '8.846608333332378559519969249427806733522940068']


def mm(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def transpose(a): return [list(row) for row in zip(*a)]
def eye(n): return [[D(i==j) for j in range(n)] for i in range(n)]
def sub(a,b): return [[x-y for x,y in zip(r,s)] for r,s in zip(a,b)]
def maxabs(a): return max(abs(x) for row in a for x in row)


def solve(a,b):
    """Small dense pivoted elimination, in the supplied scalar arithmetic."""
    n=len(a); c=[list(row)+[v] for row,v in zip(a,b)]
    for k in range(n):
        p=max(range(k,n),key=lambda i:abs(c[i][k]));c[k],c[p]=c[p],c[k]
        if abs(c[k][k])==0: raise ArithmeticError('Singular linear system')
        pivot=c[k][k];c[k]=[v/pivot for v in c[k]]
        for i in range(n):
            if i==k:continue
            factor=c[i][k];c[i]=[x-factor*y for x,y in zip(c[i],c[k])]
    return [row[-1] for row in c]


def inverse(a): return transpose([solve(a,col) for col in eye(len(a))])


def nullspace(a):
    rows=[list(r) for r in a]; pivots=[];n=len(a[0])
    threshold=D(10)**(-(getcontext().prec-15))
    for k in range(n):
        j=len(pivots)
        if j==len(rows):break
        p=max(range(j,len(rows)),key=lambda i:abs(rows[i][k]))
        if abs(rows[p][k])<threshold:continue
        rows[j],rows[p]=rows[p],rows[j];div=rows[j][k];rows[j]=[v/div for v in rows[j]]
        for i in range(len(rows)):
            if i==j:continue
            factor=rows[i][k];rows[i]=[x-factor*y for x,y in zip(rows[i],rows[j])]
        pivots.append(k)
    free=[k for k in range(n) if k not in pivots]
    columns=[]
    for k in free:
        col=[ZERO]*n;col[k]=ONE
        for j,p in enumerate(pivots):col[p]=-rows[j][k]
        columns.append(col)
    return transpose(columns),pivots,free


class C:
    """Decimal complex arithmetic for a reciprocal quartic and Krein forms."""
    def __init__(self,r=0,i=0):
        if isinstance(r,C):self.r,self.i=r.r,r.i
        else:self.r,self.i=D(r),D(i)
    def __add__(self,b):b=C(b);return C(self.r+b.r,self.i+b.i)
    __radd__=__add__
    def __neg__(self):return C(-self.r,-self.i)
    def __sub__(self,b):return self+-C(b)
    def __rsub__(self,b):return C(b)+-self
    def __mul__(self,b):b=C(b);return C(self.r*b.r-self.i*b.i,self.r*b.i+self.i*b.r)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=C(b);den=b.r*b.r+b.i*b.i
        return C((self.r*b.r+self.i*b.i)/den,(self.i*b.r-self.r*b.i)/den)
    def __rtruediv__(self,b):return C(b)/self
    def __abs__(self):return (self.r*self.r+self.i*self.i).sqrt()
    def conjugate(self):return C(self.r,-self.i)
    def sqrt(self):
        norm=abs(self)
        # Compute small component through the product identity to avoid cancellation.
        if self.r>=0:
            a=((norm+self.r)/2).sqrt();b=self.i/(2*a) if a else ZERO
        else:
            b=((norm-self.r)/2).sqrt();b=-b if self.i<0 else b
            a=self.i/(2*b) if b else ZERO
        return C(a,b)
    def json(self):return [str(self.r),str(self.i)]


def state(mu,p):
    u,v,_=map(D,p);mu=D(mu)
    return [D(-1),ZERO,ZERO,ZERO,ONE,ZERO,u,v,-2*u/mu,-2*v/mu,u,v]


def rhs(x,mu):
    m=[ONE,D(mu),ONE];a=[ZERO]*6
    for i,j in [(0,1),(0,2),(1,2)]:
        d=[x[2*j+k]-x[2*i+k] for k in range(2)]
        s=sum(v*v for v in d);g=ONE/(s*s.sqrt())
        for k in range(2):a[2*i+k]+=m[j]*d[k]*g;a[2*j+k]-=m[i]*d[k]*g
    return x[6:]+a


def invariants(x,mu):
    m=[ONE,D(mu),ONE];total=sum(m)
    e=sum(m[i]*(x[6+2*i]**2+x[7+2*i]**2)/2 for i in range(3))
    distances=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        r=sum((x[2*j+k]-x[2*i+k])**2 for k in range(2)).sqrt()
        e-=m[i]*m[j]/r;distances.append(r)
    angular=sum(m[i]*(x[2*i]*x[7+2*i]-x[2*i+1]*x[6+2*i]) for i in range(3))
    com=[sum(m[i]*x[2*i+k] for i in range(3))/total for k in range(2)]
    momentum=[sum(m[i]*x[6+2*i+k] for i in range(3)) for k in range(2)]
    return e,angular,com,momentum,min(distances)


def integrate(mu,p,digits,order,tol,columns,path):
    command=[str(ROOT/'tools/taylor_mpfr'),str(digits),str(order),str(tol),str(mu),
             *map(str,p),str(columns),str(path),'1']
    started=time.perf_counter()
    run=subprocess.run(command,capture_output=True,text=True,check=True)
    path.with_suffix('.log').write_text(run.stdout+run.stderr)
    lines=path.read_text().splitlines()
    if not lines[-1].startswith('# success=1'):raise RuntimeError('Incomplete integration')
    numeric=[list(map(D,line.split(','))) for line in lines if not line.startswith('#')]
    info={'digits':digits,'order':order,'tolerance':str(tol),'columns':columns,
          'runtime_seconds':time.perf_counter()-started,'success':True,
          'backend_header':lines[0],'completion':lines[-1],
          'raw_file':str(path.relative_to(ROOT)),
          'raw_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    return numeric,info


def correct(mu,p,digits,order,tol,out,prefix,target='1e-35'):
    p=list(map(D,p));history=[]
    for iteration in range(8):
        samples,info=integrate(mu,p,digits,order,tol,2,out/f'{prefix}_newton_{iteration}.csv')
        initial=state(mu,p);final=samples[-1][1:13];res=[b-a for a,b in zip(initial,final)]
        norm=max(map(abs,res));history.append({'iteration':iteration,'param':list(map(str,p)),
            'max_residual':str(norm),'return_vector':list(map(str,res)),'integration':info})
        print(f'{prefix} Newton {iteration}: residual={norm:.4E}',flush=True)
        if norm<D(target):return p,history
        dd=[[ZERO,ZERO] for _ in range(12)]
        dd[6]=[ONE,ZERO];dd[7]=[ZERO,ONE];dd[8]=[-2/D(mu),ZERO];dd[9]=[ZERO,-2/D(mu)]
        dd[10]=[ONE,ZERO];dd[11]=[ZERO,ONE]
        j=[[samples[-1][13+k]-dd[k][0],samples[-1][25+k]-dd[k][1],rhs(final,mu)[k]] for k in range(12)]
        jt=transpose(j);normal=mm(jt,j);right=[-sum(a*b for a,b in zip(row,res)) for row in jt]
        step=solve(normal,right)
        if max(abs(step[i])/D(['.002','.002','.04'][i]) for i in range(3))>1:
            raise RuntimeError('MP Newton step leaves local branch; improve continuation seed')
        p=[a+b for a,b in zip(p,step)]
    raise RuntimeError('MP Newton target not reached')


def reduce_floquet(phi,x,mu):
    masses=[ONE,ONE,D(mu),D(mu),ONE,ONE];omega=[[ZERO]*12 for _ in range(12)]
    for i,m in enumerate(masses):omega[i][6+i]=m;omega[6+i][i]=-m
    flow=rhs(x,mu)
    rotation=sum(([-x[i+1],x[i]] for i in range(0,12,2)),[])
    conserved=[]
    for start in [0,6]:
        for k in range(2):
            row=[ZERO]*12
            for i,m in enumerate([ONE,D(mu),ONE]):row[start+2*i+k]=m
            conserved.append(row)
    conserved+=mm([flow,rotation],omega)
    gauges=[flow,rotation];basis,pivots,free=nullspace(conserved+gauges)
    if len(free)!=4:raise RuntimeError('Reduction is not four-dimensional')
    neutral=transpose([flow,rotation]);gn_inverse=inverse(mm(gauges,neutral))
    projector=sub(eye(12),mm(mm(neutral,gn_inverse),gauges))
    image=mm(mm(projector,phi),basis)
    left_inverse=mm(inverse(mm(transpose(basis),basis)),transpose(basis))
    reduced=mm(left_inverse,image);omega_r=mm(mm(transpose(basis),omega),basis)
    checks={'full_symplectic_defect':str(maxabs(sub(mm(mm(transpose(phi),omega),phi),omega))),
        'reduced_symplectic_defect':str(maxabs(sub(mm(mm(transpose(reduced),omega_r),reduced),omega_r))),
        'constraint_basis_defect':str(maxabs(mm(conserved+gauges,basis))),
        'constraint_transport_defect':str(maxabs(mm(mm(conserved,phi),basis))),
        'quotient_reconstruction_defect':str(maxabs(sub(image,mm(basis,reduced)))),
        'flow_transport_defect':str(maxabs(sub(mm(phi,transpose([flow])),transpose([flow])))),
        'rotation_transport_defect':str(maxabs(sub(mm(phi,transpose([rotation])),transpose([rotation]))))}
    return reduced,omega_r,checks,basis


def eigenvector(a,lam):
    best=None
    for fixed in range(4):
        cols=[k for k in range(4) if k!=fixed]
        matrix=[[C(a[i][j])- (lam if i==j else C(0)) for j in range(4)] for i in range(4)]
        for excluded in range(4):
            rows=[k for k in range(4) if k!=excluded]
            try:
                small=[[matrix[i][j] for j in cols] for i in rows]
                values=solve(small,[-matrix[i][fixed] for i in rows])
            except (ArithmeticError,ZeroDivisionError):continue
            vec=[C(0)]*4;vec[fixed]=C(1)
            for j,v in zip(cols,values):vec[j]=v
            norm=sum((abs(z)**2 for z in vec),ZERO).sqrt();vec=[z/norm for z in vec]
            residual=max(abs(sum((matrix[i][j]*vec[j] for j in range(4)),C(0))) for i in range(4))
            if best is None or residual<best[0]:best=residual,vec
    if best is None:raise RuntimeError('Eigenvector solve failed')
    return best


def spectrum(reduced,omega):
    a=sum(reduced[i][i] for i in range(4));a2=mm(reduced,reduced)
    tr2=sum(a2[i][i] for i in range(4));b=(a*a-tr2)/2
    discr=2*tr2+8-a*a
    root=C(discr).sqrt();pair_sums=[(C(a)+root)/2,(C(a)-root)/2]
    vals=[]
    for s in pair_sums:
        delta=(s*s-4).sqrt();vals.extend([(s+delta)/2,(s-delta)/2])
    diagnostics=[]
    for lam in vals:
        error,vec=eigenvector(reduced,lam)
        form=sum((vec[i].conjugate()*omega[i][j]*vec[j] for i in range(4) for j in range(4)),C(0))
        # Convention: kappa = real(i * conjugate(z)^T Omega_reduced z).
        diagnostics.append({'lambda':lam.json(),'modulus':str(abs(lam)),
            'eigen_residual':str(error),'krein_form':str(-form.i),
            'krein_signature':int((-form.i>0)-(-form.i<0)) if abs(abs(lam)-1)<D('1e-25') else None,
            'normalized_eigenvector':[z.json() for z in vec]})
    reciprocal=max(min(abs(lam*other-1) for other in vals) for lam in vals)
    conjugate=max(min(abs(lam.conjugate()-other) for other in vals) for lam in vals)
    a3=mm(a2,reduced);tr3=sum(a3[i][i] for i in range(4))
    c3=(a**3-3*a*tr2+2*tr3)/6
    return {'D':str(discr),'spectral_radius':str(max(map(abs,vals))),
        'multipliers':diagnostics,'reciprocal_defect':str(reciprocal),'conjugate_defect':str(conjugate),
        'reciprocal_characteristic_coefficient_defect':str(abs(c3-a))}


def verify_orbit(mu,p,history,digits,order,tol,out,prefix,floquet=True):
    x0=state(mu,p);samples,info=integrate(mu,p,digits,order,tol,12 if floquet else 0,out/f'{prefix}_verification.csv')
    final=samples[-1][1:13];res=[b-a for a,b in zip(x0,final)]
    inv0=invariants(x0,mu);drift=[ZERO]*4;minimum=D('Infinity');t_near=ZERO
    for row in samples:
        inv=invariants(row[1:13],mu)
        ds=[abs(inv[0]-inv0[0]),abs(inv[1]),max(map(abs,inv[2])),max(map(abs,inv[3]))]
        drift=[max(a,b) for a,b in zip(drift,ds)]
        if inv[4]<minimum:minimum,t_near=inv[4],row[0]
    dat={'mu':str(mu),'u':str(p[0]),'v':str(p[1]),'T':str(p[2]),
        'E':str(inv0[0]),'S':str(p[2]*abs(inv0[0])*abs(inv0[0]).sqrt()),
        'initial_state':list(map(str,x0)),'final_state':list(map(str,final)),
        'return_vector':list(map(str,res)),'max_closure':str(max(map(abs,res))),
        'sampled_max_energy_drift':str(drift[0]),'sampled_max_abs_angular_momentum':str(drift[1]),
        'sampled_max_com':str(drift[2]),'sampled_max_momentum':str(drift[3]),
        'sampled_min_distance_upper_bound':str(minimum),'sampled_t_nearest':str(t_near),
        'minimum_distance_note':'Minimum over accepted step nodes, an upper estimate; not a certified continuous-time lower bound.',
        'correction_history':history,'verification_integration':info}
    if floquet:
        phi=transpose([samples[-1][13+12*j:25+12*j] for j in range(12)])
        reduced,omega,checks,basis=reduce_floquet(phi,x0,mu)
        dat['floquet']={**spectrum(reduced,omega),**checks,'full_monodromy':[[str(v) for v in row] for row in phi],
            'reduced_monodromy':[[str(v) for v in row] for row in reduced],
            'reduced_symplectic_form':[[str(v) for v in row] for row in omega],
            'reduction_basis':[[str(v) for v in row] for row in basis]}
        print(f"{prefix}: closure={dat['max_closure']} D={dat['floquet']['D']} rho={dat['floquet']['spectral_radius']}",flush=True)
    return dat


def seed_for(mu):
    if D(mu)==D('.613'):return SEED
    convergence=json.loads((ROOT/'results/boundary_convergence_double.json').read_text())
    candidates=[r for r in convergence if abs(D(str(r['mu']))-D(mu))<D('1e-15')]
    if candidates:return list(map(str,candidates[0]['param']))
    raw=json.loads((ROOT/'results/mass_sweep_raw_double.json').read_text())
    same=[r for r in raw if abs(D(str(r['mu']))-D(mu))<D('1e-15')]
    if same:return [str(same[0][k]) for k in ['u','v','T']]
    ordered=sorted(raw,key=lambda r:r['mu']);lower=max((r for r in ordered if D(str(r['mu']))<=D(mu)),key=lambda r:r['mu'])
    upper=min((r for r in ordered if D(str(r['mu']))>=D(mu)),key=lambda r:r['mu'])
    alpha=(D(mu)-D(str(lower['mu'])))/(D(str(upper['mu']))-D(str(lower['mu'])))
    return [str((1-alpha)*D(str(lower[k]))+alpha*D(str(upper[k]))) for k in ['u','v','T']]


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--masses',default='0.613,0.5990,0.5991,0.6178,0.6179,0.597,0.605,0.618')
    parser.add_argument('--digits',type=int,default=65)
    parser.add_argument('--order',type=int,default=40)
    parser.add_argument('--correction-tolerance',default='1e-42')
    parser.add_argument('--verification-tolerance',default='1e-46')
    parser.add_argument('--output',type=Path,default=ROOT/'results/high_precision_orbits.json')
    parser.add_argument('--raw-dir',type=Path,default=ROOT/'results/mp_raw')
    args=parser.parse_args();getcontext().prec=args.digits
    args.raw_dir.mkdir(parents=True,exist_ok=True);args.output.parent.mkdir(parents=True,exist_ok=True)
    result={'status':'provisional numerical verification; not interval-certified or publication-ready',
        'equations':'q_i_ddot=sum_{j!=i} m_j*(q_j-q_i)/|q_j-q_i|^3; G=1; masses=(1,mu,1)',
        'state_order':['q1x','q1y','q2x','q2y','q3x','q3y','v1x','v1y','v2x','v2y','v3x','v3y'],
        'precision_digits':args.digits,'linear_algebra':'Python Decimal, same decimal working precision',
        'integrator':'MPFR adaptive Taylor; local last-two-term estimate; no rigorous global error bound',
        'reduction':'Fix COM and momentum, restrict dH=dL=0, quotient time and rotation using Euclidean gauges.',
        'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'high_precision.py',ROOT/'tools/taylor_mpfr.cpp',ROOT/'tools/mpfr_abi.hpp']},
        'orbits':[]}
    for mu in args.masses.split(','):
        prefix=mu.replace('.','p')+f'_d{args.digits}_o{args.order}'
        p,history=correct(mu,seed_for(mu),args.digits,args.order,args.correction_tolerance,args.raw_dir,prefix)
        dat=verify_orbit(mu,p,history,args.digits,args.order+4,args.verification_tolerance,args.raw_dir,prefix)
        result['orbits'].append(dat);args.output.write_text(json.dumps(result,indent=2))


if __name__=='__main__':main()
