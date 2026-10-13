"""Reproducible exploratory continuation of Golden Pancake Drama orbit.
Units: G=1, q_i=(-1,0),(0,0),(1,0); v1=v3=(u,v); v2=-2*(u,v)/mu.
SciPy double precision; numbers at ~1e-9 precision require independent multiprecision verification.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
try:
    from numba import njit
    FORCE_BACKEND = "numba"
except ImportError:
    FORCE_BACKEND = "numpy (Numba unavailable)"
    def njit(**kwargs):
        return lambda fn: fn

r0=np.array([[-1.,0.],[0.,0.],[1.,0.]])


def xinit(mu,u,v):
    vel=np.array([[u,v],[-2*u/mu,-2*v/mu],[u,v]])
    return np.r_[r0.ravel(),vel.ravel()]

@njit(cache=True)
def f(t,x,m):
    pos=x[:6].reshape(3,2); acc=np.zeros((3,2))
    for i in range(3):
        for j in range(i+1,3):
            d=pos[j]-pos[i]
            r2=d@d
            inv3=1.0/(r2*np.sqrt(r2))
            acc[i]+=m[j]*d*inv3
            acc[j]-=m[i]*d*inv3
    return np.concatenate((x[6:],acc.ravel()))

@njit(cache=True)
def jac(x,m):
    pos=x[:6].reshape(3,2)
    A=np.zeros((12,12));A[:6,6:]=np.eye(6)
    for i in range(3):
        for j in range(i+1,3):
            d=pos[j]-pos[i];r2=d@d
            B=(np.eye(2)-3*d.reshape(2,1)*d.reshape(1,2)/r2)/(r2*np.sqrt(r2))
            ii=slice(6+2*i,8+2*i);jj=slice(6+2*j,8+2*j)
            pi=slice(2*i,2*i+2);pj=slice(2*j,2*j+2)
            A[ii,pi]-=m[j]*B; A[ii,pj]+=m[j]*B
            A[jj,pi]+=m[i]*B; A[jj,pj]-=m[i]*B
    return A


@njit(cache=True)
def ext_numba(t,y,m):
    xx=y[:12];M=y[12:].reshape((12,12))
    return np.concatenate((f(t,xx,m), (jac(xx,m)@M).ravel()))

def integrate(mu,u,v,T,variational=False,rtol=2e-13,atol=3e-15,dense=False,max_step=np.inf):
    mass=np.array([1.,mu,1.]); x=xinit(mu,u,v)
    if variational:
        y=np.r_[x, np.eye(12).ravel()]
        fun=lambda t,y:ext_numba(t,y,mass)
    else:
        y=x;fun=lambda t,y:f(t,y,mass)
    return solve_ivp(fun,(0,T),y,method='DOP853',rtol=rtol,atol=atol,dense_output=dense,max_step=max_step)

def legacy_correction(mu,guess,rtol=3e-13,atol=3e-15):
    mass=np.array([1.,mu,1.]);cache={};history=[]
    def compute(z):
        key=tuple(z)
        if key not in cache:
            sol=integrate(mu,*z,variational=True,rtol=rtol,atol=atol)
            pair_min=min(float(np.min(np.linalg.norm(sol.y[2*j:2*j+2]-sol.y[2*i:2*i+2],axis=0))) for i,j in [(0,1),(0,2),(1,2)])
            history.append({"param":list(map(float,z)),"success":bool(sol.success),"status":int(sol.status),"nfev":int(sol.nfev),"last_time":float(sol.t[-1]),"sampled_min_distance":pair_min})
            if not sol.success:
                error=RuntimeError(sol.message);error.trial_history=history
                raise error
            x0=xinit(mu,z[0],z[1]);xf=sol.y[:12,-1]
            m=sol.y[12:,-1].reshape(12,12)
            dd=np.zeros((12,2));dd[6:,0]=[1,0,-2/mu,0,1,0];dd[6:,1]=[0,1,0,-2/mu,0,1]
            J=np.column_stack([((m-np.eye(12))@dd), f(0,xf,mass)])
            cache.clear();cache[key]=(xf-x0,J,sol)
        return cache[key]
    res=least_squares(lambda z:compute(z)[0],np.asarray(guess,dtype=float),jac=lambda z:compute(z)[1],xtol=2e-14,ftol=2e-14,gtol=2e-14,max_nfev=12)
    rr,jj,sol=compute(res.x)
    sol.legacy_trial_history=history
    return res.x, np.max(np.abs(rr)), sol

def energy(mu,x):
    m=[1,mu,1];r=x[:6].reshape(3,2);v=x[6:].reshape(3,2)
    kinetic=sum(.5*m[i]*(v[i]@v[i]) for i in range(3))
    potential=-sum(m[i]*m[j]/np.linalg.norm(r[j]-r[i]) for i in range(3) for j in range(i+1,3))
    return kinetic+potential

def stats(mu,p,hi_variational=False):
    import time
    started=time.perf_counter()
    u,v,T=p
    sol=integrate(mu,*p,variational=hi_variational,rtol=1.5e-13,atol=2e-15,dense=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    # Minimum separation at 1600 samples; refine with scalar optimizer in time
    from scipy.optimize import minimize_scalar
    times=np.linspace(0,T,1601)
    xyz=sol.sol(times)[:6].reshape(3,2,-1)
    def rdists(t):
        p=sol.sol(t)[:6].reshape(3,2)
        return [np.linalg.norm(p[j]-p[i]) for i,j in [(0,1),(0,2),(1,2)]]
    minimum=(float('inf'),None,None)
    for k,(i,j) in enumerate([(0,1),(0,2),(1,2)]):
        ds=np.linalg.norm(xyz[j]-xyz[i],axis=0)
        indexes=np.r_[0,np.flatnonzero((ds[1:-1]<ds[:-2]) & (ds[1:-1]<ds[2:]))+1,len(times)-1]
        for index in indexes:
            t0=times[max(0,index-1)];t1=times[min(len(times)-1,index+1)]
            opt=minimize_scalar(lambda t:rdists(t)[k],bounds=(t0,t1),method='bounded',options={'xatol':2e-12}) if t1>t0 else None
            t=float(opt.x if opt else times[index]);val=rdists(t)[k]
            if val<minimum[0]:minimum=(val,t,f'{i+1}-{j+1}')
    dat={'mu':mu,'u':u,'v':v,'T':T,'E':energy(mu,xinit(mu,u,v)),'S': T*abs(energy(mu,xinit(mu,u,v)))**1.5,
         'min_distance':minimum[0],'t_nearest':minimum[1],'pair_nearest':minimum[2],
         'closure':float(np.max(np.abs(sol.y[:12,-1]-xinit(mu,u,v))))}
    dat.update(integration_metadata(sol))
    dat['reporting_settings']={'rtol':1.5e-13,'atol':2e-15,'max_step':None}
    dat['reporting_runtime_seconds']=time.perf_counter()-started
    dat['return_vector']=(sol.y[:12,-1]-xinit(mu,u,v)).tolist()
    if hi_variational:
        M=sol.y[12:,-1].reshape(12,12)
        vals=np.linalg.eigvals(M)
        # sort by |lambda-1|, first 8 trivial, remaining 4 should contain interesting eigenvalues
        idx=np.argsort(abs(vals-1))[::-1]
        nontrivial=vals[idx[:4]]
        dat['mult_real']=[float(z.real) for z in nontrivial]
        dat['mult_imag']=[float(z.imag) for z in nontrivial]
        dat['multipliers_all']=[[float(z.real),float(z.imag)] for z in vals]
    return dat



def integration_metadata(sol):
    return {"integration_success": bool(sol.success), "integration_status": int(sol.status),
            "integration_message": sol.message, "nfev": int(sol.nfev),
            "accepted_steps": len(sol.t)-1, "precision": "IEEE binary64", "method": "DOP853",
            "force_backend": FORCE_BACKEND}


def shooting_system(mu, z, rtol=3e-13, atol=3e-15, max_step=np.inf):
    sol = integrate(mu, *z, variational=True, rtol=rtol, atol=atol, max_step=max_step)
    if not sol.success:
        raise RuntimeError(sol.message)
    phi = sol.y[12:, -1].reshape(12, 12)
    dd = np.zeros((12, 2))
    dd[6:, 0] = [1, 0, -2/mu, 0, 1, 0]
    dd[6:, 1] = [0, 1, 0, -2/mu, 0, 1]
    residual = sol.y[:12, -1] - xinit(mu, *z[:2])
    derivative = np.column_stack(((phi-np.eye(12)) @ dd, f(0, sol.y[:12, -1], np.array([1., mu, 1.]))))
    return residual, derivative, sol


def correction(mu, guess, rtol=3e-13, atol=3e-15, max_step=np.inf):
    """Full-state least-squares Newton with trust cap and rejecting failed trial integrations.

    Unconstrained trial steps can leave the isola and enter collision trajectories.
    Stopping at the measured binary64 floor avoids chasing roundoff with a large step.
    """
    import time
    started = time.perf_counter()
    z = np.array(guess, dtype=float)
    history = []
    rr, jj, sol = shooting_system(mu, z, rtol, atol, max_step)
    for iteration in range(24):
        norm = float(np.max(np.abs(rr)))
        history.append({"iteration": iteration, "param": z.tolist(), "residual": norm})
        if norm < 5e-11:
            break
        step = np.linalg.lstsq(jj, -rr, rcond=None)[0]
        step /= max(1., np.max(np.abs(step)/np.array([.0015, .0015, .025])))
        accepted = False
        for power in range(12):
            trial = z + step * 2.**(-power)
            if trial[2] <= 0:
                continue
            try:
                rnew, jnew, snew = shooting_system(mu, trial, rtol, atol, max_step)
            except RuntimeError as exc:
                history.append({"rejected_param": trial.tolist(), "integration_failure": str(exc)})
                continue
            if np.linalg.norm(rnew) < np.linalg.norm(rr):
                z, rr, jj, sol = trial, rnew, jnew, snew
                accepted = True
                break
        if not accepted:
            if norm > 5e-9:
                raise RuntimeError(f"Newton stalled at residual {norm:.3e}; use a closer mass seed")
            break
    residual = float(np.max(np.abs(rr)))
    if residual > 5e-9:
        raise RuntimeError(f"Newton did not converge: residual {residual:.3e}")
    sol.correction_info = {"converged": residual < 5e-9, "target_met": residual < 5e-11,
        "history": history, "runtime_seconds": time.perf_counter()-started,
        "rtol": rtol, "atol": atol, "max_step": None if np.isinf(max_step) else max_step}
    return z, residual, sol


def continue_to(mu_start, mu_target, param, step=.00025, **kwargs):
    """Small mass increments retain the branch instead of making large initial jumps."""
    count = max(1, int(np.ceil(abs(mu_target-mu_start)/step)))
    for mu in np.linspace(mu_start, mu_target, count+1)[1:]:
        param, residual, sol = correction(float(mu), param, **kwargs)
    return param, residual, sol


if __name__ == "__main__":
    import argparse, json, time
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument("--mass-list", default="0.613,0.612,0.611,0.610,0.608,0.605")
    parser.add_argument("--stability", action="store_true")
    parser.add_argument("--out", type=Path, default=Path("results/branch_recomputed_double.json"))
    args = parser.parse_args()
    param = np.array([.151481541079005834, .140303019570684087, 8.846608333332379])
    rows = []; previous = .613; started = time.perf_counter()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    for mu in map(float, args.mass_list.split(",")):
        param, residual, sol = continue_to(previous, mu, param)
        dat = stats(mu, param, hi_variational=args.stability)
        dat.update(correction_residual=residual, correction=sol.correction_info)
        rows.append(dat); previous = mu
        args.out.write_text(json.dumps(rows, indent=2))
        print(f"mu={mu:.7f} closure={dat['closure']:.3e} elapsed={time.perf_counter()-started:.1f}s", flush=True)
