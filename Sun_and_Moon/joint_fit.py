"""Joint Poisson fit: global Galactic amplitudes + per-bin effective isotropic.

No instrument-response calculation takes place here. Input Galactic and catalog
count templates must ALREADY include the requested energy-dispersion response.
The effective isotropic template must have energy dispersion disabled.
"""
import numpy as np
from scipy.optimize import minimize, brentq
from scipy.special import xlogy


class JointPoisson:
    def __init__(self, counts, fixed, templates):
        self.n, self.fixed, self.t = [np.asarray(v, float) for v in (counts, fixed, templates)]
        if (self.n.ndim != 2 or self.fixed.shape != self.n.shape or self.t.ndim != 3
                or self.t.shape[1:] != self.n.shape or len(self.t) < 1):
            raise ValueError('Expected counts/fixed [energy,pixel], templates [component,energy,pixel].')
        if any(np.any(~np.isfinite(v)) or np.any(v < 0) for v in (self.n,self.fixed,self.t)):
            raise ValueError('Invalid counts/templates.')
        if np.any((self.fixed+self.t.sum(axis=0)==0) & (self.n>0)):
            raise ValueError('Counts without model support.')
        self.ne = self.n.shape[0]
        self.ng = len(self.t)-1
        self.npar = self.ne+self.ng
        self.scale = max(float(self.n.sum()),1.)
        if np.any(self.t[0].sum(axis=1)<=0) or np.any(self.t[1:].sum(axis=(1,2))<=0):
            raise ValueError('Fitted template without support.')

    def model(self,a):
        return self.fixed + a[:self.ne,None]*self.t[0] + np.einsum('j,jkp->kp',a[self.ne:],self.t[1:])

    def objective(self,a):
        mu = self.model(a)
        if np.any((mu<=0)&(self.n>0)):
            return np.inf, np.full(self.npar,-1e100)
        # Half the Poisson deviance differs from NLL by a constant.
        value = np.sum(mu-self.n+xlogy(self.n,self.n)-xlogy(self.n,mu))
        r = 1.-np.divide(self.n,mu,out=np.zeros_like(mu),where=mu>0)
        g = np.r_[(r*self.t[0]).sum(axis=1), np.einsum('kp,jkp->j',r,self.t[1:])]
        return float(value),g

    def covariance(self,a):
        mu = self.model(a)
        w = np.divide(self.n,mu**2,out=np.zeros_like(mu),where=mu>0)
        h = np.zeros((self.npar,self.npar))
        h[np.arange(self.ne),np.arange(self.ne)] = (w*self.t[0]**2).sum(axis=1)
        cross = np.einsum('kp,kp,jkp->kj',w,self.t[0],self.t[1:])
        h[:self.ne,self.ne:] = cross
        h[self.ne:,:self.ne] = cross.T
        h[self.ne:,self.ne:] = np.einsum('kp,ikp,jkp->ij',w,self.t[1:],self.t[1:])
        if np.linalg.cond(h)>1e12:
            return np.full(h.shape,np.nan)
        return np.linalg.inv(h)

    def solve(self,start=None,fixed_parameter=None):
        start = np.ones(self.npar) if start is None else np.asarray(start,float).copy()
        free = np.ones(self.npar,bool)
        if fixed_parameter is not None:
            i,v = fixed_parameter
            free[i] = False
            start[i] = v
        def expand(x):
            a=start.copy(); a[free]=x; return a
        def fg(x):
            v,g=self.objective(expand(x))
            return v/self.scale,g[free]/self.scale
        if not np.any(free):
            return start,self.objective(start)[0]
        solutions=[]
        for x in [start[free],np.ones(free.sum())]:
            r=minimize(fg,x,jac=True,bounds=[(0,None)]*len(x),method='L-BFGS-B',
                       options={'ftol':1e-15,'gtol':2e-11,'maxiter':5000,'maxls':60,
                                'maxcor':min(30,len(x))})
            g=fg(r.x)[1]
            pg=np.where((r.x<1e-8)&(g>0),0.,g)
            if r.success and np.isfinite(r.fun) and np.max(np.abs(pg))<2e-7:
                solutions.append(r)
        if not solutions:
            raise RuntimeError('Joint optimizer failed convergence/KKT checks.')
        r=min(solutions,key=lambda r:r.fun)
        a=expand(r.x)
        return a,self.objective(a)[0]

    def profile(self,index,best,minimum):
        def delta(v):
            # A zero component can make a positive-count pixel impossible.
            a=best.copy(); a[index]=v
            support=self.fixed.copy()
            if index<self.ne:
                support+=self.t[0]
                support[index]=self.fixed[index]+v*self.t[0,index]
                support+=self.t[1:].sum(axis=0)
            else:
                support+=self.t[0]+self.t[1:].sum(axis=0)
                support+=(v-1.)*self.t[index-self.ne+1]
            if np.any((support<=0)&(self.n>0)): return np.inf
            _,f=self.solve(best,(index,v))
            return f-minimum
        a=best[index]; d0=delta(0.)
        if d0 < -1e-3:
            raise RuntimeError('Profile found a better minimum; global fit must be refined.')
        def upper(level):
            hi=max(1.,2*a+.1)
            for _ in range(60):
                if delta(hi)>=level:
                    return brentq(lambda v:delta(v)-level,a,hi,xtol=1e-7)
                hi*=2
            raise RuntimeError('Unbounded profile interval / non-identifiable component.')
        lo=0. if d0<=.5 or a<1e-10 else brentq(lambda v:delta(v)-.5,0.,a,xtol=1e-7)
        return lo,upper(.5),upper(1.352771727),max(0.,2*d0)
