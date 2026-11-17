"""Profile linear nuisance amplitudes while optimizing physical spectral shapes."""
import copy
import numpy as np
from scipy.optimize import minimize

PARAMS={'powerlaw':['index'], 'cutoff_powerlaw':['index','cutoff_MeV'],
        'logparabola':['alpha','beta']}
DEFAULTS={'index':(2.4,.5,5.),'alpha':(2.4,.5,5.),'beta':(.1,0.,1.),
          'cutoff_MeV':(30000.,1000.,1e7)}

def parse_models(raw):
    if not isinstance(raw,dict) or not raw: raise ValueError('spectral_sources must be a nonempty mapping')
    result={}
    for name,s in raw.items():
        if not isinstance(s,dict): raise ValueError('Source specification must be a mapping')
        if set(s)-{'model','seed_index','seed_flux','pivot_MeV','parameters'}: raise ValueError('Unknown source keys')
        model=s.get('model','powerlaw')
        if model not in PARAMS: raise ValueError('Use powerlaw, cutoff_powerlaw, or logparabola')
        r={'model':model,'seed_index':float(s.get('seed_index',2.4)),
           'seed_flux':float(s.get('seed_flux',1e-9)),'pivot_MeV':float(s.get('pivot_MeV',1000.)),
           'parameters':{}}
        if not np.isfinite(r['seed_index']) or not np.isfinite(r['seed_flux']) or r['seed_flux']<=0 or not np.isfinite(r['pivot_MeV']) or r['pivot_MeV']<=0: raise ValueError('Invalid seed/pivot')
        pars=s.get('parameters',{})
        if set(pars)-set(PARAMS[model]): raise ValueError('Parameters do not match model')
        for key in PARAMS[model]:
            v,lo,hi=DEFAULTS[key]; item=pars.get(key,{})
            if not isinstance(item,dict) or set(item)-{'value','free','min','max'}: raise ValueError('Use value/free/min/max parameter mappings')
            q={'value':float(item.get('value',v)),'free':item.get('free',True),
               'min':float(item.get('min',lo)),'max':float(item.get('max',hi))}
            if not isinstance(q['free'],bool) or not np.all(np.isfinite([q['value'],q['min'],q['max']])) or not q['min']<=q['value']<=q['max'] or q['min']>=q['max']: raise ValueError('Invalid parameter bounds/value/free')
            if key=='cutoff_MeV' and q['min']<=0: raise ValueError('Cutoff must be positive')
            r['parameters'][key]=q
        result[str(name)]=r
    return result

def ratio(energy,s):
    """Desired spectrum / baked-in power-law seed, with unit amplitude at pivot."""
    e=np.asarray(energy,float); x=np.log(e/s['pivot_MeV']); p=s['parameters']
    if s['model']=='logparabola': log=(s['seed_index']-p['alpha']['value'])*x-p['beta']['value']*x*x
    else:
        log=(s['seed_index']-p['index']['value'])*x
        if s['model']=='cutoff_powerlaw': log-=(e-s['pivot_MeV'])/p['cutoff_MeV']['value']
    return np.exp(np.maximum(log,-745.))

def native_spectrum(s):
    """Fermitools analytic multiplicative correction to the seed intensity cube."""
    def par(v,lo=-1e20,hi=1e20): return dict(value=float(v),scale=1.,min=float(min(lo,v)),max=float(max(hi,v)),free=False)
    p=s['parameters']; pivot=s['pivot_MeV']
    if s['model']=='logparabola':
        return 'LogParabola',{'norm':par(1.,0.),'alpha':par(p['alpha']['value']-s['seed_index']),
                             'beta':par(p['beta']['value']),'Eb':par(pivot,1e-6)}
    if s['model']=='cutoff_powerlaw':
        cut=p['cutoff_MeV']['value']
        if pivot/cut>100: raise ValueError('Cutoff too far below pivot for stable normalization')
        return 'PLSuperExpCutoff',{'Prefactor':par(np.exp(pivot/cut),0.),
            'Index1':par(s['seed_index']-p['index']['value']),'Index2':par(1.,.1,10.),
            'Scale':par(pivot,1e-6),'Cutoff':par(cut,1e-6)}
    return 'PowerLaw',{'Prefactor':par(1.,0.),'Index':par(s['seed_index']-p['index']['value']),
                       'Scale':par(pivot,1e-6)}

class ProfileSpectra:
    def __init__(self,counts,fixed,templates,names,free,initial,models,predict,joint_class,tied_sources=None):
        self.counts=counts; self.fixed=fixed; self.base=templates
        self.names=list(names); self.free=np.asarray(free,bool); self.initial=np.asarray(initial,float)
        self.models=copy.deepcopy(models); self.predict=predict; self.joint_class=joint_class
        tied_sources=list(tied_sources or [])
        self.ties={name:tied_sources[0] for name in tied_sources[1:]}
        for name,leader in self.ties.items():
            if name not in self.models or leader not in self.models: raise ValueError('Unknown tied spectral source')
            a={k:v for k,v in self.models[name].items() if k!='seed_flux'}
            b={k:v for k,v in self.models[leader].items() if k!='seed_flux'}
            if a!=b: raise ValueError('Tied sources must have identical spectral shapes, bounds and free flags')
        self.keys=[]; start=[]; bounds=[]
        for name,s in self.models.items():
            if name not in self.names: raise ValueError('Spectral source absent from prepared cache: '+name)
            j=self.names.index(name)
            if not tied_sources and not self.free[j] and self.initial[j]==0 and any(p['free'] for p in s['parameters'].values()):
                raise ValueError('Cannot fit spectral shape of a source fixed to zero')
            for key,q in s['parameters'].items():
                if q['free'] and name not in self.ties:
                    self.keys.append((name,key)); trans=np.log if key=='cutoff_MeV' else float
                    start.append(trans(q['value'])); bounds.append((trans(q['min']),trans(q['max'])))
        self.start=np.array(start); self.bounds=bounds; self.calls=0
        self.best=None; self.last_x=None; self.last_result=None

    def evaluate(self,x):
        x=np.asarray(x,float)
        if self.last_x is not None and np.array_equal(x,self.last_x): return self.last_result
        for (name,key),v in zip(self.keys,x):
            self.models[name]['parameters'][key]['value']=float(np.exp(v) if key=='cutoff_MeV' else v)
        for name,leader in self.ties.items():
            self.models[name]['parameters']=copy.deepcopy(self.models[leader]['parameters'])
        t=self.base.copy()
        for name,s in self.models.items(): t[self.names.index(name)+1]=self.predict(name,s)
        fixed=self.fixed+np.einsum('j,jkp->kp',self.initial[~self.free],t[1:][~self.free])
        selected=np.concatenate([t[:1],t[1:][self.free]])
        # Balance linear-template scales as shapes change, especially near a
        # low cutoff. This is only a numerical reparameterization.
        factors=selected[1:].sum(axis=(1,2))/max(self.counts.sum()/max(len(selected),1),1.)
        if np.any(factors<=0): raise ValueError('Zero-support free component')
        balanced=selected.copy(); balanced[1:]/=factors[:,None,None]
        problem=self.joint_class(self.counts,fixed,balanced)
        amplitudes,value=problem.solve()
        amplitudes[self.counts.shape[0]:]/=factors
        self.calls+=1
        result=(float(value),amplitudes,t,problem)
        if self.best is None or value<self.best[0]: self.best=(float(value),x.copy())
        self.last_x=x.copy(); self.last_result=result
        return result

    def solve(self):
        if not len(self.start): return self.start,self.evaluate(self.start)
        def f(x): return self.evaluate(x)[0]/max(float(self.counts.sum()),1.)
        def grad(x):
            g=[]
            for i,(lo,hi) in enumerate(self.bounds):
                xp=x.copy(); xm=x.copy(); xp[i]=min(hi,x[i]+1e-4); xm[i]=max(lo,x[i]-1e-4)
                g.append((f(xp)-f(xm))/(xp[i]-xm[i]))
            return np.array(g)
        starts=[self.start,np.array([(a+b)/2 for a,b in self.bounds])]
        accepted=[]
        for start in starts:
            r=minimize(f,start,jac=grad,method='L-BFGS-B',bounds=self.bounds,
                       options={'ftol':1e-14,'gtol':1e-8,'maxiter':500,'maxls':40})
            if r.success and np.isfinite(r.fun): accepted.append(r)
        if not accepted: raise RuntimeError('Spectral optimization did not converge')
        r=min(accepted,key=lambda q:q.fun); x=r.x
        result=self.evaluate(x)
        if self.best[0]<result[0]-1e-3: raise RuntimeError('A better point was found during optimization; refine starts/bounds')
        # Independent finite-difference projected-gradient check in the optimizer coordinates.
        g=[]
        for i,(lo,hi) in enumerate(self.bounds):
            xp=x.copy(); xm=x.copy(); xp[i]=min(hi,x[i]+1e-4); xm[i]=max(lo,x[i]-1e-4)
            d=(f(xp)-f(xm))/(xp[i]-xm[i])
            if (x[i]<=lo+1e-4 and d>0) or (x[i]>=hi-1e-4 and d<0): d=0.
            g.append(d)
        if max(abs(np.asarray(g)))>1e-6: raise RuntimeError('Spectral projected gradient check failed')
        return x,self.evaluate(x)
