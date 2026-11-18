"""Loop I configuration and optional native spectral fitting for fit_isotropic."""
from pathlib import Path
import json
import numpy as np
from spectral_core import parse_models

ASSETS=Path(__file__).resolve().parents[1]/'LoopI'


def asset_directory(settings):
    """Default assets live beside Isotropic; read_config resolves overrides."""
    return Path(settings.get('directory',ASSETS))


def input_files(c):
    """Inputs underlying the prepared cube and the exported spectral coefficient."""
    s=c.get('loopI',{})
    if not s.get('enabled',False): return []
    assets=asset_directory(s)
    suffix='empirical' if s['spectrum']=='empirical' else 'powerlaw_seed'
    files=[assets/f'loopI_{suffix}_mapcube.fits',assets/'model_metadata.json']
    if s['spectrum']=='empirical': files.append(assets/'empirical_intensity.txt')
    if s.get('mode','combined')=='separate':
        files += [assets/f'loopI_{name}_spatial.fits' for name in ('combined','S1','S2')]
    return files

def source_specs(c):
    """Active source names and independent normalization controls."""
    s=c.get('loopI',{})
    if not s.get('enabled',False): return {}
    if s.get('mode','combined')=='combined':
        return {'LoopI':{'free_norm':s.get('free_norm',False),'normalization':s.get('normalization',1.)}}
    return {f'LoopI_{shell}':{key:s.get('shells',{}).get(shell,{}).get(key,s.get(key,default))
            for key,default in [('free_norm',False),('normalization',1.)]} for shell in ('S1','S2')}


def configure(c):
    if 'loopI' not in c: return
    s=dict(c['loopI'])
    allowed={'enabled','mode','spectrum','free_norm','normalization','parameters','shells','directory'}
    if set(s)-allowed: raise ValueError('Unknown loopI options: '+str(set(s)-allowed))
    s.setdefault('enabled',True); s.setdefault('spectrum','empirical'); s.setdefault('free_norm',False); s.setdefault('normalization',1.)
    # Missing mode preserves the original combined-component configuration/cache.
    if s.get('mode','combined') not in ('combined','separate'): raise ValueError('loopI.mode must be combined or separate')
    if not isinstance(s['enabled'],bool) or not isinstance(s['free_norm'],bool): raise ValueError('Loop I enabled/free_norm must be booleans')
    if not np.isfinite(s['normalization']) or s['normalization']<0: raise ValueError('Invalid Loop I normalization')
    shells=s.get('shells',{})
    if not isinstance(shells,dict) or set(shells)-{'S1','S2'}: raise ValueError('loopI.shells accepts S1 and S2 only')
    if shells and s.get('mode','combined')!='separate': raise ValueError('loopI.shells requires mode: separate')
    for shell,q in shells.items():
        if not isinstance(q,dict) or set(q)-{'free_norm','normalization'}: raise ValueError('Each shell accepts free_norm and normalization only; the spectral shape is shared')
        if 'free_norm' in q and not isinstance(q['free_norm'],bool): raise ValueError('Shell free_norm must be boolean')
        if 'normalization' in q and (not np.isfinite(q['normalization']) or q['normalization']<0): raise ValueError('Invalid shell normalization')
    c['loopI']=s
    if not s['enabled']: return
    reserved={'LoopI','LoopI_S1','LoopI_S2'}
    if reserved & (set(c.setdefault('additional_diffuse',{}))|set(c['galprop'])): raise ValueError('Loop I source already exists; do not add it twice')
    if s['spectrum']=='empirical':
        if s.get('parameters'): raise ValueError('Empirical curve has no adjustable shape; use free_norm')
        suffix='empirical'
    else:
        model=analytic_model(s)
        s['parameters']=model['parameters']
        shape_free=any(q['free'] for q in s['parameters'].values())
        if c.get('compute_profiles',False) and shape_free:
            raise ValueError('Free Loop I shape currently supports best fits only; set compute_profiles: false. Fixed-shape normalization profiling is supported.')
        if shape_free and all(q['normalization']==0 and not q['free_norm'] for q in source_specs(c).values()):
            raise ValueError('Cannot fit spectral shape with all Loop I sources fixed at zero')
        suffix='powerlaw_seed'
    for name,q in source_specs(c).items():
        filename=f'loopI_{suffix}_mapcube.fits' if name=='LoopI' else f'loopI_{name.split("_")[-1]}_{suffix}_mapcube.fits'
        c['additional_diffuse'][name]={'file':str(asset_directory(s)/filename),'free':q['free_norm'],'normalization':q['normalization']}


def spatial_integral(name,assets=ASSETS):
    """Solid-angle integral of the dimensionless map, in sr."""
    meta=json.loads((assets/'model_metadata.json').read_text())
    if name=='LoopI': return meta['allsky_map_integral_sr']
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(assets/f'loopI_{name.split("_")[-1]}_spatial.fits') as h:
        a=np.asarray(h[0].data,float); w=WCS(h[0].header).celestial
        yy,xx=np.indices(a.shape); _,b=w.all_pix2world(xx,yy,0)
        dx,dy=abs(w.wcs.cdelt[0]),abs(w.wcs.cdelt[1])
        om=np.deg2rad(dx)*(np.sin(np.deg2rad(b+dy/2))-np.sin(np.deg2rad(b-dy/2)))
        return float(np.sum(a*om))


def analytic_model(s,name='LoopI'):
    assets=asset_directory(s)
    meta=json.loads((assets/'model_metadata.json').read_text())
    raw={'model':s['spectrum'],'seed_index':meta['seed_index'],
         'seed_flux':meta['intensity_at_1000_MeV']*spatial_integral(name,assets),
         'pivot_MeV':1000.,'parameters':s.get('parameters',{})}
    return parse_models({name:raw})[name]


def ensure_cubes(c):
    """Build selected separate-shell cubes during prepare; no extra CLI stage.

    Preserve the existing common spatial normalization: S1+S2 at unit
    amplitudes equals the combined model, rather than normalizing each to one.
    """
    s=c.get('loopI',{})
    if not s.get('enabled',False): return
    from astropy.io import fits
    import os,tempfile
    suffix='empirical' if s['spectrum']=='empirical' else 'powerlaw_seed'
    assets=asset_directory(s)
    base=assets/f'loopI_{suffix}_mapcube.fits'
    combined=assets/'loopI_combined_spatial.fits'
    missing=[str(p) for p in input_files(c) if not p.is_file()]
    if missing:
        raise FileNotFoundError('Missing Loop I assets: '+', '.join(missing)+
            '. Run LoopI/build_templates.py --config <isotropic-config.yaml> --outdir '+str(assets)+
            ' before prepare (see Isotropic/README.md).')
    if s.get('mode','combined')!='separate': return
    for name in source_specs(c):
        shell=name.split('_')[-1]; spatial=assets/f'loopI_{shell}_spatial.fits'
        target=Path(c['additional_diffuse'][name]['file'])
        if target.exists() and target.stat().st_mtime_ns>=max(p.stat().st_mtime_ns for p in (base,combined,spatial)): continue
        print('Generating separate Loop I input cube:',target,flush=True)
        with fits.open(base,memmap=True) as h, fits.open(spatial) as sh, fits.open(combined) as co:
            frac=np.divide(np.asarray(sh[0].data,float),co[0].data,out=np.zeros_like(sh[0].data,dtype=float),where=co[0].data>0)
            header=h[0].header.copy(); header['SHELL']=shell
            header['HISTORY']='Separate contribution; common S1+S2 spatial normalization retained.'
            data=np.empty(h[0].data.shape,dtype='float32')
            for k,plane in enumerate(h[0].data): data[k]=plane*frac
            hdul=fits.HDUList([fits.PrimaryHDU(data,header),h['ENERGIES'].copy()])
            fd,tmp=tempfile.mkstemp(prefix='loopI_build_',suffix='.fits',dir=target.parent); os.close(fd)
            try:
                hdul.writeto(tmp,overwrite=True); os.replace(tmp,target)
            finally:
                if Path(tmp).exists(): Path(tmp).unlink()


def analytic_enabled(c):
    s=c.get('loopI',{})
    return s.get('enabled',False) and s.get('spectrum','empirical')!='empirical'


def optimize_templates(c,z,keep,driver,output):
    """Reuse prepared source maps; refold analytic shapes BEFORE the linear fit."""
    from native_spectra import NativePredictor,check_seed
    from spectral_core import ProfileSpectra
    from fermipy.gtanalysis import GTAnalysis
    from scipy.special import xlogy
    import yaml,copy
    models={name:analytic_model(c['loopI'],name) for name in source_specs(c)}
    for name,s in models.items(): check_seed(c['additional_diffuse'][name]['file'],s)
    counts=z['counts']; fixed=z['fixed']; templates=z['templates'].copy(); names=list(z['names'])
    free,initial=driver.galprop_controls(c,z['names'][1:])
    cfg=yaml.safe_load((Path(c['outdir'])/'fermipy_config.yaml').read_text())
    gta=GTAnalysis(cfg); gta.setup(init_sources=False); gta.free_sources(free=False)
    if len(gta.components)!=1 or not gta.components[0].like.logLike.use_edisp(): raise ValueError('Expected one component with energy dispersion')
    comp=gta.components[0]
    np.testing.assert_allclose(comp.counts_map().data,counts,rtol=0,atol=0)
    np.testing.assert_allclose(gta.energies,z['edges'])
    for n in names: gta.set_norm_scale(n,1.); gta.set_norm(n,1.,update_source=False)
    total=fixed+templates.sum(axis=0)
    np.testing.assert_allclose(comp.model_counts_map(use_mask=False).data,total,rtol=2e-4,atol=1e-5)
    predictor=NativePredictor(gta,keep,counts.shape)
    gates=[]
    for name,s in models.items():
        neutral=parse_models({name:{'model':'powerlaw','seed_index':s['seed_index'],'seed_flux':s['seed_flux'],
            'pivot_MeV':s['pivot_MeV'],'parameters':{'index':{'value':s['seed_index'],'free':False}}}})[name]
        old=templates[names.index(name)]
        np.testing.assert_allclose(predictor.full(name,neutral),old,rtol=2e-4,atol=1e-5)
        baseline=float(gta.like()); probe=copy.deepcopy(neutral); probe['parameters']['index']['value']+=.2
        for label,trial in [('index_perturbation',probe),('requested_model',s)]:
            changed=predictor.full(name,trial)
            if label=='index_perturbation' and np.allclose(changed,old,rtol=1e-5,atol=1e-10): raise RuntimeError('Native spectrum update did not change prediction')
            mu=total+changed-old
            if np.any((mu<=0)&(counts>0)): raise RuntimeError('Invalid spectral probe model')
            explicit=float(np.sum(mu-total-xlogy(counts,np.divide(mu,total,out=np.ones_like(mu),where=total>0))))
            native=float(gta.like())-baseline
            if abs(native-explicit)>.05+.002*abs(native): raise RuntimeError('Native spectral likelihood validation failed')
            gates.append(dict(name=name,probe=label,native_delta_nll=native,explicit_delta_nll=explicit))
            predictor.full(name,neutral)
    problem=ProfileSpectra(counts[:,keep],fixed[:,keep],templates[:,:,keep],names[1:],free,initial,models,predictor,driver.JointPoisson,tied_sources=list(models))
    shape,(minimum,best,_,_)=problem.solve()
    for name in models: templates[names.index(name)]=predictor.full(name,problem.models[name])
    report={'models':problem.models,'shape_free':bool(problem.keys),'minimum_half_deviance':minimum,
        'native_response_checks':gates,'shared_shape_sources':list(models),'uncertainties':'Linear covariance conditional on best-fit shape; no shape errors or jointly profiled intervals.'}
    (output/'loopI_spectral_fit.json').write_text(json.dumps(report,indent=2))
    np.savez_compressed(output/'spectral_best_fit_parameters.npz',shape=shape,linear=best,minimum=minimum)
    return templates,report


def export_spectrum(c,energies,normalizations,spectral_report,output):
    if not c.get('loopI',{}).get('enabled',False): return
    from spectral_core import ratio
    s=c['loopI']; assets=asset_directory(s)
    norms={r['name']:r['normalization'] for r in normalizations}
    if s['spectrum']=='empirical':
        table=np.loadtxt(assets/'empirical_intensity.txt')
        values=np.exp(np.interp(np.log(energies),np.log(table[:,0]),np.log(table[:,1])))
    else:
        meta=json.loads((assets/'model_metadata.json').read_text())
        model=next(iter(spectral_report['models'].values()))
        values=meta['intensity_at_1000_MeV']*(energies/1000.)**(-model['seed_index'])*ratio(energies,model)
    for name in source_specs(c):
        np.savetxt(output/f'{name.lower().replace("loopi","loopI")}_spectrum.txt',np.c_[energies,norms[name]*values],
            header=f'Energy_MeV spectral_coefficient_ph_cm-2_s-1_sr-1_MeV-1; multiply by corresponding {name} spatial contribution; NOT its own sky-mean intensity')


def export_correlations(c,parameter_names,covariance,output):
    active=source_specs(c)
    if not active: return
    names=list(parameter_names); cov=np.asarray(covariance,float)
    valid=bool(np.all(np.isfinite(cov)) and np.all(np.diag(cov)>0))
    report={'covariance_valid':valid,'conditional_on_shared_shape':bool(analytic_enabled(c) and any(q['free'] for q in c['loopI'].get('parameters',{}).values())),
            'free_LoopI_parameters':[name for name in active if name in names],
            'correlations':{}}
    if valid:
        corr=cov/np.sqrt(np.diag(cov)[:,None]*np.diag(cov)[None,:])
        for source in active:
            if source in names:
                j=names.index(source)
                report['correlations'][source]={name:float(corr[j,i]) for i,name in enumerate(names) if i!=j}
    (output/'loopI_correlations.json').write_text(json.dumps(report,indent=2))
