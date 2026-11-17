"""Rebuild supplied intensity cubes; not needed for normal prepare/run/fit/plot."""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import PchipInterpolator
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import yaml
from shell_geometry import path_length,SHELLS,basic_mask

HERE=Path(__file__).resolve().parent

def curve(energy):
    tab=np.loadtxt(HERE/'leptonic_curve.txt'); e,y=tab[:,:2].T
    energy=np.asarray(energy,float)
    if np.any(energy<e[0]): raise ValueError('Below extracted curve')
    f=PchipInterpolator(np.log(e),np.log(y),extrapolate=False)
    out=np.empty_like(energy); inside=energy<=e[-1]
    out[inside]=np.exp(f(np.log(energy[inside])))
    # Anchor an exponential-cutoff continuation at the last vector point.
    # Fit its tail slope/curvature; this is padding, not published information.
    x=np.log(e[-6:]/e[-1]); z=e[-6:]/e[-1]-1
    a,b=np.linalg.lstsq(np.c_[x,z],np.log(y[-6:]/y[-1]),rcond=None)[0]
    if b>=0: raise ValueError('Tail does not have a decaying exponential cutoff')
    out[~inside]=y[-1]*np.exp(a*np.log(energy[~inside]/e[-1])+b*(energy[~inside]/e[-1]-1))
    return out

def main():
    w=WCS(naxis=2); w.wcs.ctype=['GLON-CAR','GLAT-CAR']; w.wcs.crval=[0,0]
    w.wcs.crpix=[360.5,180.5]; w.wcs.cdelt=[-.5,.5]
    yy,xx=np.indices((360,720)); l,b=w.all_pix2world(xx,yy,0); l=(l+180)%360-180
    om=np.deg2rad(.5)*(np.sin(np.deg2rad(b+.25))-np.sin(np.deg2rad(b-.25)))
    raw=np.array([path_length(l,b,s) for s in SHELLS]); ref=abs(b)>10
    # One shared volume emissivity. The SUM has mean one over this stated domain.
    scale=np.sum(raw.sum(axis=0)[ref]*om[ref])/om[ref].sum()
    shells=raw/scale; combined=shells.sum(axis=0)
    hdr=w.to_header(); hdr['BUNIT']='1'; hdr['HISTORY']='S1+S2 uniform-volume LOS lengths / combined mean at |b|>10 deg.'
    for name,t in [('combined',combined),('S1',shells[0]),('S2',shells[1])]:
        fits.writeto(HERE/f'loopI_{name}_spatial.fits',t.astype('float32'),hdr,overwrite=True)
    energies=np.geomspace(10,1e6,161); pivot=float(curve(np.array([1000.]))[0])
    np.savetxt(HERE/'empirical_intensity.txt',np.c_[energies,curve(energies)],header='MeV ph_cm-2_s-1_sr-1_MeV-1; reference-region mean intensity')
    for name,intensity in [('empirical',curve(energies)),('powerlaw_seed',pivot*(energies/1000)**-2.4)]:
        h=hdr.copy(); h['BUNIT']='ph cm-2 s-1 sr-1 MeV-1'; h['CTYPE3']='ENERGY'; h['CUNIT3']='MeV'; h['CRPIX3']=1.; h['CRVAL3']=energies[0]; h['CDELT3']=1.
        h['HISTORY']='Actual energy samples: ENERGIES extension. See MODEL_NOTES.md.'
        eh=fits.BinTableHDU.from_columns([fits.Column(name='Energy',format='D',unit='MeV',array=energies)],name='ENERGIES')
        fits.HDUList([fits.PrimaryHDU((intensity[:,None,None]*combined).astype('float32'),h),eh]).writeto(HERE/f'loopI_{name}_mapcube.fits',overwrite=True)
    config=yaml.safe_load((HERE.parent/'config.yaml').read_text()); spec=config['mask']
    dec=SkyCoord(l*u.deg,b*u.deg,frame='galactic').icrs.dec.deg
    keep=basic_mask(l,b,dec,spec); roi=(abs(l)<=20)&(abs(b)<=20)
    report={'intensity_at_1000_MeV':pivot,'seed_index':2.4,'allsky_map_integral_sr':float(np.sum(combined*om)),
        'reference_region':'All longitudes, |b|>10 deg; our explicit convention, not a recovered paper footprint',
        'reference_solid_angle_sr':float(om[ref].sum()),'reference_mean':float(np.sum(combined[ref]*om[ref])/om[ref].sum()),
        'curve_max_MeV':float(np.loadtxt(HERE/'leptonic_curve.txt')[-1,0]),
        'mask':spec,'roi_assumption':'Galactic-centre aligned rectangle |l|<=20, |b|<=20; confirm if your nominal ROI has a different centre',
        'intrinsic_flux_fractions':{}}
    for name,t in [('S1',shells[0]),('S2',shells[1]),('combined',combined)]:
        weight=t*om
        report['intrinsic_flux_fractions'][name]={'isotropic_fit':float(weight[keep].sum()/weight.sum()),'GC_40x40':float(weight[roi].sum()/weight.sum()),'both':float(weight[keep&roi].sum()/weight.sum())}
    (HERE/'model_metadata.json').write_text(json.dumps(report,indent=2))
    fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained'); tab=np.loadtxt(HERE/'leptonic_curve.txt')
    for ax in axs:
        ax.loglog(tab[:,0],tab[:,2],color='#9467bd',ls='--',label='Extracted purple leptonic curve')
        ax.loglog(energies,energies**2*curve(energies),color='#173f5f',lw=1,label='MapCube spectrum (PCHIP + tail)')
        ax.axvspan(300,100000,color='grey',alpha=.1,label='Configured fit range')
        ax.axvline(tab[-1,0],color='grey',ls=':',label='Last extracted point')
        ax.set(xlabel='Energy [MeV]',ylabel=r'$E^2 I(E)$ [MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]'); ax.grid(alpha=.2)
    axs[0].set(xlim=(.1,1e6),ylim=(1e-10,.003)); axs[1].set(xlim=(300,3e5),ylim=(1e-12,.003)); axs[0].legend(fontsize=8)
    fig.savefig(HERE/'loopI_spectrum.png',dpi=170); plt.close(fig)
    fig,axs=plt.subplots(3,2,figsize=(13,11),layout='constrained')
    for row,(name,t) in enumerate([('S1',shells[0]),('S2',shells[1]),('S1 + S2',combined)]):
        for col in range(2):
            ax=axs[row,col]; im=ax.imshow(t,extent=[180,-180,-90,90],origin='lower',aspect='auto',cmap='magma')
            ax.contourf(l,b,(~keep).astype(int),levels=[.5,1.5],colors='none',hatches=['///'],alpha=.2)
            ax.add_patch(Rectangle((-20,-20),40,40,fill=False,edgecolor='cyan',lw=2))
            ax.set(xlabel='Galactic longitude [deg]',ylabel='Galactic latitude [deg]',title=name)
            if col: ax.set(xlim=(60,-60),ylim=(-40,70))
            fig.colorbar(im,ax=ax,label='Dimensionless intensity multiplier')
    fig.suptitle('Uploaded isotropic mask: hatching = excluded\nCyan = assumed GC-centred 40 x 40 degree ROI')
    fig.savefig(HERE/'loopI_overlap.png',dpi=150); plt.close(fig)
    print(json.dumps(report,indent=2))
if __name__=='__main__': main()
