"""Diagnostics for full-sky Galactic CAR differential-intensity templates."""
from pathlib import Path
import argparse
import json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u


def car_solid_angles(header, shape):
    """Exact rectangular pixel areas for the unrotated Galactic CAR output grid."""
    ny,nx=shape
    w=WCS(header).celestial
    if list(w.wcs.ctype)!=['GLON-CAR','GLAT-CAR'] or not np.isclose(w.wcs.crval[1],0):
        raise ValueError('Diagnostics require equatorial-reference Galactic CAR maps.')
    matrix=w.pixel_scale_matrix
    if not np.allclose(matrix-np.diag(np.diag(matrix)),0):
        raise ValueError('Diagnostics require unrotated CAR pixels.')
    dx,dy=abs(matrix[0,0]),abs(matrix[1,1])
    lat=w.pixel_to_world(np.full(ny,(nx-1)/2),np.arange(ny)).galactic.b.deg
    row=np.deg2rad(dx)*(np.sin(np.deg2rad(np.clip(lat+dy/2,-90,90)))-
                       np.sin(np.deg2rad(np.clip(lat-dy/2,-90,90))))
    omega=np.broadcast_to(row[:,None],shape)
    if not np.isclose(omega.sum(),4*np.pi,rtol=1e-5):
        raise ValueError('Input grid does not cover exactly the full sky.')
    return omega


def plot_templates(sun, moon, outdir, energies=(300.,1000.,10000.)):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots(figsize=(7,5),layout='constrained')
    maps,mapaxes=plt.subplots(2,len(energies),figsize=(5*len(energies),6),
                             subplot_kw={'projection':'mollweide'},layout='constrained',squeeze=False)
    report={}
    for row,(name,path) in enumerate([('Sun',sun),('Moon',moon)]):
        with fits.open(path,memmap=True) as hdus:
            a=hdus[0].data;e=np.asarray(hdus['ENERGIES'].data.field(0),float).ravel()
            omega=car_solid_angles(hdus[0].header,a.shape[1:])
            flux=np.array([np.sum(np.asarray(plane,float)*omega) for plane in a])
            np.savetxt(outdir/f'{name.lower()}_template_spectrum.csv',
                       np.column_stack([e,flux,flux/(4*np.pi),e*e*flux]),delimiter=',',
                       header='energy_MeV,full_sky_effective_flux_ph_cm-2_s-1_MeV-1,mean_intensity_ph_cm-2_s-1_sr-1_MeV-1,E2_flux_MeV_cm-2_s-1',comments='')
            good=flux>0
            ax.loglog(e[good],e[good]**2*flux[good],label=name)
            w=WCS(hdus[0].header).celestial
            xx=np.linspace(-np.pi,np.pi,721);yy=np.linspace(-np.pi/2,np.pi/2,361)
            x,y=np.meshgrid((xx[1:]+xx[:-1])/2,(yy[1:]+yy[:-1])/2)
            px,py=w.world_to_pixel(SkyCoord(l=-np.degrees(x)*u.deg,b=np.degrees(y)*u.deg,frame='galactic'))
            ix=np.rint(px).astype(int)%a.shape[2];iy=np.clip(np.rint(py).astype(int),0,a.shape[1]-1)
            actual=[]
            for col,target in enumerate(energies):
                k=int(np.argmin(abs(np.log(e/target))));actual.append(float(e[k]))
                values=np.asarray(a[k],float)[iy,ix]
                positive=values[values>0];axis=mapaxes[row,col]
                if len(positive):
                    low=max(float(np.percentile(positive,1)),float(positive.max())*1e-6)
                    high=max(float(positive.max()),low*(1+1e-6))
                    im=axis.pcolormesh(xx,yy,np.ma.masked_less_equal(values,0),
                                      norm=LogNorm(low,high),cmap='magma',shading='flat',rasterized=True)
                    maps.colorbar(im,ax=axis,orientation='horizontal',pad=.09,
                                  label=r'ph cm$^{-2}$ s$^{-1}$ sr$^{-1}$ MeV$^{-1}$')
                axis.set_title(f'{name}: {e[k]:.4g} MeV')
                axis.set_xticks(np.radians([-120,-60,0,60,120]))
                axis.set_xticklabels(['120°','60°','0°','300°','240°'])
                axis.grid(alpha=.3);axis.set_xlabel('Galactic longitude (increases left)')
            report[name]={'path':str(Path(path).resolve()),'solid_angle_sr':float(omega.sum()),
                          'map_energies_MeV':actual,'energy_range_MeV':[float(e[0]),float(e[-1])]}
    ax.set(xlabel='Energy [MeV]',ylabel=r'$E^2\int I_{\rm eff}(E,\Omega)\,d\Omega$ [MeV cm$^{-2}$ s$^{-1}$]',
           title='Exposure-weighted all-sky template spectra\nReference normalizations = 1')
    ax.set_ylim(1e-5, 1e-3)
    ax.grid(alpha=.25,which='both');ax.legend()
    for ext in ['png','pdf']:
        fig.savefig(outdir/f'sun_moon_spectra.{ext}',dpi=160)
        maps.savefig(outdir/f'sun_moon_skymaps.{ext}',dpi=160)
    plt.close(fig);plt.close(maps)
    report['convention']='Effective stationary-sky template intensity integrated with exact CAR pixel solid angles; not instantaneous source-frame flux or a fitted spectrum.'
    (outdir/'solar_diagnostics.json').write_text(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directory',type=Path)
    p.add_argument('--energies',nargs='+',type=float,default=[300,1000,10000])
    args=p.parse_args()
    if any(e<=0 for e in args.energies):p.error('Energies must be positive MeV.')
    plot_templates(args.directory/'sun_template.fits',args.directory/'moon_template.fits',
                   args.directory,args.energies)
