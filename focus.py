import glob
import os
import copy
import numpy as np
from astropy.io import fits
from astropy.time import Time
from pyvista import centroid, stars, tv, imred, image
import matplotlib
import matplotlib.pyplot as plt
import pdb
from photutils.detection import find_peaks, DAOStarFinder
from holtztools import plots, html
import time
import database
import logging
import yaml
import logging.config
import aposong
import eshel
try :
    with open('logging.yml', 'rt') as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)
except :
    print("can't open logging.yml")
logger=logging.getLogger(__name__)


def focus(files, apers=np.arange(0.3,6,0.2), thresh=100, fwhm=2, skyrad=[8,12], 
          pixscale=0.5, display=None, plot=False,max=None,red=None, root='',hard=False) :
    """ Get best focus position from a set of focus run images

        For each image, find stars and do concentric aperture photometry to 
        determine half-flux radius. Take median of these across all stars in
        image and consider result as a function of focus. Find the minimum
        half-flux image and do a quadratic fit around this point to find minimum
        best fit focus. Report best fit and minimum focus and half flux diameters
        in arcsec.
    """

    pixapers=apers/pixscale
    hfmed=np.zeros([len(files)])
    foc=np.zeros([len(files)])
    if plot : fig,ax=plots.multi(1,len(files),hspace=0.001,figsize=(4,6))
    if display is not None :
        display.plotax2.cla()
    for ifile,file in enumerate(files) :
        # read file and get focus value
        if red is None :
            a=fits.open(root+file)[0]
        else :
            a=red.rd(root+file)
        im=a.data.astype(float)
        foc[ifile]=a.header['FOCUS']
        logger.info('focus: {:s} {:f}'.format(file,foc[ifile]))

        # set saturated pixels to NaN
        bd=np.where(im>60000)
        im[bd] = np.nan

        # display
        if display is not None :
            display.tvclear()
            display.tv(a,max=max)

        # find stars
        mad=np.nanmedian(np.abs(im-np.nanmedian(im)))
        daofind=DAOStarFinder(fwhm=fwhm/pixscale,threshold=thresh*mad,exclude_border=True)
        tab=daofind(im[50:-50,50:-50]-np.nanmedian(im))
        if tab is None :
            logger.error('no sources found...')
            continue
        tab['xcentroid']+=50
        tab['ycentroid']+=50
        tab.rename_column('xcentroid','x')
        tab.rename_column('ycentroid','y')
        
        for i,star in enumerate(tab) :
            center=centroid.rasym_centroid(im,star['x'],star['y'],int(pixapers[-1]),skyrad=skyrad)
            tab['x'][i]=center.x
            tab['y'][i]=center.y

        # aperture photometry through range of apertures
        phot=stars.photom(im,tab,skyrad=np.array(skyrad)/pixscale,rad=pixapers,mag=False)
        gd = np.where(np.isfinite(phot['aper{:.1f}'.format(pixapers[-1])]))[0]
        if len(gd) < 1 :
            logger.error('no sources with finite photometry')
            continue
        phot=phot[gd]
        if display is not None :
            display.tvclear()
            stars.mark(display,phot,exit=True)

        #only take stars within factor of 3 in flux of brightest star
        tot = phot['aper{:.1f}'.format(pixapers[-1])]
        gd = np.where(tot > tot.max()/3.)[0]
        phot=phot[gd]
        if display is not None : stars.mark(display,phot,exit=True,color='g')

        # get cumulative photometry curve
        cum=np.zeros([len(phot),len(pixapers)-1])
        for i,aper in enumerate(pixapers[0:-1]) :
            cum[:,i] = phot['aper{:.1f}'.format(aper)]/phot['aper{:.1f}'.format(pixapers[-1])]

        # identify half flux point
        hf=np.zeros(len(phot))
        tot=np.zeros(len(phot))
        for istar in range(len(phot)) :
            if cum[istar].max() > 1.1 : continue
            if plot : ax[ifile].plot(pixapers[0:-1],cum[istar,:].T)
            j=np.where(cum[istar]>0.5)[0]
            if len(j) > 0 :
                j0=j[0]-1
                if j0<0 : continue
                j1=j[0]
                hf[istar] = pixapers[j0] + (0.5-cum[istar,j0])/(cum[istar,j1]-cum[istar,j0])*(pixapers[j1]-pixapers[j0])
                tot[istar] = cum[istar,-1]
                #hfold=(pixapers[j[0]-1]+pixapers[j[0]])/2.

        hfmed[ifile]=np.median(hf)
        if plot :
            ax[ifile].scatter(hf,[0.5]*len(hf),c='k')
            ax[ifile].scatter(np.median(hf),0.5,c='b',s=75)
            ax[ifile].set_ylim(0,1.2)
            ax[ifile].text(0.05,0.8,str(a.header['FOCUS']),transform=ax[ifile].transAxes)
            plt.figure(fig.number)
            plt.draw()
        if display is not None and len(hf) > 0 :
            plots.plotp(display.plotax2,np.array([foc[ifile]]*len(hf)),hf*pixscale*2,
                        xt='Focus',yt='R(half total)',size=5)

    gd=np.where(hfmed>0)[0]
    if len(hfmed) > 3 :
        # get smallest median half-flux radius
        besthf=np.min(hfmed[gd])
        # convert to diameter in arcsec
        besthf *= 2*pixscale
        bestind=np.argmin(hfmed[gd])
        bestfoc=foc[gd[bestind]]
        logger.info('bestfoc: {:f}'.format( bestfoc))
        try: 
            # attempt a quadratic fit around the minimum
            poly=np.polyfit(foc[gd[bestind-2:bestind+3]],hfmed[gd[bestind-2:bestind+3]],2)
            bestfitfoc = -poly[1]/2/poly[0]
            if bestfitfoc < foc[gd[bestind-2]] or bestfitfoc > foc[gd[bestind+2]] :
                logger.warning('best focus out of fit range!')
                bestfitfoc=-1
                bestfithf=-1
            else :
                bestfithf = poly[0]*bestfitfoc**2+poly[1]*bestfitfoc+poly[2]
                bestfithf *= 2*pixscale
            logger.info('bestfoc, bestfitfoc: {:f} {:d} '.format( bestfitfoc,bestfithf))
        except: 
            logger.warning('focus fit failed...')
            bestfitfoc=-1
            bestfithf=-1

        # plot and/or display
        if plot :
            fig2,ax2 = plots.multi(1,1)
            plots.plotp(ax2,foc[gd],hfmed[gd]*2*pixscale,xt='Focus',yt='R(half total)',size=50)
            ax2.set_title('{:s} {:s}'.format(files[0],a.header['DATE-OBS']))
        if display is not None :
            plots.plotp(display.plotax2,foc[gd],hfmed[gd]*pixscale*2,
                        xt='Focus',yt='R(half total)',size=50)
            if bestind>1 and bestind < len(gd)-3 :
                x=np.linspace(foc[gd[bestind-2]],foc[gd[bestind+2]],100)
                plots.plotl(display.plotax2,x,(poly[0]*x**2+poly[1]*x+poly[2])*pixscale*2)
            plt.draw()
    else :
        raise Exception('focus run failed')

    if plot :
        if hard :
            fig2.savefig(root+files[0].replace('.fits','.png'))
        else :
            print('Hit RETURN key to close plot windows and continue: ')
            inp=input()
        plt.close(fig)
        plt.close(fig2)

    return bestfitfoc, bestfithf,  bestfoc, besthf
    # return best fit values and minimum values

def mkplots(mjd,display=None,root='/data/1m/') :

    matplotlib.use('Agg')
    d=database.DBSession()
    out=d.query('obs.focus',fmt='list')
    i =np.where(np.array(out[0]) == 'mjd')[0][0]
    ifile =np.where(np.array(out[0]) == 'files')[0][0]
    mjds=[]
    files=[]
    for o in out[1:] : 
        mjds.append(o[i])
        files.append(o[ifile])

    gd=np.where(np.array(mjds).astype(int) == mjd)[0]

    grid=[]
    row=[]
    for i,seq in enumerate(gd) :
        print(files[seq])
        try : 
            #focus(files[seq],pixscale=0.22,display=display,root=root,plot=True,hard=True)
            if i>0 and i%5 == 0 :
                grid.append(row)
                row=[]
            row.append(os.path.basename(files[seq][0].replace('.fits','.png')))
        except : continue
    dir=files[seq][0].split('/')[0]
    html.htmltab(grid,file=root+dir+'/focus.html',size=250)

def calfocus(foc0=34700,display=None) :
    """ Focus run for calibration channel
    """
    foc=aposong.foc()
    aposong.calstage_in()
    eshel.lamps(quartz=True,led=True)
    time.sleep(3)
    if foc0 == None : foc0=aposong.foc()
    plt.figure()
    figno=plt.gcf().number
    for df in range(-1000,1001,200) :
        print(df)
        aposong.foc(foc0+df)
        out=aposong.gexp(0.001,max=30000,display=display)
        plt.figure(figno)
        plt.plot(out.hdu.data[510],label='{:d}'.format(foc0+df))
    plt.legend()
    aposong.foc(foc0)
    eshel.lamps()
    aposong.calstage_out()
    aposong.foc(foc)


def specfocus(foc0=425000) :
    red=imred.Reducer('SONG',dir='/data/1m/UT251002')
    trace=spectra.Trace('./UT2509xx_Trace.fits')
    wav=spectra.WaveCal('./UT2509xx_WaveCal.fits')

    aposong.specfoc(foc0)

    aposong.calstage_in()
    eshel.lamps(quartz=True,led=True)
    s=aposong.sexp(50,name='flat',display=disp,max=65535)
    t=red.reduce(s.name)
    trace.retrace(t)

    eshel.lamps(thar=True)
    fiber='specfoc'
    for i,df in enumerate(range(-2500,2501,500)) :
        aposong.specfoc(foc0+df)
        s=aposong.sexp(50,name='{:s}_{:d}'.format(fiber,foc0+df),display=disp,max=65535)
        t=red.reduce(s.name)
        tec=trace.extract(t)
        wav=spectra.WaveCal('./UT2509xx_WaveCal.fits')
        wav.identify(tec)
        fig,ax=plots.multi(1,2)
        ax[0].scatter(wav.waves,wav.fwhm,c=wav.pix)
        ax[0].set_ylim(0,10)
        ax[1].scatter(wav.waves,wav.waves/(wav.fwhm*(wav.wave(pixels=[wav.pix,wav.y])-wav.wave(pixels=[wav.pix+1,wav.y]))),c=wav.pix)
        ax[1].set_ylim(0,120000)
        fig.savefig('{:s}_{:d}.png'.format(fiber,foc0+df))
    eshel.lamps()
    aposong.calstage_out()


def montage(display,red=None) :

    if red == None :
        red=imred.Reducer('SONG',dir='/data/1m/UT251016')
        files=glob.glob('/data/1m/UT251016/focus*.fits')
    else :
        files=glob.glob(red.dir+'/focus*.fits')
    dates=[]
    for file in files :
        a=red.rd(file)
        dates.append(Time(a.header['DATE-OBS']).mjd)
    dates
    j=np.argsort(dates)

    focold=1000000
    run=[]
    runs=[]
    for file in np.array(files)[j] :
        foc=int(file.split('.')[0].split('_')[1])
        print(foc)
        if foc < focold :
            runs.append(run)
            run=[]
        run.append(file)
        focold=copy.copy(foc)
    box=image.BOX(n=512,cc=580,cr=620)
    montage=[]
    mfiles=[]
    for run in runs[1:] :
        f,s,m=image.seq(run,red=red,box=box)
        montage.append(m)
        mfiles.append(f)

    for f,m in zip(mfiles,montage) :
        n=f[0].split('.')[1]+'-'+f[-1].split('.')[1]
        tit=f[0].split('.')[0].split('_')[1]+'-'+f[-1].split('.')[0].split('_')[1]+' '+n
        display.clear()
        display.tv(m,max=3000)
        display.tvtext(0,0,tit,ha='left')
        plt.draw()
        #display.imexam()
    return mfiles,montage

def profile(ims,red) :
    fig,ax=plots.multi(1,2,hspace=0.001)
    for i,im in enumerate(ims):
        a=red.rd(im)
        tab=stars.find(a,brightest=1,thresh=500)
        cent= centroid.rasym_centroid(a,tab[0]['x'],tab[0]['y'],rad=25,skyrad=[100,150])
        tab['x'] = cent.x
        tab['y'] = cent.y
        ap=stars.photom(a,tab,rad=np.arange(1,80),skyrad=[100,150],mag=False,cum=True)
        ax[0].plot(ap[1],label='{:d}'.format(im))
        if i == 0: 
            cum0 = ap[1]
        else :
            ax[1].plot(ap[1]/cum0,label='{:d}/{:d}'.format(im,ims[0]))

    for i in range(2) :
        ax[i].grid()
        ax[i].legend()
    fig.suptitle(red.dir)

    return ap


