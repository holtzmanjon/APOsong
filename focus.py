import glob
import numpy as np
from astropy.io import fits
from pyvista import stars, tv
import matplotlib.pyplot as plt
import pdb
from photutils.detection import find_peaks
from photutils import DAOStarFinder
from holtztools import plots

def focus(files, apers=np.arange(1,15,0.3), thresh=50, fwhm=5, skyrad=[20,30], display=None) :

    hfmed=np.zeros([len(files)])
    foc=np.zeros([len(files)])
    fig,ax=plots.multi(1,len(files),hspace=0.001)
    for ifile,file in enumerate(files) :
        a=fits.open(file)[0]
        foc[ifile]=a.header['FOCUS']
        bd=np.where(a.data>60000)
        a.data[bd] = 0

        if display is not None :
            display.tv(a)
        print('find')

        mad=np.nanmedian(np.abs(a.data-np.nanmedian(a.data)))
        daofind=DAOStarFinder(fwhm=fwhm,threshold=thresh*mad)
        tab=daofind(a.data-np.nanmedian(a.data))
        if tab is None :
            print('no sources found...')
            continue
        gd=np.where((tab['xcentroid']>50)&(tab['ycentroid']>50)&
                    (tab['xcentroid']<a.data.shape[1]-50)&
                    (tab['ycentroid']<a.data.shape[0]-50))[0]
        tab=tab[gd]
        tab.rename_column('xcentroid','x')
        tab.rename_column('ycentroid','y')

        #tab=stars.find(a,thresh=500,fwhm=5,sharp=[0.5,1.2],round=[-1,1])
        #tab=find_peaks(a,3000,box_size=21)
        #tab.rename_column('x_peak','x')
        #tab.rename_column('y_peak','y')
        print('photom')
        phot=stars.photom(a.data,tab,skyrad=skyrad,rad=apers,mag=False)
        gd = np.where(np.isfinite(phot['aper{:.1f}'.format(apers[-1])]))[0]
        phot=phot[gd]
        if display is not None :
            display.tvclear()
            stars.mark(display,phot,exit=True)
        print('cum')
        cum=np.zeros([len(phot),len(apers)-1])
        for i,aper in enumerate(apers[0:-1]) :
            cum[:,i] = phot['aper{:.1f}'.format(aper)]/phot['aper{:.1f}'.format(apers[-1])]
            #ax[ifile].scatter([aper]*len(phot),phot['aper{:.1f}'.format(aper)]/phot['aper29.0'])
        print('hf')
        hf=np.zeros(len(phot))
        for istar in range(len(phot)) :
            if cum[istar].max() > 1.1 : continue
            ax[ifile].plot(apers[0:-1],cum[istar,:].T)
            j=np.where(cum[istar]>0.5)[0]
            if len(j) > 0 :
                hf[istar]=(apers[j[0]-1]+apers[j[0]])/2.
        ax[ifile].scatter(hf,[0.5]*len(phot),c='k')
        ax[ifile].text(0.05,0.9,str(a.header['FOCUS']),transform=ax[ifile].transAxes)
        hfmed[ifile]=np.median(hf)
        plt.draw()
        fig.canvas.flush_events()

    if len(hfmed) > 0 :
        plt.figure()
        plt.plot(foc,hfmed)
