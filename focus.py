import glob
import numpy as np
from astropy.io import fits
from pyvista import stars, tv
import matplotlib.pyplot as plt
import pdb
from photutils.detection import find_peaks
from photutils import DAOStarFinder
from holtztools import plots

def focus(files, apers=np.arange(0.3,10,0.2), thresh=25, fwhm=2, skyrad=[10,15], pixscale=0.5, display=None) :

    pixapers=apers/pixscale
    pdb.set_trace()
    hfmed=np.zeros([len(files)])
    foc=np.zeros([len(files)])
    fig,ax=plots.multi(1,len(files),hspace=0.001,figsize=(4,6))
    for ifile,file in enumerate(files) :
        # read file and get focus value
        a=fits.open(file)[0]
        im=a.data.astype(float)
        foc[ifile]=a.header['FOCUS']
        print(file,foc[ifile])

        # set saturated pixels to NaN
        bd=np.where(im>60000)
        im[bd] = np.nan

        # display
        if display is not None :
            display.tvclear()
            display.tv(a)

        # find stars
        mad=np.nanmedian(np.abs(im-np.nanmedian(im)))
        daofind=DAOStarFinder(fwhm=fwhm/pixscale,threshold=thresh*mad,exclude_border=True)
        tab=daofind(im-np.nanmedian(im))
        if tab is None :
            print('no sources found...')
            continue
        tab.rename_column('xcentroid','x')
        tab.rename_column('ycentroid','y')
        #tab=stars.find(a,thresh=500,fwhm=5,sharp=[0.5,1.2],round=[-1,1])
        #tab=find_peaks(a,3000,box_size=21)
        #tab.rename_column('x_peak','x')
        #tab.rename_column('y_peak','y')

        # aperture photometry through range of apertures
        phot=stars.photom(im,tab,skyrad=np.array(skyrad)/pixscale,rad=pixapers,mag=False)
        gd = np.where(np.isfinite(phot['aper{:.1f}'.format(pixapers[-1])]))[0]
        phot=phot[gd]
        if display is not None :
            display.tvclear()
            stars.mark(display,phot,exit=True)

        # get cumulative photometry curve
        cum=np.zeros([len(phot),len(pixapers)-1])
        for i,aper in enumerate(pixapers[0:-1]) :
            cum[:,i] = phot['aper{:.1f}'.format(aper)]/phot['aper{:.1f}'.format(pixapers[-1])]
            #ax[ifile].scatter([aper]*len(phot),phot['aper{:.1f}'.format(aper)]/phot['aper29.0'])

        # identify half flux point
        hf=np.zeros(len(phot))
        for istar in range(len(phot)) :
            if cum[istar].max() > 1.1 : continue
            ax[ifile].plot(pixapers[0:-1],cum[istar,:].T)
            j=np.where(cum[istar]>0.5)[0]
            if len(j) > 0 :
                j0=j[0]-1
                if j0<0 : continue
                j1=j[0]
                hf[istar] = pixapers[j0] + (0.5-cum[istar,j0])/(cum[istar,j1]-cum[istar,j0])*(pixapers[j1]-pixapers[j0])
                #hfold=(pixapers[j[0]-1]+pixapers[j[0]])/2.
        gd = np.where(hf>0)[0]
        ax[ifile].scatter(hf[gd],[0.5]*len(gd),c='k')
        ax[ifile].scatter(np.median(hf[gd]),0.5,c='b',s=75)
        ax[ifile].set_ylim(0,1.2)
        ax[ifile].text(0.05,0.8,str(a.header['FOCUS']),transform=ax[ifile].transAxes)
        hfmed[ifile]=np.median(hf[gd])
        plt.figure(fig.number)
        plt.draw()
        #fig.canvas.flush_events()

    if len(hfmed) > 0 :
        fig2,ax2 = plots.multi(1,1)
        gd=np.where(hfmed>0)[0]
        plots.plotp(ax2,foc[gd],hfmed[gd],xt='Focus',yt='R(half total)',size=50)
        if display is not None :
            display.plotax2.cla()
            plots.plotp(display.plotax2,foc[gd],hfmed[gd],xt='Focus',yt='R(half total)',size=50)

    print('Hit RETURN key to close plot windows and continue: ')
    inp=input()
    plt.close(fig)
    plt.close(fig2)
