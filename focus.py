import glob
import numpy as np
from astropy.io import fits
from pyvista import centroid, stars, tv
import matplotlib.pyplot as plt
import pdb
from photutils.detection import find_peaks
from photutils import DAOStarFinder
from holtztools import plots

def focus(files, apers=np.arange(0.3,4,0.2), thresh=100, fwhm=2, skyrad=[8,12], 
          pixscale=0.5, display=None, plot=False,max=max,red=None) :
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
            a=fits.open(file)[0]
        else :
            a=red.rd(file)
        im=a.data.astype(float)
        foc[ifile]=a.header['FOCUS']
        print(file,foc[ifile])

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
            print('no sources found...')
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
            print('no sources with finite photometry')
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

    if len(hfmed) > 0 :
        # get smallest median half-flux radius
        gd=np.where(hfmed>0)[0]
        besthf=np.min(hfmed[gd])
        # convert to diameter in arcsec
        besthf *= 2*pixscale
        bestind=np.argmin(hfmed[gd])
        bestfoc=foc[gd[bestind]]
        print('bestfoc: ', bestfoc)
        try: 
            # attempt a quadratic fit around the minimum
            poly=np.polyfit(foc[gd[bestind-2:bestind+3]],hfmed[gd[bestind-2:bestind+3]],2)
            bestfitfoc = -poly[1]/2/poly[0]
            if bestfitfoc < foc[gd[bestind-2]] or bestfitfoc > foc[gd[bestind+2]] :
                print('best focus out of fit range!')
                bestfitfoc=-1
                bestfithf=-1
            else :
                bestfithf = poly[0]*bestfitfoc**2+poly[1]*bestfitfoc+poly[2]
                bestfithf *= 2*pixscale
            print('bestfoc, bestfitfoc: ', bestfitfoc,bestfithf)
        except: 
            print('focus fit failed...')
            bestfitfoc=-1
            bestfithf=-1

        # plot and/or display
        if plot :
            fig2,ax2 = plots.multi(1,1)
            plots.plotp(ax2,foc[gd],hfmed[gd],xt='Focus',yt='R(half total)',size=50)
        if display is not None :
            plots.plotp(display.plotax2,foc[gd],hfmed[gd]*pixscale*2,
                        xt='Focus',yt='R(half total)',size=50)
            if bestind>1 and bestind < len(gd)-3 :
                x=np.linspace(foc[gd[bestind-2]],foc[gd[bestind+2]],100)
                plots.plotl(display.plotax2,x,(poly[0]*x**2+poly[1]*x+poly[2])*pixscale*2)

    if plot :
        print('Hit RETURN key to close plot windows and continue: ')
        inp=input()
        plt.close(fig)
        plt.close(fig2)

    # return best fit values and minimum values
    return bestfitfoc, bestfithf,  bestfoc, besthf
