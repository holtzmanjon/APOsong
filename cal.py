# routines for obtaining calibration data using eShel calibration source
import aposong
import time
from holtztools import plots
from pyvista import imred, spectra
import matplotlib.pyplot as plt
import pdb

def getlamps() :
    state = ['Off','On']
    for dev in range(4):
        print(dev,aposong.SW[1].GetSwitchName(dev),state[aposong.SW[1].GetSwitch(dev)])

def lamps(close=False,mirror=False,thar=False,quartz=False,led=False) :
    """ Send commands to remote socket for eShel calibration control
    """
    for dev,state in enumerate([quartz,thar,led,mirror]) :
        if aposong.SW[1].GetSwitch(dev) != state :
            aposong.SW[1].SetSwitch(dev,state)
  
def cals(display=None,flats=2,thar=2,flat_exptime=50,thar_exptime=120,cam=3,bin=2,root='',calstage=True,mirror=False) :
    """ Take series of eShel cals, flat and ThAr

    Parameters
    ----------
    flats : int, default=2
           number of flats (quartz+led) to take
    thar : int, default=2
           number of ThAr to take
    flat_exptime: float, default=50
           flat (quartz+led) exposure time
    thar_exptime: float, default=120
           ThAr exposure time
    cam : int, default=3
           which camera to use
    bin : int, default=2
           Camera binning
    root  : str, default=''
           root prefix for output file names (prepends 'flat' and 'thar')
    calstage : bool, default=True
           move calstage in before cals, and out afterwards
    mirror : bool, default=False
           move eShel internal calibration mirror, for eShel calibration feed
    display : pyvista tv, default=None
           pyvista display tool to display into if specified
    """

    if calstage :
        foc=aposong.foc()
        aposong.calstage_find()
        aposong.calstage_in()
    lamps(mirror=mirror,quartz=True,led=True)
    time.sleep(3)
    aposong.gexp(0.001,display=display,name=root+'guide',max=20000)
    if flats>0 :
        for i in range(flats) :
            aposong.expose(flat_exptime,filt=None,bin=bin,display=display,cam=cam,name=root+'flat')
    if thar>0 :
        lamps(mirror=mirror,thar=True)
        time.sleep(3)
        for i in range(thar) :
            aposong.expose(thar_exptime,filt=None,bin=bin,display=display,cam=cam,name=root+'thar')
    lamps()
    time.sleep(3)
    if calstage :
        aposong.calstage_out()
        aposong.foc(foc)

def dark(exptime,ccdtemp=None,cam=0,n=5,bin=1,filt=None,name=None) :
    """ Take one series of darks
    """

    if ccdtemp is not None :
        aposong.settemp(ccdtemp,cam=cam)
        time.sleep(300)
    for i in range(n) :
        aposong.expose(exptime,cam=cam,light=False,bin=bin,filt=filt,name=name)

def darks(temps=[-10,-5,-15],bins=[2],exps=[30,60,120,180,240,300,600],cam=3,n=11,filt=None) :
    """ Take multiple series of darks at different temps, bins, exptimes
    """
    for temp in temps :
        aposong.settemp(temp,cam=cam)
        time.sleep(300)
        for bin in bins :
            for exp in exps :
               name='dark_cam{:d}_temp{:d}_bin{:d}_t{:d}'.format(cam,temp,bin,exp)
               dark(exp,cam=cam,n=n,filt=filt,bin=bin,name=name)


def sfocrun(foc0=427500,red=None,trace=None,wav=None,display=None,plotonly=False) :

    if red == None : red=imred.Reducer('SONG',dir='/data/1m/UT251119')
    if trace == None : trace=spectra.Trace('./UT251119_Trace_fiber2.fits')
    if wav == None : wav=spectra.WaveCal('./UT251119_WaveCal_fiber2.fits')

    if not plotonly :
        aposong.specfoc(foc0)

        aposong.calstage_find()
        aposong.calstage_in()
        lamps(quartz=True,led=True)
        s=aposong.sexp(50,name='flat',display=display,max=65535)
        t=red.reduce(s.name)
        trace.retrace(t)

        lamps(thar=True)

    fiber='specfoc'
    for i,df in enumerate(range(-2500,2501,500)) :
        if not plotonly :
            aposong.specfoc(foc0+df)
            s=aposong.sexp(50,name='{:s}_{:d}'.format(fiber,foc0+df),display=display,max=65535)
            t=red.reduce(s.name)
        else :
            t=red.reduce(267+i)
            if display is not None : 
                display.tv(t)
                display.imexam()
        tec=trace.extract(t)
        wav=spectra.WaveCal('./UT251119_WaveCal_fiber2.fits')
        wav.identify(tec)
        fig,ax=plots.multi(1,2,hspace=0.001)
        plots.plotc(ax[0],wav.waves,wav.fwhm,wav.pix,yr=[0,10],yt='FWHM',colorbar=True,zt='xpixel',cmap='viridis',size=50)
        ax[0].set_ylim(0,10)
        plots.plotc(ax[1],wav.waves,wav.waves/(wav.fwhm*(wav.wave(pixels=[wav.pix,wav.y])-wav.wave(pixels=[wav.pix+1,wav.y]))),wav.pix,
                   colorbar=True,yr=[60000,100000],xt='Wavelength',yt='R',zt='xpixel',size=50,cmap='viridis')
        ax[0].grid()
        ax[1].grid()
        fig.suptitle('specfoc: {:d}'.format(foc0+df))
        fig.savefig('{:s}_{:d}.png'.format(fiber,foc0+df))
        plt.close()
    if not plotonly: 
        lamps()
        aposong.calstage_out()
        aposong.specfoc(foc0)

