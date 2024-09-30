#interactive plots
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

import numpy as np
import os
import pdb
import glob
import yaml
import time
from datetime import datetime
import multiprocessing as mp
import threading
import yaml
import remote
import subprocess

import logging
import guiding
import yaml
import logging.config
try :
    with open('logging.yml', 'rt') as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)
except:
    print("can't open logging.yml")

from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table

from astroquery.vo_conesearch import ConeSearch, conesearch
from astroquery.vizier import Vizier
import astroquery.utils
astroquery.utils.suppress_vo_warnings()

import eshel

# alpaca imports, put in try/except for readthedocs
try:
  from alpaca import discovery, management
  from alpaca.telescope import *
  from alpaca.covercalibrator import *
  from alpaca.dome import *
  from alpaca.safetymonitor import *
  from alpaca.focuser import *
  from alpaca.filterwheel import *
  from alpaca.camera import *
except:
  logger.exception('no alpaca')

# pwi4 HTTP interface
import pwi4_client

from pyvista import tv, skycalc, centroid, stars, image
import focus
import status
import database

from collections import namedtuple
Exposure = namedtuple('Exposure', ['hdu', 'name', 'exptime', 'filter'])

# some global variables
dataroot = None
sync_process = None
guide_process = None
disp = None

# discovery seems to fail on 10.75.0.0, so hardcode servers
def ascom_init(svrs) :
    global D, S, T, F, Filt, C, Covers
    D, S, T, F, Filt, C, Covers = (None, None, None, None, None, [], None)
    print("Alpaca devices: ")
    if svrs is None : return
    for svr in svrs:
        print(f"  At {svr}")
        print (f"    V{management.apiversions(svr)} server")
        print (f"    {management.description(svr)['ServerName']}")
        devs = management.configureddevices(svr)
        def isconnected(dev,target,append=False) :
            try: 
                dev.Connected
                if append : target.append(dev)
                else :target = dev
            except :
                pass
            return target
                
        # open Alpaca devices
        for dev in devs:
            print(f"      {dev['DeviceType']}[{dev['DeviceNumber']}]: {dev['DeviceName']}")
            if dev['DeviceType'] == 'Telescope' :
                T =isconnected(Telescope(svr,dev['DeviceNumber']),T)
            elif dev['DeviceType'] == 'Dome' :
                D = isconnected(Dome(svr,dev['DeviceNumber']),D)
            elif dev['DeviceType'] == 'CoverCalibrator' :
                Covers = isconnected(CoverCalibrator(svr,dev['DeviceNumber']),Covers)
            elif dev['DeviceType'] == 'Focuser' :
                F = isconnected(Focuser(svr,dev['DeviceNumber']),F)
            elif dev['DeviceType'] == 'FilterWheel' :
                Filt = isconnected(FilterWheel(svr,dev['DeviceNumber']),Filt)
            elif dev['DeviceType'] == 'Camera' :
                C = isconnected(Camera(svr,dev['DeviceNumber']),C,append=True)
            elif dev['DeviceType'] == 'SafetyMonitor' :
                S = isconnected(SafetyMonitor(svr,dev['DeviceNumber']),S)

    C[getcam(0)].Magnification=1.33
    print()
    print("All ASCOM commands available through devices: ")
    print('    T : telescope commands')
    print('    C : camera commands')
    print('    F : focusser command')
    print('    Filt : filter wheel commands')
    print('    D : dome commands')

def getcam(camera) :
    """ Get correct list index for specified camera (ASCOM doesn't always deliver them in order!)
    """
    for index,c in enumerate(C) :
        if c.device_number == camera : 
            return index 
    print('no such camera!')
    return -1

# Camera commands
def qck(exptime,filt='current') :
    """ (shorthand) Take exposure without saving to disk, current filter by default
    """
    expose(exptime,filt)

def expose(exptime=1.0,filt='current',bin=3,box=None,light=True,display=None,name=None,
           min=None, max=None, cam=0, insert=True) :
    """ Take an exposure with camera

    Parameters
    ----------
    exptime : float
           Exposure time in seconds
    filt : str
           Name of filter
    bin : int, default=3
           binning factor (both x and y)
    box : pyvista BOX, default=None
           if specified, window to box parameters 
    light : bool, default=True
           open shuter for exposure?
    display : pyvista tv, default=None
           pyvista display tool to display into if specified
    name  : str, default=None
           root file to save image to 

    Returns
    -------
           HDU of image [,output file name if name!=None]
    """
    exposure = Exposure(None,None,None,None)
    if filt is not None and filt != 'current': 
        pos = np.where(np.array(Filt.Names) == filt)
        if len(pos) == 0 :
            logger.warning('no such filter')
            logger.warning('available filters: ', Filt.Names)
            return exposure
        Filt.Position=pos[0]
    elif filt == 'current' :
        filt = Filt.Names[Filt.Position]

    try :
        icam = getcam(cam) 
        C[icam].BinX=bin
        C[icam].BinY=bin
        if box is not None :
            C[icam].StartX = box.xmin
            C[icam].StartY = box.ymin
            nx = box.ncol()
            ny = box.nrow()
        else :
            C[icam].StartX = 0
            C[icam].StartY = 0
            nx = C[icam].CameraXSize//bin
            ny = C[icam].CameraYSize//bin
        C[icam].NumX = nx
        C[icam].NumY = ny
        t = Time.now()
        C[icam].StartExposure(exptime,light)
        for i in range(int(exptime),0,-1) :
            print('{:<4d}'.format(i),end='\r')
            time.sleep(0.99)
        while not C[icam].ImageReady or C[icam].CameraState != 0:
            time.sleep(1.0)
        data = np.array(C[icam].ImageArray).T
    except :
        logger.exception('ERROR : exposure failed')
        return exposure

    if display is not None :
        display.tv(data,min=min,max=max)
        display.fig.canvas.flush_events()

    hdu=fits.PrimaryHDU(data)
    hdu.header['DATE-OBS'] = t.fits
    hdu.header['JD'] = t.jd
    hdu.header['MJD'] = t.mjd
    hdu.header['EXPTIME'] = exptime 
    hdu.header['FILTER'] = filt 
    try : 
        hdu.header['FOCUS'] = getfoc()
    except : pass
    try :
        stat = pwi.status()
        hdu.header['RA'] = stat.mount.ra_j2000_hours
        hdu.header['DEC'] = stat.mount.dec_j2000_degs
        hdu.header['AZ'] = stat.mount.azimuth_degs
        hdu.header['ALT'] = stat.mount.altitude_degs
        hdu.header['ROT'] = stat.rotator.mech_position_degs
        hdu.header['TELESCOP'] = 'APO SONG 1m'
    except : pass
    hdu.header['CCD-TEMP'] = C[icam].CCDTemperature
    hdu.header['XBINNING'] = C[icam].BinX
    hdu.header['YBINNING'] = C[icam].BinY
    if light: hdu.header['IMAGTYP'] = 'LIGHT'
    else: hdu.header['IMAGTYP'] = 'DARK'

    tab=Table()
    cards = ['DATE-OBS','MJD','EXPTIME','FILTER','FOCUS','CCD-TEMP','XBINNING','YBINNING','RA','DEC','AZ','ALT','ROT'] 
    cols = ['dateobs','mjd','exptime','filter','focus','ccdtemp','xbin','ybin','ra','dec','az','alt','rot'] 
    for card,col in zip(cards,cols) :
        tab[col] = [hdu.header[card]]
    tab['camera'] = [cam]

    if name is not None :
        y,m,d,hr,mi,se = t.ymdhms
        dirname = os.path.dirname('{:s}/UT{:d}{:02d}{:02d}/{:s}'.format(dataroot,y-2000,m,d,name))
        try: os.makedirs(dirname)
        except : pass
        files = glob.glob(dirname+'/*.fits')
        exts = []
        for f in files :
            exts.append(int(f.split('.')[-2]))
        if len(exts) > 0 : ext = np.array(exts).max() + 1
        else : ext=1
        outname = '{:s}/{:s}.{:04d}.fits'.format(dirname,os.path.basename(name),ext)
        hdu.writeto(outname)
        exposure = Exposure(hdu, outname, exptime, filt)
        tmp=outname.split('/')
        tab['file'] = [tmp[-2]+'/'+tmp[-1]]
    else :
        outname=''
        exposure = Exposure(hdu, None, exptime, filt)
        tab['file'] = ''

    if insert :
        d=database.DBSession()
        d.ingest('obs.exposure',tab,onconflict='update')
        d.close()

    return exposure

def settemp(temp,cam=0) :
    """ Change detector temperature set point and turn cooler on
    """
    icam=getcam(cam)
    C[icam].SetCCDTemperature = temp
    C[icam].CoolerOn = True

def cooler(state=True,cam=0) :
    """ Set detector cooler state on/off
    """
    icam=getcam(cam)
    C[icam].CoolerOn = state

def focrun(cent,step,n,exptime=1.0,filt='V',bin=3,box=None,display=None,
           max=30000, thresh=25,cam=0,plot=False) :
    """ Obtain a focus run

    Parameters
    ----------
    cent : int
           Middle focus value
    step : int
           Focus step size
    n : int
           Number of steps to take, centered on middle value
    exptime : float
           Exposure time
    filt : str
           Filter to use
    bin : int, default=3
           Binning factor in x and y
    box : pyvista BOX, default=None
           if specified, window to box parameters 
    disp : pyvista tv, default=None
           If given, display each image as it is being taken into specified device

    Returns
    -------
           List of file names taken for focus run
    """

    foc0=getfoc()
    files=[]
    names=[]
    images=[]
    focvals= np.arange(int(cent)-n//2*int(step),int(cent)+n//2*int(step)+1,int(step))
    for i,focval in enumerate(focvals) :
        foc(int(focval))
        if i==0 : time.sleep(5)
        while True :
            try :
                moving=F.IsMoving
                if not moving : break
            except : 
                logger.exception('error with F.IsMoving')
            time.sleep(2)
        
        logger.info('position: {:d} {}'.format(getfoc(),moving))
        exp = expose(exptime,filt,box=box,bin=bin,display=display,
                          max=max,name='focus_{:d}'.format(focval),cam=cam)
        hdu=exp.hdu
        tmp=exp.name.split('/')
        name=tmp[-2]+'/'+tmp[-1]
        names.append(name)
        files.append(exp.name)
        images.append(hdu)
        #if i == 0 :
        #    nr,nc=hdu.data.shape
        #    mosaic = np.zeros([nr,n*nc])
        #mosaic[:,i*nc:(i+1)*nc] = hdu.data
    #plt.figure()
    #plt.imshow(mosaic,vmin=0,vmax=max,cmap='gray')
    #plt.axis('off')
    try :
        bestfitfoc, bestfithf,  bestfoc, besthf = focus.focus(files,pixscale=pixscale(cam,bin=bin),
                                                display=display,max=max,thresh=thresh,plot=plot)
        if bestfitfoc > 0 :
            logger.info('setting focus to best fit focus : {:.1f} with hf diameter {:.2f}'.format(
                  bestfitfoc,bestfithf))
            f=foc(int(bestfitfoc))
        else :
            logger.info('setting focus to minimum image focus : {:.1f} with hf diameter {:.2f}'.format(
              bestfoc,besthf))
            f=foc(int(bestfoc))
    except :
        bestfitfoc, bestfithf,  bestfoc, besthf = -1, -1, -1, -1
        logger.exception('focus failed')
        f=foc(foc0)

    tab=Table()
    tab['mjd'] = [Time.now().mjd]
    tab['exptime'] = [exptime]
    tab['filter'] = [filt]
    tab['bin'] = [bin]
    tab['camera'] = [cam]
    tab['focvals'] = [focvals]
    tab['bestfoc'] = [bestfoc]
    tab['besthf'] = [besthf]
    tab['bestfitfoc'] = [bestfitfoc]
    tab['bestfithf'] = [bestfithf]
    tab['files'] = [names]
    d=database.DBSession()
    d.ingest('obs.focus',tab,onconflict='update')
    d.close()

    return f

def pixscale(cam=0,bin=1) :
    """ Return pixscale for desired camera
    """
    icam=getcam(cam)
    scale=206265/6000*C[icam].PixelSizeX*1.e-3*bin
    try: scale /= C[icam].Magnification
    except: pass
    return scale

def slew(ra, dec,dome=True) :
    """ Slew to RA/DEC and wait for telescope/dome

    Parameters
    ----------
    ra : float or str
         RA in degrees (float), or hh:mm:ss (str)
    dec : float or str
         DEC in degrees (float), or dd:mm:ss (str)
    """
    domesync(False)
    if (isinstance(ra,float) or isinstance(ra,int) ) and  \
       (isinstance(dec,float) or isinstance(dec,int) ) :
        pwi.mount_goto_ra_dec_j2000(ra,dec)
    elif isinstance(ra,str) and isinstance(dec,str) :
        coords=SkyCoord("{:s} {:s}".format(ra,dec),unit=(u.hourangle,u.deg))
        pwi.mount_goto_ra_dec_j2000(coords.ra.value/15.,coords.dec.value)

    #tra, tdec = j2000totopocentric(ra,dec) 
    #T.SlewToCoordinatesAsync(tra,tdec)
    time.sleep(5)
    while T.Slewing :
        time.sleep(1)
    T.Tracking = True

    # Wait for dome slewing to stop
    domesync(dome)
    while True :
        if not T.Slewing :
            time.sleep(10)
            if not D.Slewing : break

def altaz(az,alt) :
    """ Slew to specified az / alt
    Parameters
    ----------
    az : float
         az in degrees 
    alt : float
         alt in degrees
    """
    pwi.mount_goto_alt_az(alt, az)

def usno(ra=None,dec=None,rad=1*u.degree,mag='Rmag',magmin=0,magmax=20,goto=True,
         cat= 'I/252') :
    """ Find/goto nearest USNO stars from specified position (default current telescope) and mag range

    Parameters
    ----------
    ra : str or float, default=None
         RA to search around, default current telescope position
    dec : str or float, default=None
         DEC to search around, default current telescope position
    rad : angular quantity, default = 1*u.degree
         Radius to use for sear4ch
    bmin, bmax : float, default=(0,15)
         Minimum, maximum Bmag
    rmin, rmax : float, default=(0,15)
         Minimum, maximum Rmag
    goto : bool, default=True
         If True, slew to closest returned object
    cat : str, default= 'The USNO-A2.0 Catalogue 1'
         astroquery conesearch catalog to search from
         hr=V/50,  sao=I/131A/sao (but doesn't have J2000)
    """
    if ra is None :
        stat = pwi.status()
        ra = stat.mount.ra_j2000_hours
        dec = stat.mount.dec_j2000_degs
        coords=SkyCoord("{:f} {:f}".format(ra,dec),unit=(u.hourangle,u.deg))
    elif (isinstance(ra,float) or isinstance(ra,int) ) and  \
         (isinstance(dec,float) or isinstance(dec,int) ) :
        coords=SkyCoord("{:f} {:f}".format(ra,dec),unit=(u.hourangle,u.deg))
    else :
        coords=SkyCoord("{:s} {:s}".format(ra,dec),unit=(u.hourangle,u.deg))

    try :
      viz=Vizier(columns=['*','+_r'],catalog=cat)
      if mag is not None :
        result=viz.query_region(coords,radius=rad,
                column_filters={mag:'<{:f}'.format(magmax)})[0]
        gd=np.where(result[mag]>magmin)[0]
        result=result[gd]
      else :
        result=viz.query_region(coords,radius=rad)[0]
    except :
        logger.warning('usno: no stars found within specified rad and mag range')
        return

    best = result['_r'].argmin()
    logger.info(result[best])
    if goto: 
        try: 
            slew(result['RAJ2000'][best]/15., result['DEJ2000'][best])
        except : 
            slew(result['RA_ICRS'][best]/15., result['DE_ICRS'][best])

def center(x0=None,y0=None,exptime=5,bin=1,filt=None,settle=3,cam=0) :
    """ Center star
    """
    exp=expose(exptime,filt=filt,bin=bin)
    disp.tv(exp.hdu)
    print('mark star on display: ')
    k,x,y=disp.tvmark()
    icam=getcam(cam)
    if x0 is None : x0=C[icam].NumX//2
    if y0 is None : y0=C[icam].NumY//2
    offsetxy((x-x0),(y-y0),scale=pixscale(),pa=rotator()-T.Altitude+100)
    time.sleep(settle)

def guide(start=True,x0=774,y0=466,rad=25,exptime=5,bin=1,filt=None,data=None,maskrad=7,
          thresh=100,fwhm=1.5,
          display=None,prop=0.7,settle=3,vmax=10000,rasym=True,name='guide',inter=False) :
    """ Start guiding

        Parameters
        ----------
        start   bool default=True
                Start or stop guiding
        x0, y0  float, default=(646,494)
                Location to guide to
        rad     int, default=25
                radius to use for centering algorithm
        exptime float, default=5
                guide exposure time
        bin     integer default=1
                guider binning factor
        filt    str, default=None
                filter to use
        display default=None
                if pyvista TV, display guider images, and don't run in thread
        vmax    float, default=10000
                maximum value for display
        prop    float, default=0.7
                coefficient for proportional term 
        settle  float, default=3
                time (seconds) to wait after offset before starting next image
        rasym   bool, default=True
                use radial asymmetry minimum centroider, otherwise marginal_gfit
        data    array-like default=None
                for testing, supplied image rather than acquire one
    """

    global guide_process

    if start and guide_process is None :
        # refine position of hole
        eshel.lamps(mirror=True,quartz=True)
        time.sleep(3)
        a=expose(0.01,bin=1,filt=None,cam=0,max=50000,display=disp,name='guide/hole')
        time.sleep(3)
        eshel.lamps(close=False)
        center=centroid.rasym_centroid(a.hdu.data,x0,y0,rad=12)
        x0=center.x
        y0=center.y
        logger.info('hole center : {:.2f} {:.2f}'.format(x0,y0))

        if data is None :
            exp=expose(exptime,filt=filt,bin=bin,display=disp,name='guide/acquire')
            hdu=exp.hdu
        else :
            hdu = data
        if inter :
            disp.tv(hdu,max=vmax)
            print('mark star on display: ')
            k,x,y=disp.tvmark()
        else :
            # find brightest star automatically
            # adjust exposure time to get good counts
            niter=0
            while True :
                mad=np.nanmedian(np.abs(hdu.data-np.nanmedian(hdu.data)))
                objs=stars.find(hdu.data,thresh=thresh*mad,fwhm=fwhm/pixscale(),brightest=1)
                if objs is None : raise RuntimeError('no objects found')
                peak=objs[0]['peak']
                if (peak < 60000 and peak>5000 ) or niter>10 or exptime>9.99: break
                if peak>60000 : exptime*=0.8
                else: exptime = np.max([0.01,np.min([exptime*60000/peak,10])])
                logger.info('new exptime: {:.2f}'.format(exptime))
                exp=expose(exptime,filt=filt,bin=bin,display=disp)
                hdu=exp.hdu
                niter+=1

            x=objs[0]['x']
            y=objs[0]['y']

        exp=expose(exptime,display=disp,filt=filt,bin=bin,name='guide/acquire')
        mad=np.nanmedian(np.abs(hdu.data-np.nanmedian(hdu.data)))
        objs=stars.find(hdu.data,thresh=thresh*mad,fwhm=fwhm/pixscale(),brightest=1)
        x=objs[0]['x']
        y=objs[0]['y']
        logging.info('offsetting {:.1f} {:.1f}'.format(x-x0,y-y0))
        offsetxy((x-x0),(y-y0),scale=pixscale())
        time.sleep(settle)

        # create mask
        mask=np.zeros_like(hdu.data)
        yg,xg=np.mgrid[0:hdu.data.shape[0],0:hdu.data.shape[1]]
        r2=(xg-x0)**2+(yg-y0)**2
        bd=np.where(r2<maskrad**2)
        mask[bd]=1

        # subwindow for guiding
        box=image.BOX(cr=int(y0),cc=int(x0),n=int(7*rad))
        x0-=box.xmin
        y0-=box.ymin
        mask=image.window(mask,box=box)
        exp=expose(exptime,display=disp,filt=filt,bin=bin,box=box,name='guide/acquire')
        center=centroid.rasym_centroid(exp.hdu.data,x0,y0,rad,mask=mask,skyrad=[35,40],plot=disp)
        if center.x>0 and center.tot>10000:
            x = center.x
            y = center.y
        else :
            objs=stars.find(exp.hdu.data,thresh=thresh*mad,fwhm=fwhm/pixscale(),brightest=1)
            x=objs[0]['x']
            y=objs[0]['y']
        logging.info('offsetting {:.1f} {:.1f}'.format(x-x0,y-y0))
        offsetxy((x-x0),(y-y0),scale=pixscale())
        time.sleep(settle)

        if inter :
            disp.tv(mask)
            disp.tvcirc(x0,y0,rad=rad)
            pdb.set_trace()

        guiding.run_guide = True
        # navg no greater than 5, no less than 1, otherwise 5/exptime
        navg=np.min([5,np.max([1,int(5/exptime)])])
        logger.info('starting guiding: {:.2f} {:.2f} exptime: {:.2f} navg: {:d} prop: {:.2f}'.format(
                     x0,y0,exptime, navg, prop))
        if display is None :
            guide_process=threading.Thread(target=guiding.doguide,args=(x0,y0),
                kwargs={'exptime' :exptime,'navg' :navg,'mask': mask ,'disp' :None,
                        'filt' : filt, 'vmax': vmax, 'box': box, 'prop' : prop})
            guide_process.start()
        else :
            guiding.doguide(exptime=exptime,navg=navg,x0=x0,y0=y0,mask=mask,disp=disp,filt=filt,vmax=vmax,box=box,prop=prop)
    elif not start and guide_process is not None :
        logger.info('stopping guiding')
        guiding.run_guide = False
        guide_process.join()
        guide_process = None


def offsetxy(dx,dy,sign=-1,scale=0.16,pa=None) :
    """
    Offset in detector coordinates
    """
    pa = sign*rotator()*np.pi/180.
    dra =  -dx*np.cos(pa) - dy*np.sin(pa)
    ddec = -dx*np.sin(pa) + dy*np.cos(pa) 
    offset(dra*scale,ddec*scale)
    #offset(dra*scale/np.cos(T.Declination*np.pi/180.),ddec*scale)

def offset(dra, ddec) :
    """
    One or more of the following offsets can be specified as a keyword argument:

    AXIS_reset: Clear all position and rate offsets for this axis. Set this to any value to issue the command.
    AXIS_stop_rate: Set any active offset rate to zero. Set this to any value to issue the command.
    AXIS_add_arcsec: Increase the current position offset by the specified amount
    AXIS_set_rate_arcsec_per_sec: Continually increase the offset at the specified rate

    As of PWI 4.0.11 Beta 7, the following options are also supported:
    AXIS_stop: Stop both the offset rate and any gradually-applied commands
    AXIS_stop_gradual_offset: Stop only the gradually-applied offset, and maintain the current rate
    AXIS_set_total_arcsec: Set the total accumulated offset at the time the command is received to the specified value. Any in-progress rates or gradual offsets will continue to be applied on top of this.
    AXIS_add_gradual_offset_arcsec: Gradually add the specified value to the total accumulated offset. Must be paired with AXIS_gradual_offset_rate or AXIS_gradual_offset_seconds to determine the timeframe over which the gradual offset is applied.
    AXIS_gradual_offset_rate: Paired with AXIS_add_gradual_offset_arcsec; Specifies the rate at which a gradual offset should be applied. For example, if an offset of 10 arcseconds is to be applied at a rate of 2 arcsec/sec, then it will take 5 seconds for the offset to be applied.
    AXIS_gradual_offset_seconds: Paired with AXIS_add_gradual_offset_arcsec; Specifies the time it should take to apply the gradual offset. For example, if an offset of 10 arcseconds is to be applied over a period of 2 seconds, then the offset will be increasing at a rate of 5 arcsec/sec.

    Where AXIS can be one of:

    ra: Offset the target Right Ascension coordinate
    dec: Offset the target Declination coordinate
    axis0: Offset the mount's primary axis position 
           (roughly Azimuth on an Alt-Az mount, or RA on In equatorial mount)
    axis1: Offset the mount's secondary axis position 
           (roughly Altitude on an Alt-Az mount, or Dec on an equatorial mount)
    path: Offset along the direction of travel for a moving target
    transverse: Offset perpendicular to the direction of travel for a moving target
    """
    pwi.mount_offset(ra_add_arcsec=dra)
    pwi.mount_offset(dec_add_arcsec=ddec)

def j2000totopocentric(ra,dec) :
    """ Routine to convert J2000 coordinates to topocentric RA/DEC
    """
    #from coordio import sky, site, time
    coords=sky.ICRS(np.array([[15*ra,dec]]))
    s=site.Site('APO')
    s.set_time(time.Time())    # defaults to now
    obs=sky.Observed(coords,site=s)
    print('topo: ',obs.ra/15.,obs.dec)
    return obs.ra/15., obs.dec

    # tpm routines
    #apo=EarthLocation.of_site('APO')
    #d=datetime.now()
    #utc = tpm.gcal2j(d.year,d.month,d.day)
    #tt = tpm.utc3tdb(utc)
    #v6 = convert.cat2v6(ra*np.pi/180,dec*np.pi/180)
    #v6_app = convert.convertv6(s1=6,s2=16,lon=apo.lon.value,lat=apo.lat.value,alt=apo.height.value)

def rotator(offset=100.) :
    """ Get current rotator position angle, from parallactic angle and altitude for port 2  and telescope status for port 1

    Parameters
    ----------
    offset : float, default=100
             Angle constant (degrees) that needs to be set depending on camera/fixed rotator position for port 2
    """
    stat = pwi.status()
    if stat.m3.port == 2 :
        t=Time(Time.now(),location=EarthLocation.of_site('APO'))
        lst=t.sidereal_time('mean').value
        ha = lst - T.RightAscension
        pa = skycalc.parang(ha,T.Declination,'APO')-T.Altitude+offset
    else :
        inst = stat.rotator.mech_position_degs
        pa = stat.rotator.field_angle_degs

    return pa 

def tracking(tracking) :
    """ Set telescope tracking on (True) or off (False)

    Parameters
    ----------
    tracking : bool
               if True, turn tracking on, if False, turn tracking off
    """
    T.Tracking= tracking

def park() :
    """ Park telescope and dome
    """
    domesync(False)
    try: T.Park()
    except: logger.error('telescope.Park raised an exception')
    D.Park()

def foc(val, relative=False) :
    """ Change focus, depending on port
    """
    stat = pwi.status()
    if stat.m3.port == 1 :
        return pfoc(val, relative=relative)
    else :
        return zfoc(val, relative=relative)

def getfoc() :
    """ Get focus, depending on port
    """
    stat = pwi.status()
    if stat.m3.port == 1 :
        return F.Position
    else :
        return zfoc()

def pfoc(val, relative=False) :
    """ Change Planewave port 1 focus, absolute (default) or relative

    Parameters
    ----------
    val : int
          new focus value, absolute value unless relative=True
    relative : bool, default=False
          if True, then move amount relative to current position
    """
    if relative :
        val += F.Position
    F.Move(val)
    return val

def zfoc(val=None, relative=False) :
    """ Change camera focus, absolute (default) or relative
    """
    if val is None :
        val= remote.client(remote_srv,'zaber position') 
        return int(float(val)*1000)
    if relative :
        val += int(float(remote.client(remote_srv,'zaber position') )*1000)
    val = remote.client(remote_srv,'zaber position {:f}'.format(val/1000.)) 
    return int(float(val)*1000)

def iodine_tset(val=None) :
    """ Get/set iodine cell set temperature
    """
    if val is not None :
        tset = remote.client(remote_srv,'tc300 tset {:d}'.format(val*10))
    else :
        tset = remote.client(remote_srv,'tc300 tset')
    return tset

def iodine_tget() :
    """ Get/set iodine cell actual temperature
    """
    tact = remote.client(remote_srv,'tc300 tact')
    return tact

def iodine_in(val=21.,focoffset=-4625) :
    """ Move iodine cell into beam
    """
    iodine_position(val)
    foc(focoffset,relative=True)

def iodine_out(val=141.,focoffset=4625) :
    """ Move iodine cell out of beam
    """
    iodine_position(val)
    foc(focoffset,relative=True)

def iodine_position(val=None) :
    """ Get/set iodine stage position
    """
    if val is not None :
        pos = remote.client(remote_srv,'lts position {:f}'.format(val))
    else :
        pos = remote.client(remote_srv,'lts position')
    return pos

def domehome() :
    """ Home dome
    """
    D.FindHome()

def fans_on(roles=None):
    """
    roles: if None, turn on all fans
    Otherwise, can be a CSV string of one or more fan roles to turn on:
        m1rear: Primary mirror fans (rear fans only)
        m1side: Primary mirror fans (side fans only)
        m2: Secondary mirror fans
        m3: M3 mirror fans
        m1heaters: Primary mirror heat distribution fans
    """
    logger.info("Turning fans on ...")
    pwi.fans_on(roles)

def fans_off(roles=None):
    """ Turn off PWI fans
    """
    logger.info("Turning fans off ...")
    pwi.fans_off(roles)

def port(port) :
    """ Change PWI M3 port
    """
    pwi.m3_goto(port)

def mirror_covers(open=False) :
    """ Open/close mirror covers
    """
    current = Covers.CoverState.value
    if open and current == 3 : return
    if not open and current == 1 : return
    altaz(T.Azimuth,85.)
    logger.info('waiting for telescope to slew to high altitude...')
    while T.Slewing :
        time.sleep(1)
    if open and current != 3 :
        Covers.OpenCover()
    elif not open and current != 1 :
        Covers.CloseCover()
        logger.info('waiting 20 seconds for mirror covers to close...')
        time.sleep(20)

def louvers(open=False) :
    """ Open louvers
    """
    if open :
        subprocess.run("ms on 9",shell=True)
    else :
        subprocess.run("ms off 9",shell=True)
    subprocess.run("ms list 9",shell=True)

def coverstate() :
    print(Covers.CoverState.value)

def domeopen(dome=True,covers=True,fans=True,louvers=False) :
    """ Open dome and mirror covers
    """
    if dome : D.OpenShutter()
    if covers : 
        # Wait for shutter open before opening mirror covers
        logger.info('waiting for shutter to open...')
        while D.ShutterStatus.name != 'shutterOpen' :
            time.sleep(1)
        mirror_covers(True) 
        logger.info('waiting for mirror covers to open...')
        while Covers.CoverState.value != 3 :
            time.sleep(1)
    if fans :
        fans_on()

def domeclose(dome=True,covers=True,fans=True, closelouvers=True) :
    """ Close mirror covers and dome
    """
    guide(False)
    if closelouvers : louvers(False)
    if fans : fans_off()
    if covers : 
        mirror_covers(False)
    if dome : 
        # don't wait for mirror covers to report closed, in case they don't!
        park()
        D.CloseShutter()

def domesync(dosync=True,manual=False) :
    """ Start/stop domesync thread
    """
    global sync_process, run_sync
    def sync(update=5) :
        while run_sync :
            while T.Slewing or D.Slewing : 
                 time.sleep(1)
            # Get telescope az 3x and median to avoid glitches
            az=[]
            for i in range(3) :
                az.append(T.Azimuth)
                time.sleep(0.2)

            D.SlewToAzimuth(np.median(az))
            time.sleep(update)

    if manual :
        print('Rotate dome manually....then continue')
        pdb.set_trace()
        return
    if dosync and sync_process is None :
        logger.info('starting dome sync')
        sync_process=threading.Thread(target=sync)
        run_sync = True
        sync_process.start()
    elif dosync == False and sync_process is not None :
        logger.info('stopping dome sync')
        run_sync = False
        sync_process.join()
        sync_process = None

def start_status(camera=True) :
    """ Start status window thread
    """
    global proc
    # allow for no camera state to avoid command collisions
    if camera :
        Cam = C
    else : 
        Cam = None
    proc = mp.Process(target=status.status,
                 kwargs={'pwi' : pwi, 'T' : T, 'D' : D, 'F' : F, 'Filt' : Filt, 'C' : Cam, 'Covers': Covers})
    proc.daemon = True
    proc.start()

def stop_status() :
    """ Stop status window thread
    """
    proc.terminate()

def commands() :
    print()
    print("Observatory commands")
    print("  domeopen(): open dome, mirror covers")
    print("  domeclose(): close mirror covers, dome, and park")
    print()
    print("Dome commands")
    print("  domehome(): move dome to home position")
    print()
    print("Telescope commands")
    print("  slew(ra,dec): slew to coordinates")
    print("  offset(ra,dec): offset telescope")
    print("  usno([ra,dec]) : find/slew to USNO A2.0 star")
    print("  altaz(az,alt): slew to az/alt coordinates")
    print("  foc(focus) : set focus to specified value")
    print("  tracking(True|False): turn tracking on/off")
    print("  mirror_covers(True|False): control mirror covers")
    print("  port(port): move tertiary to requested port")
    print("  park(): park telescope and dome")
    print()
    print("Camera commands")
    print("  expose(exptime,filt,**kwargs: take an exposure")
    print("  focrun(cent,step,nsteps,exptime,filt,**kwargs): take series of exposures at different focus positions")
    print("  settemp(temp): set camera temperature set point")
    print("  cooler(state): set camera cooler state on (True) or off (False)")
    print()
    print("Iodine commands")
    print("  iodine_in : move iodine cell into beam")
    print("  iodine_out : move iodine cell out of beam")
    print("  iodine_position([val]) : get or set (with val) iodine stage position")
    print("  iodine_tset(val) : set iodine temperature (both channels)")
    print("  iodine_tget(): get actual iodine temperatures")
    print()
    print("Status commands")
    print("  start_status(): start status window")
    print("  stop_status(): stop status window ")
    print()
    print("Use help(command) for more details")

def init() :
    """ Start ascom and pwi connections and pyvista display
    """
    global disp, dataroot, pwi_srv, remote_srv
    try :
        with open('aposong.yml','r') as config_file :
            config = yaml.safe_load(config_file) 
    except:
        print('no configuration file found')
        config={}
        config['devices']={}
        config['devices']['ascom_search'] = False
        config['devices']['ascom_srvs'] = []
        config['devices']['pwi_srv'] = None
        config['dataroot'] = './'

    if config['devices']['ascom_search'] :
        svrs=discovery.search_ipv4(timeout=30,numquery=3)
    else :
        svrs=config['devices']['ascom_srvs']
    dataroot=config['dataroot']
    try : updatecamera=config['updatecamera']
    except : updatecamera=True
    ascom_init(svrs)
    print('pwi_init...')
    pwi_srv = config['devices']['pwi_srv']
    pwi_init(pwi_srv)
    try : remote_srv = config['devices']['remote_srv']
    except : remote_srv = None
    #print('start_status...')
    #start_status(updatecamera)
    disp=tv.TV(figsize=(9.5,6))
    commands()

def pwi_init(pwi_srv) :
    global pwi
    print('Starting PWI4 client...')
    if pwi_srv is not None :
        pwi=pwi4_client.PWI4(host=pwi_srv)
        for axis in [0,1] : pwi.mount_enable(axis)
    else :
        pwi = None

logger=logging.getLogger(__name__)
init()

