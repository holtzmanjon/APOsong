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

from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.io import fits
from astropy.time import Time

from astroquery.vo_conesearch import ConeSearch, conesearch
import astroquery.utils
astroquery.utils.suppress_vo_warnings()

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
  print('no alpaca')

# pwi4 HTTP interface
import pwi4_client

from pyvista import tv
import focus
import status

# some global variables
dataroot = None
sync_process = None
disp = None

# discovery seems to fail on 10.75.0.0, so hardcode servers
def ascom_init(svrs) :
    global D, S, T, F, Filt, C, Covers
    D, S, T, F, Filt, C, Covers = (None, None, None, None, None, None, None)
    print("Alpaca devices: ")
    if svrs is None : return
    for svr in svrs:
        print(f"  At {svr}")
        print (f"    V{management.apiversions(svr)} server")
        print (f"    {management.description(svr)['ServerName']}")
        devs = management.configureddevices(svr)
        def isconnected(dev,target) :
            try: 
                dev.Connected
                target = dev
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
                C = isconnected(Camera(svr,dev['DeviceNumber']),C)
            elif dev['DeviceType'] == 'SafetyMonitor' :
                S = isconnected(SafetyMonitor(svr,dev['DeviceNumber']),S)

    print()
    print("All ASCOM commands available through devices: ")
    print('    T : telescope commands')
    print('    C : camera commands')
    print('    F : focusser command')
    print('    Filt : filter wheel commands')
    print('    D : dome commands')
 

# Camera commands
def qck(exptime,filt='current') :
    """ (shorthand) Take exposure without saving to disk, current filter by default
    """
    expose(exptime,filt)

def expose(exptime=1.0,filt='current',bin=3,box=None,light=True,display=None,name=None,
           min=None, max=None) :
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
    if filt is not None and filt != 'current': 
        pos = np.where(np.array(Filt.Names) == filt)
        if len(pos) == 0 :
            print('no such filter')
            print('available filters: ', Filt.Names)
            return
        Filt.Position=pos[0]
    elif filt == 'current' :
        filt = Filt.Names[Filt.Position]

    C.BinX=bin
    C.BinY=bin
    if box is not None :
        C.StartX = box.xmin
        C.StartY = box.ymin
        nx = box.ncol()
        ny = box.nrow()
    else :
        C.StartX = 0
        C.StartY = 0
        nx = C.CameraXSize//bin
        ny = C.CameraYSize//bin
    C.NumX = nx
    C.NumY = ny
    t = Time.now()
    C.StartExposure(exptime,light)
    while not C.ImageReady or C.CameraState != 0:
        time.sleep(1.0)
    try : data = np.array(C.ImageArray).T
    except :
        print('error reading ImageArray')
        pdb.set_trace()
    if disp is not None :
        disp.tv(data,min=min,max=max)
        disp.fig.canvas.flush_events()

    hdu=fits.PrimaryHDU(data)
    hdu.header['DATE-OBS'] = t.fits
    hdu.header['JD'] = t.jd
    hdu.header['MJD'] = t.mjd
    hdu.header['EXPTIME'] = exptime 
    hdu.header['FILTER'] = filt 
    try : 
        hdu.header['FOCUS'] = F.Position
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
    hdu.header['CCD-TEMP'] = C.CCDTemperature
    hdu.header['XBINNING'] = C.BinX
    hdu.header['YBINNING'] = C.BinY
    if light: hdu.header['IMAGTYP'] = 'LIGHT'
    else: hdu.header['IMAGTYP'] = 'DARK'

    if name is not None :
        y,m,d,hr,mi,se = t.ymdhms
        dirname = '{:s}/UT{:d}{:02d}{:02d}'.format(dataroot,y-2000,m,d)
        try: os.mkdir(dirname)
        except : pass
        files = glob.glob(dirname+'/*.fits')
        exts = []
        for f in files :
            exts.append(int(f.split('.')[-2]))
        if len(exts) > 0 : ext = np.array(exts).max() + 1
        else : ext=1
        outname = '{:s}/{:s}.{:04d}.fits'.format(dirname,name,ext)
        hdu.writeto(outname)
        return hdu, outname
    else :
        return hdu

def settemp(temp) :
    """ Change detector temperature set point and turn cooler on
    """
    C.SetCCDTemperature = temp
    C.CoolerOn = True

def cooler(state=True) :
    """ Set detector cooler state on/off
    """
    C.CoolerOn = state

def focrun(cent,step,n,exptime=1.0,filt='V',bin=3,box=None,display=None,
           max=30000, thresh=25) :
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

    files=[]
    images=[]
    for i,focval in enumerate(
             np.arange(int(cent)-n//2*int(step),int(cent)+n//2*int(step)+1,int(step))) :
        F.Move(int(focval))
        if i==0 : time.sleep(5)
        while F.IsMoving :
            time.sleep(1)
        print('position: ',F.Position,F.IsMoving)
        hdu,name = expose(exptime,filt,box=box,bin=bin,display=display,
                          max=max,name='focus_{:d}'.format(focval))
        if i == 0 :
            nr,nc=hdu.data.shape
            mosaic = np.zeros([nr,n*nc])
        mosaic[:,i*nc:(i+1)*nc] = hdu.data
        files.append(name)
        images.append(hdu)
    plt.figure()
    plt.imshow(mosaic,vmin=0,vmax=max,cmap='gray')
    plt.axis('off')
    pixscale=206265/6000*C.PixelSizeX*1.e-3*bin
    bestfitfoc, bestfithf,  bestfoc, besthf = focus.focus(files,pixscale=pixscale,
                                                display=display,max=max,thresh=thresh)
    if bestfitfoc > 0 :
        print('setting focus to best fit focus : {:.1f} with hf diameter {:.2f}',
              bestfitfoc,bestfithf)
        foc(int(bestfitfoc))
    else :
        print('setting focus to minimum image focus : {:.1f} with hf diameter {:.2f}',
              bestfitfoc,besthf)
        foc(int(bestfoc))
    return images,files

def slew(ra, dec) :
    """ Slew to RA/DEC

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
    while T.Slewing :
        time.sleep(1)
    T.Tracking = True
    domesync()

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

def usno(ra=None,dec=None,rad=1*u.degree,rmin=0,rmax=15,bmin=0,bmax=15,goto=True,
         cat= 'The USNO-A2.0 Catalogue 1') :
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
         astroquery conesearch catgalog to search from
    """
    if ra is None :
        stat = pwi.status()
        ra = stat.mount.ra_j2000_hours
        dec = stat.mount.dec_j2000_degs
        print(ra,dec)
        coords=SkyCoord("{:f} {:f}".format(ra,dec),unit=(u.hourangle,u.deg))
    elif (isinstance(ra,float) or isinstance(ra,int) ) and  \
         (isinstance(dec,float) or isinstance(dec,int) ) :
        coords=SkyCoord("{:f} {:f}".format(ra,dec),unit=(u.hourangle,u.deg))
    else :
        coords=SkyCoord("{:s} {:s}".format(ra,dec),unit=(u.hourangle,u.deg))
    result = conesearch.conesearch(coords,rad,catalog_db=cat,verbose=True)
    print(result)
    gd=np.where((result['Bmag']<bmax) & (result['Bmag']>bmin) &
                (result['Rmag']<rmax) & (result['Rmag']>rmin))[0]
    if len(gd) <= 0 :
        print('no USNO A2.0 stars found within specified rad and mag range')
        return
    best = result['_r'][gd].argmin()
    print(result[gd[best]])
    if goto: slew(result['RAJ2000'][gd[best]]/15., result['DEJ2000'][gd[best]])

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
    T.Park()
    D.Park()

def foc(val, relative=False) :
    """ Change focus, absolute (default) or relative

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

def domehome() :
    """ Home dome
    """
    D.FindHome()

def mirror_covers(open=False) :
    """ Open/close mirror covers
    """
    current = Covers.CoverState.value
    if open and current == 3 : return
    if not open and current == 1 : return
    altaz(T.Azimuth,85.)
    print('waiting for telescope to slew to high altitude...')
    while T.Slewing :
        time.sleep(1)
    if open and current != 3 :
        Covers.OpenCover()
    elif not open and current != 1 :
        Covers.CloseCover()

def coverstate() :
    print(Covers.CoverState.value)

def domeopen(dome=True,covers=True) :
    """ Open dome and mirror covers
    """
    if dome : D.OpenShutter()
    if covers : 
        # Wait for shutter open before opening mirror covers
        print('waiting for shutter to open...')
        while D.ShutterStatus.name != 'shutterOpen' :
            time.sleep(1)
        mirror_covers(True) 
        print('waiting for mirror covers to open...')
        while Covers.CoverState.value != 3 :
            time.sleep(1)

def domeclose(dome=True,covers=True) :
    """ Close mirror covers and dome
    """
    if covers : mirror_covers(False)
    if dome : 
        print('waiting 20 seconds for mirror covers to close...')
        time.sleep(20)
        # don't wait for mirror covers to report closed, in case they don't!
        try: park()
        except : pass
        D.CloseShutter()

def domesync(dosync=True) :

    global sync_process, run_sync
    def sync(update=5) :
        while run_sync :
            while T.Slewing or D.Slewing : 
                 print(T.Slewing, D.Slewing)
                 time.sleep(1)
            D.SlewToAzimuth(T.Azimuth)
            time.sleep(update)

    if dosync and sync_process is None :
        print('starting dome sync')
        sync_process=threading.Thread(target=sync)
        run_sync = True
        sync_process.start()
    elif dosync == False and sync_process is not None :
        print('stopping dome sync')
        run_sync = False
        sync_process.join()
        sync_process = None

def start_status() :
    """ Start status window thread
    """
    global proc
    proc = mp.Process(target=status.status,
                 kwargs={'pwi' : pwi, 'T' : T, 'D' : D, 'F' : F, 'Filt' : Filt, 'C' : C, 'Covers': Covers})
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
    print("  park(): park telescope and dome")
    print()
    print("Camera commands")
    print("  expose(exptime,filt,**kwargs: take an exposure")
    print("  focrun(cent,step,nsteps,exptime,filt,**kwargs): take series of exposures at different focus positions")
    print("  settemp(temp): set camera temperature set point")
    print("  cooler(state): set camera cooler state on (True) or off (False)")
    print()
    print("Status commands")
    print("  start_status(): start status window")
    print("  stop_status(): stop status window ")
    print()
    print("Use help(command) for more details")

def init() :
    """ Start ascom and pwi connections and pyvista display
    """
    global disp, dataroot, pwi_srv
    try :
        with open('aposong.yml','r') as config_file :
            config = yaml.safe_load(config_file) 
    except:
        print('no configuration file found')
        config={}
        config['devices']={}
        config['devices']['ascom_search'] = False
        config['devices']['ascom_svrs'] = []
        config['devices']['pwi_srv'] = None

    if config['devices']['ascom_search'] :
        svrs=discovery.search_ipv4(timeout=30,numquery=3)
    else :
        svrs=config['devices']['ascom_srvs']
    dataroot=config['dataroot']
    ascom_init(svrs)
    print('pwi_init...')
    pwi_srv = config['devices']['pwi_srv']
    pwi_init(pwi_srv)
    print('start_status...')
    start_status()
    disp=tv.TV(figsize=(8,6))
    commands()

def pwi_init(pwi_srv) :
    global pwi
    print('Starting PWI4 client...')
    if pwi_srv is not None :
        pwi=pwi4_client.PWI4(host=pwi_srv)
        for axis in [0,1] : pwi.mount_enable(axis)
    else :
        pwi = None

init()

