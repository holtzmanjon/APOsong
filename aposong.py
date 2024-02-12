#interactive plots
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

import numpy as np
import os
import pdb
import glob

# alpaca imports
try:
  from alpaca import discovery, management
  from alpaca.telescope import *
  from alpaca.dome import *
  from alpaca.safetymonitor import *
  from alpaca.focuser import *
  from alpaca.filterwheel import *
  from alpaca.camera import *
except:
  print('no alpaca')

# pwi4 HTTP interface
import pwi4_client

import time
from datetime import datetime

from pyvista import tv
import focus

from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from astropy.io import fits
from astropy.time import Time

from astroquery.vo_conesearch import ConeSearch, conesearch
import astroquery.utils
astroquery.utils.suppress_vo_warnings()

import status
import multiprocessing as mp

dataroot='/data/1m/'

# discovery seems to fail on 10.75.0.0, so hardcode servers
#svrs=discovery.search_ipv4(timeout=30,numquery=3)
global D, S, T, F, Filt, C
def ascom_init() :
    svrs=['10.75.0.21:32227','10.75.0.22:11111']
    print("Alpaca devices: ")
    for svr in svrs:
        print(f"  At {svr}")
        print (f"    V{management.apiversions(svr)} server")
        print (f"    {management.description(svr)['ServerName']}")
        devs = management.configureddevices(svr)
        for dev in devs:
            print(f"      {dev['DeviceType']}[{dev['DeviceNumber']}]: {dev['DeviceName']}")

    # open Alpaca devices
    print('Opening Alpaca devices...')
    global D, S, T, F, Filt, C
    D=Dome(svrs[0],0)
    S=SafetyMonitor(svrs[0],0)
    T=Telescope(svrs[1],0)
    F=Focuser(svrs[1],0)
    Filt=FilterWheel(svrs[1],0)
    C=Camera(svrs[1],0)
    print()
    print("All ASCOM commands available through devices: ")
    print('    T : telescope commands')
    print('    C : camera commands')
    print('    F : focusser command')
    print('    Filt : filter wheel commands')
    print('    D : dome commands')
 

def pwi_init() :
    global pwi
    print('Starting PWI4 client...')
    pwi=pwi4_client.PWI4(host='pwi1m')
    for axis in [0,1] : pwi.mount_enable(axis)

# Camera commands
def qck(exptime,filt='current') :
    """ (shorthand) Take exposure without saving to disk, current filter by default
    """
    expose(exptime,filt)

def expose(exptime=1.0,filt='current',bin=3,box=None,light=True,display=None,name=None) :
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
        C.NumX = box.ncol()
        C.NumY = box.nrow()
    else :
        C.StartX = 0
        C.StartY = 0
        C.NumX = C.CameraXSize//bin
        C.NumY = C.CameraYSize//bin
    t = Time.now()
    C.StartExposure(exptime,light)
    while not C.ImageReady :
        time.sleep(0.5)
    data = np.array(C.ImageArray).T
    if display is not None :
        display.tv(data)
        display.fig.canvas.flush_events()

    hdu=fits.PrimaryHDU(data)
    hdu.header['DATE-OBS'] = t.fits
    hdu.header['JD'] = t.jd
    hdu.header['MJD'] = t.mjd
    hdu.header['EXPTIME'] = exptime 
    hdu.header['FILTER'] = filt 
    hdu.header['FOCUS'] = F.Position
    stat = pwi.status()
    hdu.header['RA'] = stat.mount.ra_j2000_hours
    hdu.header['DEC'] = stat.mount.dec_j2000_degs
    hdu.header['AZ'] = stat.mount.azimuth_degs
    hdu.header['ALT'] = stat.mount.altitude_degs
    hdu.header['ROT'] = stat.rotator.mech_position_degs
    hdu.header['TELESCOP'] = 'APO SONG 1m'
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
    """ Change CCD temperature set point
    """
    C.SetCCDTemperature = temp

def focrun(cent,step,n,extime=1.0,filt='V',bin=3,box=None,display=None) :
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
    display : pyvista tv, default=None
           If given, display each image as it is being taken

    Returns
    -------
           List of file names taken for focus run
    """

    files=[]
    images=[]
    for foc in np.arange(int(cent)-n//2*int(step),int(cent)+n//2*int(step)+1,int(step)) :
        F.Move(int(foc))
        while F.IsMoving :
            time.sleep(1)
        print('position: ',F.Position)
        hdu,name = expose(exptime,filt,box=box,bin=bin,display=display,name='focus_{:d}'.format(foc))
        files.append(name)
        images.append(hdu)
    pixscale=206265/6000*C.PixelSizeX*1.e-3*bin
    focus.focus(files,pixscale=pixscale,display=display)
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
    if isinstance(ra,float) and isinstance(dec,float) :
        pwi.mount_goto_ra_dec_j2000(ra,dec)
    elif isinstance(ra,str) and isinstance(dec,str) :
        coords=SkyCoord("{:s} {:s}".format(ra,dec),unit=(u.hourangle,u.deg))
        pwi.mount_goto_ra_dec_j2000(coords.ra.value/15.,coords.dec.value)

    #tra, tdec = j2000totopocentric(ra,dec) 
    #T.SlewToCoordinatesAsync(tra,tdec)
    while T.Slewing :
        time.sleep(1)
    T.Tracking = True
    D.Slaved = True
    #domesync()

def usno(ra=None,dec=None,rad=1*u.degree,rmin=0,rmax=15,bmin=0,bmax=15,goto=True) :
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
    """
    if ra is None :
        stat = pwi.status()
        ra = stat.mount.ra_j2000_hours
        dec = stat.mount.dec_j2000_degs
        coords=SkyCoord("{:f} {:f}".format(ra,dec),unit=(u.hourangle,u.deg))
    else :
        coords=SkyCoord("{:s} {:s}".format(ra,dec),unit=(u.hourangle,u.deg))
    cat= 'The USNO-A2.0 Catalogue 1'
    result = conesearch.conesearch(coords,rad,catalog_db=cat,verbose=False)
    gd=np.where((result['Bmag']<bmax) & (result['Bmag']>bmin) &
                (result['Rmag']<rmax) & (result['Rmag']>rmin))[0]
    if len(gd) <= 0 :
        print('no USNO A2.0 stars found within specified rad and mag range')
        return
    best = result['_r'][gd].argmin()
    print(result[gd[best]])
    if goto: slew(result['RAJ2000'][gd[best]], result['DEJ2000'][gd[best]])

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

    Paramters
    ---------
    tracking : bool
               if True, turn tracking on, if False, turn tracking off
    """
    T.Tracking= tracking

def park() :
    """ Park telescope
    """
    T.Park()
    while T.Slewing :
        time.sleep(1)

    D.Slaved = False
    D.Park()

def open() :
    """ Open dome
    """
    D.OpenShutter()

def close() :
    """ Close dome
    """
    D.CloseShutter()

def domesync(update=60) :

    def sync(update=60) :
        while D.Slaved :
            D.SlewToAzimuth(T.Azimuth)
            time.sleep(update)

    t=mp.Process(target=sync)
    t.start()

def start_status() :   
    """ Start status window thread
    """
    global proc
    proc = mp.Process(target=status.status,kwargs={'pwi' : pwi})
    proc.start()

def stop_status() :
    """ Stop status window thread
    """
    proc.terminate()

def commands() :
    print()
    print("Dome commands")
    print("  open(): open dome")
    print("  close(): close dome")
    print()
    print("Telescope commands")
    print("  slew(ra,dec): slew to coordinates")
    print("  offset(ra,dec): offset telescope")
    print("  usno([ra,dec]) : find/slew to USNO A2.0 star")
    print("  park(): park telescope")
    print("  tracking(True|False): turn tracking on/off")
    print()
    print("Camera commands")
    print("  expose(exptime,filt,**kwargs: take an exposure")
    print("  focrun(cent,step,nsteps,exptime,filt,**kwargs): take series of exposures at different focus positions")
    print("  settemp(temp): set camera temperature set point")
    print()
    print("Status commands")
    print("  start_status(): start status window")
    print("  stop_status(): stop status window ")
    print()
    print("Use help(command) for more details")

def start() :
    global D, S, T, F, Filt, C
    ascom_init()
    pwi_init()
    start_status()
    commands()

ascom_init()
pwi_init()

if __name__ == '__main__' :
    start_status()
    commands()
