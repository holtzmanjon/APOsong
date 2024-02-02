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
from alpaca import discovery, management
from alpaca.telescope import *
from alpaca.dome import *
from alpaca.safetymonitor import *
from alpaca.focuser import *
from alpaca.filterwheel import *
from alpaca.camera import *

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

import status
import multiprocessing as mp

dataroot='/data/1m/'

# discovery seems to fail on 10.75.0.0, so hardcode servers
#svrs=discovery.search_ipv4(timeout=30,numquery=3)
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
    print('    FILT : filter wheel commands')
    print('    D : dome commands')

def pwi_init() :
    global pwi
    print('Starting PWI4 client...')
    pwi=pwi4_client.PWI4(host='pwi1m')
    for axis in [0,1] : pwi.mount_enable(axis)

# Camera commands
def expose(exptime,filt,bin=3,box=None,light=True,display=None,name=None) :
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
    pos = np.where(np.array(Filt.Names) == filt)
    if len(pos) == 0 :
        print('no such filter')
        print('available filters: ', Filt.Names)
        return

    Filt.Position=pos[0]
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

def focrun(cent,step,n,exptime,filt,bin=3,box=None,display=None) :
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
    for foc in np.arange(int(cent)-n//2*int(step),int(cent)+n//2*int(step)+1,int(step)) :
        F.Move(int(foc))
        print('position: ',F.Position)
        hdu,name = expose(exptime,filt,box=box,bin=bin,display=display,name='focus_{:d}'.format(foc))
        files.append(name)
    pixscale=206265/6000*C.PixelSizeX*1.e-3*bin
    focus.focus(files,pixscale=pixscale,display=display)
    return files

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

ascom_init()
pwi_init()
start_status()
commands()
