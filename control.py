import matplotlib
matplotlib.use('TkAgg')
import numpy as np

from alpaca import discovery, management
from alpaca.telescope import *
from alpaca.dome import *
from alpaca.safetymonitor import *
from alpaca.focuser import *
from alpaca.filterwheel import *
from alpaca.camera import *

import time
from datetime import datetime

from pyvista import tv

from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u


#svrs=discovery.search_ipv4(timeout=30,numquery=3)
svrs=['10.75.0.21:32227','10.75.0.22:11111']
for svr in svrs:
    print(f"At {svr}")
    print (f"  V{management.apiversions(svr)} server")
    print (f"  {management.description(svr)['ServerName']}")
    devs = management.configureddevices(svr)
    for dev in devs:
        print(f"    {dev['DeviceType']}[{dev['DeviceNumber']}]: {dev['DeviceName']}")

D=Dome(svrs[0],0)
S=SafetyMonitor(svrs[0],0)
T=Telescope(svrs[1],0)
F=Focuser(svrs[1],0)
Filt=FilterWheel(svrs[1],0)
C=Camera(svrs[1],0)


def expose(exptime,filt,binx=3,biny=3,light=True,display=None,file=None) :
    Filt.Position=filt
    C.BinX=binx
    C.BinY=biny
    C.StartExposure(exptime,light)
    while not C.ImageReady :
        time.sleep(0.5)
    data = np.array(C.ImageArray).T
    if display is not None :
        display.tv(data)



def j2000totopocentric(ra,dec) :

    #apo=EarthLocation.of_site('APO')
    #d=datetime.now()
    #utc = tpm.gcal2j(d.year,d.month,d.day)
    #tt = tpm.utc2tdb(utc)
    #v6 = convert.cat2v6(ra*np.pi/180,dec*np.pi/180)
    #v6_app = convert.convertv6(s1=6,s2=16,lon=apo.lon.value,lat=apo.lat.value,alt=apo.height.value)


    from coordio import sky, site, time
    coords=sky.ICRS(np.array([[ra,dec]]))
    s=site.Site('APO')
    s.set_time(time.Time())    # defaults to now
    obs=sky.Observed(coords,site=s)
    return obs.ra, obs.dec


