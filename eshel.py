# routines for use with Shelyak eShel spectograph and calibration channel
import aposong
import time

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
  
#def cals(display=None,flats=15,thar=15,flat_exptime=8,thar_exptime=60,cam=3,bin=2,root='') :
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

