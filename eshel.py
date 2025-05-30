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
  
def cals(display=None,flats=15,thar=15,flat_exptime=8,thar_exptime=60,cam=3,bin=2) :
    """ Take series of eShel cals
    """
    if flats>0 :
        lamps(mirror=True,quartz=True,led=True)
        time.sleep(3)
        for i in range(flats) :
            aposong.expose(flat_exptime,filt=None,bin=bin,display=display,cam=cam,name='flat')
    if thar>0 :
        lamps(mirror=True,thar=True)
        time.sleep(3)
        for i in range(thar) :
            aposong.expose(thar_exptime,filt=None,bin=bin,display=display,cam=cam,name='thar')
    lamps()
    time.sleep(3)

