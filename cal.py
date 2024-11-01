import time
import aposong

def dark(exptime,ccdtemp=None,cam=0,n=5,bin=1,filt=None,name=None) :
    """ Take one series of darks
    """

    if ccdtemp is not None :
        aposong.settemp(ccdtemp,cam=cam)
        time.sleep(300)
    for i in range(n) :
        aposong.expose(exptime,cam=cam,light=False,bin=bin,filt=filt,name=name)

def darks(temps=[-15,-10,-5],bins=[1],exps=[0,90,120,900],cam=0,n=11,filt=None) :
    """ Take multiple series of darks at different temps, bins, exptimes
    """
    for temp in temps :
        aposong.settemp(temp,cam=cam)
        time.sleep(300)
        for bin in bins :
            for exp in exps :
               name='dark_cam{:d}_temp{:d}_bin{:d}_t{:d}'.format(cam,temp,bin,exp)
               dark(exp,cam=cam,n=n,filt=filt,bin=bin,name=name)

