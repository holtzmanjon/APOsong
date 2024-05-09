import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()
from pyvista import imred, tv, centroid
import time
import pdb
from astropy.io import fits
import numpy as np

display=None
def show(n,date='UT240503',sleep=0.5,max=30000) :
    global display

    if display is None :
        display=tv.TV(figsize=(8,4))

    red=imred.Reducer(dir='/data/1m/{:s}/guide/'.format(date))
    while True :
      print(n)
      try :
        a=red.rd(n)
        print(a.header['EXPTIME'])
        display.tv(a,max=max)
        if a.shape[0] > 200 :
            display.imexam()
        n+=1
      except: pass
      time.sleep(sleep)

run_guide = True
guide_process = None

import logging
import time
import yaml
import logging.config
with open('logging.yml', 'rt') as f:
    config = yaml.safe_load(f.read())
logging.config.dictConfig(config)
logger=logging.getLogger(__name__)

import aposong

def doguide(x0,y0,rad=25,exptime=5,filt=None,bin=1,n=1,navg=1,mask=None,disp=None,name='guide',
            prop=0.7, rasym=True,settle=3,vmax=10000,box=None,date='UT240503',weight=False) :
    """ Actual guide loop
    """
    # accumulators for navg offsets 
    nseq=1 
    xtot=0
    ytot=0
    x,y = x0,y0
    red=imred.Reducer(dir='/data/1m/{:s}/guide/'.format(date))

    while run_guide :
        logger.debug('guide start: {:.1f} {:.1f} {:.1f} {:.1f} {:.2f} {:d}'.format(x,y,prop,bin,exptime,navg))
        if disp is not None : disp.tvclear()
        if exptime > 0 :
            exp=aposong.expose(exptime,display=disp,bin=bin,filt=filt,max=vmax,box=box,name='guide/guide'.format(n))
            hdu=exp.hdu
        else :
            print(n)
            try: hdu=red.rd(n)
            except:
                time.sleep(0.5)
                continue
            if disp is not None : disp.tv(hdu,max=vmax)

        if rasym :
            center=centroid.rasym_centroid(hdu.data,x,y,rad,mask=mask,skyrad=[35,40],plot=disp,weight=weight)
            if center.x<0 :
                try :
                    # if rasync_centroid fails, try marginal_gfit
                    center=centroid.marginal_gfit(hdu.data,x,y,rad)
                    logger.debug('marginal: {:.2f} {:.2f} {:.2f}'.format(center.x,center.y,center.tot))
                except:
                    # if marginal_gfit fails, use peak
                    yp,xp = np.unravel_index(np.argmax(hdu.data),hdu.data.shape)
                    center.x = xp
                    center.y = yp
        else :
            center==centroid.marginal_gfit(hdu.data,x,y,rad)
            tot=-9

        if center.x>0 and center.y>0 :
            if disp is not None: 
                disp.tvcirc(center.x,center.y,rad=2)
            xtot+=center.x
            ytot+=center.y
            if nseq < navg :
                logger.debug('  instantaneous offset: {:d} {:.1f} {:.1f} {:.1f}'.format(nseq,center.x-x0,center.y-y0,center.tot))
                logger.debug('  accumulated offset:   {:.1f} {:.1f}'.format(xtot/nseq-x0,ytot/nseq-y0))
                nseq+=1
            else :
                logger.debug('  APPLIED OFFSET: {:.1f} {:.1f} {:.1f}'.format(xtot/nseq-x0,ytot/nseq-y0,center.tot))
                if exptime>0 : 
                    aposong.offsetxy(prop*(xtot/nseq-x0),prop*(ytot/nseq-y0),scale=aposong.pixscale())
                    time.sleep(settle)
                else :
                    time.sleep(0.5)
                nseq=1
                xtot=0
                ytot=0
            x=center.x
            y=center.y
        n+=1


def mkmask(data,x0,y0,maskrad=7) :
    """ Create a mask
    """
    mask=np.zeros_like(data)
    yg,xg=np.mgrid[0:data.shape[0],0:data.shape[1]]
    r2=(xg-x0)**2+(yg-y0)**2
    bd=np.where(r2<maskrad**2)
    mask[bd]=1
    return mask

