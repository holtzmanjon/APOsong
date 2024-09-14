import matplotlib
import glob
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()
from pyvista import imred, tv, centroid
import time
import pdb
from astropy.io import fits
import numpy as np
from holtztools import html

display=None
def show(n,date='UT240503',sleep=0.5,max=30000,pause=False) :
    global display

    if display is None :
        display=tv.TV(figsize=(8,4))

    red=imred.Reducer(dir='/data/1m/{:s}/guide/'.format(date))
    while True :
      try :
        a=red.rd(n)
        print(a.header['EXPTIME'])
        display.tv(a,max=max)
        if pause and a.shape[0] > 200 :
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
try:
    with open('logging.yml', 'rt') as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)
except:
    print("can't open logging.yml")

logger=logging.getLogger(__name__)

import aposong

def doguide(x0,y0,rad=25,exptime=5,filt=None,bin=1,n=1,navg=1,mask=None,disp=None,name='guide/guide',
            prop=0.7, rasym=True,settle=3,vmax=10000,box=None,date='UT240503',weight=False) :
    """ Actual guide loop
    """
    # accumulators for navg offsets 
    nseq=1 
    xtot=0
    ytot=0
    x,y = x0,y0
    red=imred.Reducer(dir='/data/1m/{:s}/guide/'.format(date))

    pixscale=aposong.pixscale()

    while run_guide :
        logger.debug('guide start: {:.1f} {:.1f} {:.1f} {:.1f} {:.2f} {:d}'.format(x,y,prop,bin,exptime,navg))
        if disp is not None : disp.tvclear()
        if exptime > 0 :
            exp=aposong.expose(exptime,display=disp,bin=bin,filt=filt,max=vmax,box=box,
                               name=name,insert=False)
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
            logger.debug('asym: {:.2f} {:.2f} {:.2f}'.format( center.x,center.y,center.tot))
            if center.x<0 or center.tot < 10000:
                try :
                    # if rasync_centroid fails, try marginal_gfit
                    center=centroid.marginal_gfit(hdu.data,x,y,rad)
                    logger.debug('marginal: {:.2f} {:.2f} {:.2f}'.format(center.x,center.y,center.tot))
                    if center.tot < 10000 : continue
                except:
                    # if marginal_gfit fails, use peak
                    center=centroid.peak(hdu.data,x,y)
                    if center.tot < 1000 :
                        logger.info('peak<1000, no offset')
                        continue
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
                    dx=xtot/nseq-x0
                    xrem=(1-prop)*dx
                    if abs(xrem)*pixscale > 2 : xoff=dx
                    elif abs(xrem)*pixscale > 0.5 : xoff=prop*dx
                    else : xoff=0
                    dy=ytot/nseq-y0
                    yrem=(1-prop)*dy
                    if abs(yrem)*pixscale > 2 : yoff=dy
                    elif abs(yrem)*pixscale > 0.5 : yoff=prop*dy
                    else : yoff=0.
                    aposong.offsetxy(xoff,yoff,scale=pixscale)
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

def mkmovie(ut,root='/data/1m/') :
    matplotlib.use('Agg')
    dir=root+ut+'/guide'
    red=imred.Reducer(dir=dir)
    files=glob.glob(dir+'/acquire*.fits')
    ims=[]
    seqs=[]
    for file in files :
        im=int(file.split('.')[-2])
        ims.append(im)
    for i,im in enumerate(ims[1:]) :
        if ims[i]-ims[i-1] > 1 :
            seqs.append([ims[i-1]+1,ims[i]])

    print(ims)
    print(seqs)
    grid=[]
    row=[]
    for i,seq in enumerate(seqs):
        plt.close('all')
        t=tv.TV()
        out='{:s}/{:d}.gif'.format(dir,seq[0])
        #red.movie(range(seq[0],seq[1]),display=t,max=10000,out=out)
        if i>0 and i%5 == 0 :
            grid.append(row)
            row=[]
        row.append('guide/{:d}.gif'.format(seq[0]))
    html.htmltab(grid,file=root+ut+'/guide.html',size=200)
