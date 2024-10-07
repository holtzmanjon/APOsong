import matplotlib
import glob
import os
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()
from pyvista import imred, tv, centroid
import time
import pdb
from astropy.io import fits
import numpy as np
from holtztools import html

import influxdb_client
from influxdb_client.client.write_api import SYNCHRONOUS

bucket = "guider"
org = "NMSU"
def setup_influx() :
    url="http://localhost:8086"
    os.environ['INFLUX_TOKEN'] = 'v-RuHY6T1pyOIL1SU9lrWYKYEU_SDZ0VWkPHOIU9hMECF7axu2wiFzY1u8N7J6s9KCbOreQKI43mJUi9pj5BbA=='
    # Store the URL of your InfluxDB instance
    url="http://localhost:8086"
    token = os.environ['INFLUX_TOKEN']
    client = influxdb_client.InfluxDBClient(
       url=url,
       token=token,
       org=org
    )
    write_api = client.write_api(write_options=SYNCHRONOUS)
    return write_api


display=None
def show(n,date='UT240503',sleep=0.5,max=30000,pause=False,weight=False,rad=25,maskrad=7,prop=0.7) :
    global display

    if display is None :
        display=tv.TV(figsize=(8,4))

    red=imred.Reducer(dir='/data/1m/{:s}/guide/'.format(date))
    red.formstr=['*.{:04d}.f*']
    pixscale=aposong.pixscale()

    write_api = setup_influx()
    while True :
      try :
        hdu=red.rd(n)
        print(hdu.header['EXPTIME'])
        mask=np.zeros_like(hdu.data)
        yg,xg=np.mgrid[0:hdu.data.shape[0],0:hdu.data.shape[1]]
        x0=hdu.data.shape[1]//2
        y0=hdu.data.shape[0]//2
        r2=(xg-x0)**2+(yg-y0)**2
        bd=np.where(r2<maskrad**2)
        mask[bd]=1


        display.tv(hdu,max=max)
        center=centroid.rasym_centroid(hdu.data,x0,y0,rad,mask=mask,skyrad=[35,40],plot=display,weight=weight)
        dx=(center.x-x0)
        dy=(center.y-y0)
        print(dx,dy,center.tot)
        if pause and hdu.shape[0] > 200 :
            display.imexam()

        if abs(dx)*pixscale > 2 : xoff=dx-dx/abs(dx)
        elif abs(dx)*pixscale > 0.5 : xoff=prop*dx
        else : xoff=0
        if abs(dy)*pixscale > 2 : yoff=dy-dy/abs(dy)
        elif abs(dy)*pixscale > 0.5 : yoff=prop*dy
        else : yoff=0.


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
            prop=0.7, ki=0.2, maxint=10,rasym=True,settle=3,vmax=10000,box=None,date='UT240503',weight=False) :
    """ Actual guide loop
    """
    # accumulators for navg offsets 
    nseq=1 
    xtot=0
    ytot=0
    x,y = x0,y0
    red=imred.Reducer(dir='/data/1m/{:s}/guide/'.format(date))

    pixscale=aposong.pixscale()

    write_api = setup_influx()
    if box is not None :
        x0off = box.xmin
        y0off = box.ymin
    else : 
        x0off = 0
        y0off = 0
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
                    logger.info('marginal: {:.2f} {:.2f} {:.2f} {:d}'.format(x,y,rad,n))
                    center=centroid.marginal_gfit(hdu.data,x,y,rad)
                    logger.debug('marginal: {:.2f} {:.2f} {:.2f}'.format(center.x,center.y,center.tot))
                    if center.tot < 10000 : continue
                except:
                    # if marginal_gfit fails, use peak
                    logger.info('peak: {:.2f} {:.2f} {:.2f} {:d}'.format(x,y,rad,n))
                    try: 
                        center=centroid.peak(hdu.data,x,y,rad)
                    except:
                        logger.debug('error in peak')
                        continue
                    if center.tot < 1000 :
                        logger.info('peak<1000, no offset')
                        continue
        else :
            center==centroid.marginal_gfit(hdu.data,x,y,rad)
            tot=-9

        xint=np.zeros(maxint)
        yint=np.zeros(maxint)
        nint=0
        if center.x>0 and center.y>0 :
            if disp is not None: 
                disp.tvcirc(center.x,center.y,rad=2)
                disp.tvcirc(x0,y0,rad=1,color='g')
            xtot+=center.x
            ytot+=center.y
            if nseq < navg :
                logger.debug('  instantaneous offset: {:d} {:.1f} {:.1f} {:.1f}'.format(nseq,center.x-x0,center.y-y0,center.tot))
                logger.debug('  accumulated offset:   {:.1f} {:.1f}'.format(xtot/nseq-x0,ytot/nseq-y0))
                dx=center.x-x0
                dy=center.y-y0
                p = [influxdb_client.Point("my_measurement").tag("location", "APO").field("dx", float(dx*pixscale)),
                     influxdb_client.Point("my_measurement").tag("location", "APO").field("dy", float(dy*pixscale)),
                     influxdb_client.Point("my_measurement").tag("location", "APO").field("nseq", nseq)]
                write_api.write(bucket=bucket, org=org, record=p)
                nseq+=1
            else :
                logger.debug('  APPLIED OFFSET: {:.1f} {:.1f} {:.1f}'.format(xtot/nseq-x0,ytot/nseq-y0,center.tot))
                if exptime>0 : 
                    # average offsets over navg points
                    x=xtot/nseq
                    dx=xtot/nseq-x0
                    y=ytot/nseq
                    dy=ytot/nseq-y0
                    # integral offsets over maxint points
                    xint[nint//maxint]=dx
                    yint[nint//maxint]+=dy
                    nint+=1
                    # if offset>2 arcsec, apply full offset less one pixel
                    # elif >0.1 arcsec, apply prop*offset
                    # else no offset 
                    if abs(dx)*pixscale > 2 : xoff=dx-dx/abs(dx)
                    elif abs(dx)*pixscale > 0.1 : xoff=prop*dx + ki*xint.mean()
                    else : xoff=0
                  
                    if abs(dy)*pixscale > 2 : yoff=dy-dy/abs(dy) + ki*yint.mean()
                    elif abs(dy)*pixscale > 0.1 : yoff=prop*dy
                    else : yoff=0.
                    aposong.offsetxy(xoff,yoff,scale=pixscale)
                    p = [influxdb_client.Point("my_measurement").tag("location", "APO").field("x0", float(x0+x0off)),
                         influxdb_client.Point("my_measurement").tag("location", "APO").field("y0", float(y0+y0off)),
                         influxdb_client.Point("my_measurement").tag("location", "APO").field("x", float(x)),
                         influxdb_client.Point("my_measurement").tag("location", "APO").field("y", float(y)),
                         influxdb_client.Point("my_measurement").tag("location", "APO").field("dx", float(dx*pixscale)),
                         influxdb_client.Point("my_measurement").tag("location", "APO").field("dy", float(dy*pixscale)),
                         influxdb_client.Point("my_measurement").tag("location", "APO").field("xoff", float(xoff*pixscale)),
                         influxdb_client.Point("my_measurement").tag("location", "APO").field("yoff", float(yoff*pixscale)),
                         influxdb_client.Point("my_measurement").tag("location", "APO").field("nseq", nseq)]
                    write_api.write(bucket=bucket, org=org, record=p)
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
