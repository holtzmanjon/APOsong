import matplotlib
import glob
import os
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()
from pyvista import imred, tv, centroid, image, stars
import time
import pdb
from astropy.io import fits
import numpy as np
from holtztools import html
import re
import select
import sys
import socket

import aposong
import eshel
import influx

display=None
run_guide = True

# set up logger
import logging
import yaml
import logging.config
try: 
    with open('logging.yml', 'rt') as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)
    logger=logging.getLogger('aposong')
except FileNotFoundError :
    #trap for readthedocs
    print('logging.yml not found')


class Guider :

    def __init__(self,x0=774,y0=462,exptime=5,filt=None, bin=1,rad=25,skyrad=[35,50],mask=None,maskrad=6,
          display=None,nintegral=10, prop=0.7,ki=0.2,nint=10,settle=1,pixscale=1.) :
        self.x0 = x0
        self.y0 = y0
        self.box=image.BOX(cr=int(y0),cc=int(x0),n=int(7*rad))
        self.target = [x0-self.box.xmin,y0-self.box.ymin]
        self.guess = [x0-self.box.xmin,y0-self.box.ymin]
        self.exptime = float(exptime)
        self.filt = filt
        self.bin = bin
        self.rad = rad
        self.skyrad = skyrad
        self.mask = mask
        self.maskrad = maskrad
        self.disp = display
        self.prop = prop
        self.ki = ki
        self.navg = np.min([5,np.max([1,int(5/exptime)])])
        self.nseq = 1
        self.xtot = 0.
        self.ytot = 0.
        self.nintegral = nintegral
        self.integral = np.zeros([nintegral,2])
        self.nint = 0
        self.settle = settle
        self.pixscale = pixscale
        self.mask = mask
        if mask is not None : mask=image.window(mask,box=box)

    def get_hole_position(self) :
        """ Refine position of target hole in pickoff mirror, return position and mask
        """
        # illuminate the aperture
        eshel.lamps(mirror=True,quartz=True)
        time.sleep(3)
        # expose
        im=aposong.expose(0.01,bin=1,filt=None,cam=0,max=50000,display=self.disp,name='guide/hole')
        time.sleep(3)
        eshel.lamps(close=False)
        time.sleep(3)

        # find minimum pixel (maximum of negative image)
        box=image.BOX(cr=int(self.y0),cc=int(self.x0),n=25)
        stats=image.abx(-im.hdu.data,box)
        x0 = stats['peakx']
        y0 = stats['peaky']
   
        # get center from image and return
        center=centroid.rasym_centroid(im.hdu.data,x0,y0,rad=12,plot=self.disp,mask=self.mask)
        logger.info('hole center : {:.2f} {:.2f}'.format(center.x,center.y))

        if self.disp is not None :
            self.disp.tvclear()
            self.disp.tvcirc(center.x,center.y,rad=self.maskrad)

        self.x0 = center.x
        self.y0 = center.y

        self.mkmask(im.hdu.data,center.x, center.y)

    def mkmask(self,data,x0,y0) :
        """ Create a mask around desired position
        """
        mask=np.zeros_like(data)
        yg,xg=np.mgrid[0:data.shape[0],0:data.shape[1]]
        r2=(xg-x0)**2+(yg-y0)**2
        bd=np.where(r2<self.maskrad**2)
        mask[bd]=1
        self.mask=image.window(mask,box=self.box)


    def acquire(self,thresh=100,fwhm=1.5,data=None,inter=False) :
        """ Find brightest object, optimize exposure time, move object to desired position
        """
        if data is None :
            exp=aposong.expose(self.exptime,filt=self.filt,bin=self.bin,display=self.disp,name='guide/acquire')
            hdu=exp.hdu
        else :
            hdu = data

        if inter :
            self.disp.tv(hdu,max=vmax)
            print('mark star on display: ')
            k,x,y=self.disp.tvmark()
        else :
            # find brightest star automatically
            # adjust exposure time to get good counts
            niter=0
            while True :
                mad=np.nanmedian(np.abs(hdu.data-np.nanmedian(hdu.data)))
                objs=stars.find(hdu.data,thresh=thresh*mad,fwhm=fwhm/self.pixscale,brightest=1)
                if objs is None : 
                    raise RuntimeError('no objects found')
                    return
                peak=objs[0]['peak']
                # if peak between 20000 and 60000 DN, or exptime>=10s, accept exposure time
                if (peak < 60000 and peak>20000 ) or niter>10 or self.exptime>9.99: break
                if peak>60000 : self.exptime*=0.1
                else: self.exptime = np.max([0.01,np.min([self.exptime*50000/peak,10])])
                logger.info('new exptime: {:.2f}'.format(self.exptime))
                exp=aposong.expose(self.exptime,filt=self.filt,bin=self.bin,display=self.disp)
                hdu=exp.hdu
                niter+=1

            x=objs[0]['x']
            y=objs[0]['y']

        self.navg = np.min([5,np.max([1,int(5/self.exptime)])])

        # take full frame acquisition image to save and offset to brightest object
        exp=aposong.expose(self.exptime,display=self.disp,filt=self.filt,bin=self.bin,name='guide/acquire')
        mad=np.nanmedian(np.abs(hdu.data-np.nanmedian(hdu.data)))
        objs=stars.find(hdu.data,thresh=thresh*mad,fwhm=fwhm/self.pixscale,brightest=1)
        x=objs[0]['x']
        y=objs[0]['y']
        logging.info('offsetting {:.1f} {:.1f}'.format(x-self.x0,y-self.y0))
        aposong.offsetxy((x-self.x0),(y-self.y0),scale=self.pixscale)
        time.sleep(self.settle)
        return self.get_offset(name='guide/acquire')

    def get_offset(self,data=None,name=None) :
        """ Take exposure and get offset from star to desired center
        """
        if data is None :
            exp=aposong.expose(self.exptime,display=self.disp,filt=self.filt,bin=self.bin,box=self.box,name=name,max=5000)
            data = exp.hdu.data
        else :
            self.disp.tv(data,max=5000)
        x,y=self.guess
        center=centroid.rasym_centroid(data,x,y,self.rad,mask=self.mask,skyrad=self.skyrad,plot=self.disp)
        bad=False
        if center.x<0 or center.tot<10000:
            try :
                # if rasync_centroid fails, try marginal_gfit
                logger.debug('marginal: {:.2f} {:.2f} {:.2f}'.format(x,y,self.rad))
                center=centroid.marginal_gfit(data,x,y,self.rad)
                logger.debug('marginal: {:.2f} {:.2f} {:.2f}'.format(center.x,center.y,center.tot))
                if center.tot < 10000 : 
                    bad=True
            except    :
                # if marginal_gfit fails, use peak
                logger.debug('peak: {:.2f} {:.2f} {:.2f}'.format(x,y,self.rad))
                try: 
                    center=centroid.peak(data,x,y,self.rad)
                except:
                    logger.debug('error in peak')
                    bad=True
                if center.tot < 1000 :
                    logger.debug('peak<1000, no offset')
                    bad=True
        else :
            bad=True
        if self.disp is not None and not bad :
            self.disp.tvcirc(center.x,center.y,rad=2)
            self.disp.tvcirc(self.x0,self.y0,rad=1,color='g')
            self.guess = [center.x,center.y]

        self.center = center
        x0=self.x0-self.box.xmin
        y0=self.y0-self.box.ymin
        logger.debug('  instantaneous offset: {:d} {:.1f} {:.1f} {:.1f}'.format(self.nseq,self.center.x-x0,self.center.y-y0,self.center.tot))
        return self.center[0], self.center[1]

    def guide(self) :
        """ Take image and accumulate offset, make correction if accumulated
        """
        x,y = self.get_offset(name='guide/guide')
        if self.center.x>0 :
            self.xtot+=self.center.x
            self.ytot+=self.center.y
            x0=self.x0-self.box.xmin
            y0=self.y0-self.box.ymin
            if self.nseq < self.navg :
                logger.debug('  instantaneous offset: {:d} {:.1f} {:.1f} {:.1f}'.format(self.nseq,self.center.x-x0,self.center.y-y0,self.center.tot))
                logger.debug('  accumulated offset:   {:.1f} {:.1f}'.format(self.xtot/self.nseq-x0,self.ytot/self.nseq-y0))
                dx=self.center.x-x0
                dy=self.center.y-y0
                self.nseq+=1
                self.ingest_offset(dx,dy)
            else :
                logger.debug('  AVERAGE OFFSET: {:.1f} {:.1f} {:.1f}'.format(self.xtot/self.nseq-x0,self.ytot/self.nseq-y0,self.center.tot))
                # average offsets over navg points
                x=self.xtot/self.navg
                dx=x-x0
                y=self.ytot/self.navg
                dy=y-y0
                # integral offsets over maxint points
                self.integral[self.nint%self.nintegral]=[dx,dy]
                self.nint+=1
                # if offset>2 arcsec, apply full offset less one pixel
                # elif >0.1 arcsec, apply prop*offset
                # else no offset 
                if abs(dx)*self.pixscale > 2 : xoff=dx-dx/abs(dx)
                elif abs(dx)*self.pixscale > 0.1 : xoff=self.prop*dx + self.ki*self.integral[:,0].mean()
                else : xoff=0
                logger.debug('   APPLIED X: {:.1f} {:.1f} {:.1f}'.format(xoff,self.prop*dx,self.ki*self.integral[:,0].mean()))
                
                if abs(dy)*self.pixscale > 2 : yoff=dy-dy/abs(dy)
                elif abs(dy)*self.pixscale > 0.1 : yoff=self.prop*dy + self.ki*self.integral[:,1].mean()
                else : yoff=0.
                logger.debug('   APPLIED Y: {:.1f} {:.1f} {:.1f}'.format(yoff,self.prop*dy,self.ki*self.integral[:,1].mean()))
                aposong.offsetxy(xoff,yoff,scale=self.pixscale)
                time.sleep(self.settle)
                self.ingest_correction(x,y,dx,dy,xoff,yoff)
                # reset counter and accumulators
                self.nseq=1
                self.xtot=0
                self.ytot=0
                # set starting guess for next centroid to be guided position
                self.guess = [self.x0-self.box.xmin,self.y0-self.box.ymin]

    def ingest_offset(self,dx,dy) :
        """ Ingest single image offset into InfluxDB
        """
        idict={}
        for k,v in zip(['dx','dy','nseq'],
                       [float(dx*self.pixscale),float(dy*self.pixscale),self.nseq]) :
            idict[k] = v
        influx.write(idict,bucket='guider',measurement='my_measurement')

    def ingest_correction(self,x,y,dx,dy,xoff,yoff) :
        """ Ingest correction offset into InfluxDB
        """
        idict={}
        for k,v in zip(['x0','y0','x','y','dx','dy','xoff','yoff','nseq'],
                       [float(self.x0+self.box.xmin),float(self.y0+self.box.ymin),
                        float(x),float(y),
                        float(dx*self.pixscale), float(dy*self.pixscale),
                        float(xoff*self.pixscale), float(yoff*self.pixscale),
                        self.nseq]) :
            idict[k] = v
        influx.write(idict,bucket='guider',measurement='my_measurement')
    
def doguide(x0,y0,rad=25,exptime=5,filt=None,bin=1,n=1,navg=1,mask=None,disp=None,name='guide/guide',
            prop=0.7, ki=0.2, maxint=10,rasym=True,settle=3,vmax=10000,box=None,date='UT240503',weight=False) :
    """ OLD Actual guide loop
    """
    # accumulators for navg offsets 
    nseq=1 
    xtot=0
    ytot=0
    x,y = x0,y0
    red=imred.Reducer(dir='/data/1m/{:s}/guide/'.format(date))

    pixscale=aposong.pixscale()

    write_api = influx.setup_influx()
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
                    logger.debug('marginal: {:.2f} {:.2f} {:.2f} {:d}'.format(x,y,rad,n))
                    center=centroid.marginal_gfit(hdu.data,x,y,rad)
                    logger.debug('marginal: {:.2f} {:.2f} {:.2f}'.format(center.x,center.y,center.tot))
                    if center.tot < 10000 : continue
                except:
                    # if marginal_gfit fails, use peak
                    logger.debug('peak: {:.2f} {:.2f} {:.2f} {:d}'.format(x,y,rad,n))
                    try: 
                        center=centroid.peak(hdu.data,x,y,rad)
                    except:
                        logger.debug('error in peak')
                        continue
                    if center.tot < 1000 :
                        logger.debug('peak<1000, no offset')
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

def mkmovie(ut,root='/data/1m/') :
    """ Make guider movies for guider sequences on specifed UT
    """
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

def show(n,date='UT240503',sleep=0.5,max=30000,pause=False,weight=False,rad=25,maskrad=7,prop=0.7) :
    """ Display guider images from specified date
    """
    global display

    if display is None :
        display=tv.TV(figsize=(8,4))

    red=imred.Reducer(dir='/data/1m/{:s}/guide/'.format(date))
    red.formstr=['*.{:04d}.f*']
    pixscale=aposong.pixscale()

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


def parse_string_to_kwargs(string):
  kwargs= {}
  for word in string.split() :
    pattern = r"(.*)=(.*)"  # Matches key=value pairs
    matches = re.findall(pattern, word)

    for key, value in matches:
        kwargs[key] = value

  return kwargs

def main() :

    #fp = os.open('pipe',os.O_RDONLY | os.O_NONBLOCK)
    print('opening socket...')
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.setblocking(False)  # Make the socket non-blocking
    server_socket.bind(('localhost', 5000))
    server_socket.listen()
    conn = None
    disp=tv.TV()

    inputs=[sys.stdin,server_socket]
    outputs=[]

    guiding = False
    while True :
      try:
        readable, writable, exceptional = select.select(inputs, outputs, inputs, 0.)
 
        for s in readable : 
            #if s is fp:
            #   txt = os.read(fp,80).decode().strip('\n')
            txt=''
            if s is sys.stdin :
                txt = sys.stdin.readline().strip('\n')
                print('stdin: ',txt)
            if s is server_socket :
                conn, addr = server_socket.accept()
                inputs.append(conn)
                print('adding client')
            if s is conn :
                txt = conn.recv(1024).decode()
                print('socket: ', txt)
            if txt == 'start' :
                print('starting...')
                kwargs=parse_string_to_kwargs(txt)
                g = Guider(pixscale=aposong.pixscale(),display=disp,**kwargs)
                g.get_hole_position()
                try : 
                    cent=g.acquire()
                    guiding = True
                    if s is conn : conn.send(f'acquired {cent[0]} {cent[1]}'.encode())
                    else : print(f'acquired {cent[0]} {cent[1]}')
                except :
                    if s is conn : conn.send(b'failed acquire')
                    else : print('failed acquire')
            elif txt == 'stop' or txt == 'pause' :
                guiding = False
            elif txt == 'resume ' :
                guiding = True
            else :
                print('unknown command: ', txt, len(txt))
 
        if guiding :
              g.guide()

        time.sleep(0.25)
      except KeyboardInterrupt :
        server_socket.close()
        break

#main()
