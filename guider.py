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
from astropy.time import Time
from astropy.table import Table
import numpy as np
from holtztools import html
import re
import select
import sys
import socket

import aposong
import database
import cal
import influx

run_guide = True

# set up logger
import logging
import yaml
import logging.config
import datetime
try: 
    with open('logging.yml', 'rt') as f:
        logconfig = yaml.safe_load(f.read())
    log_rollover_time = datetime.time(hour=14)
    logging.config.dictConfig(logconfig)
    logger=logging.getLogger('guider')
except FileNotFoundError :
    #trap for readthedocs
    print('logging.yml not found')

def parse_string_to_kwargs(string):
  kwargs= {}
  for word in string.split() :
    pattern = r"(.*)=(.*)"  # Matches key=value pairs
    matches = re.findall(pattern, word)

    for key, value in matches:
        kwargs[key] = value

  return kwargs

# run aposong.init() here to get aposong.config
aposong.init()

class Guider :

    def __init__(self,x0=aposong.config['hole_pos'][1],y0=aposong.config['hole_pos'][0],findhole=True,
                 exptime=5,expavg=1,filt=None, bin=1,rad=25,skyrad=[35,50],mask=None,maskrad=6,sat=65000,
                 display=True,nintegral=10, prop=0.5,ki=0.2,nint=10,settle=1,pixscale=1.,exptime_min=0.05) :
        x0=float(x0)
        y0=float(y0)
        findhole=int(findhole)
        self.x0 = x0
        self.y0 = y0
        self.box=image.BOX(cr=int(y0),cc=int(x0),n=int(7*rad))
        self.target = [x0-self.box.xmin,y0-self.box.ymin]
        self.guess = [x0-self.box.xmin,y0-self.box.ymin]
        self.findhole = findhole
        self.exptime = float(exptime)
        self.expavg = expavg
        self.filt = filt
        self.bin = bin
        self.rad = rad
        self.skyrad = skyrad
        self.mask = mask
        self.maskrad = maskrad
        self.sat = sat
        if display :
            self.disp=tv.TV(figsize=(8,4))
            self.disp.fig.canvas.manager.window.wm_geometry('-0-0')
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
        self.display_max = 30000
        self.exptime_min = exptime_min

    def get_hole_position(self) :
        """ Refine position of target hole in pickoff mirror, return position and mask
        """
        # illuminate the aperture
        foc=aposong.foc()
        aposong.calstage_find(display=self.disp)
        cal.lamps(mirror=True,quartz=True)
        time.sleep(3)
        # expose
        im=aposong.expose(0.01,bin=1,filt=None,cam=0,max=50000,display=self.disp,name='guide/hole')
        time.sleep(3)
        cal.lamps(close=False)
        time.sleep(3)
        aposong.calstage_out()
        aposong.foc(foc)

        # find minimum pixel (maximum of negative image)
        #box=image.BOX(cr=int(self.y0),cc=int(self.x0),n=15)
        #stats=image.abx(-im.hdu.data,box)
        #x0 = stats['peakx']
        #y0 = stats['peaky']
        x0 = self.x0
        y0 = self.y0
   
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
            exp=aposong.expose(self.exptime,avg=self.expavg,filt=self.filt,bin=self.bin,display=self.disp,name='guide/acquire',max=self.display_max)
            hdu=exp.hdu
        else :
            hdu = data

        if inter :
            self.disp.tv(hdu,max=self.display_max)
            print('mark star on display: ')
            k,x,y=self.disp.tvmark()
        else :
            # find brightest star automatically
            # adjust exposure time to get good counts
            niter=0
            while True :
                oldexptime=self.exptime
                mad=np.nanmedian(np.abs(hdu.data-np.nanmedian(hdu.data)))
                objs=stars.find(hdu.data,thresh=thresh*mad,fwhm=fwhm/self.pixscale,brightest=1)
                if objs is None : 
                    raise RuntimeError('no objects found')
                    return
                peak=objs[0]['peak']
                logger.info('niter: {:d}'.format(niter))
                # if peak between 20000 and 60000 DN, or exptime>=10s, accept exposure time
                if (peak < 60000 and peak>20000 ) or niter>10 or self.exptime>9.99: break
                if peak>60000 : self.exptime = np.max([self.exptime*0.1,self.exptime_min])
                else: self.exptime = np.max([self.exptime_min,np.min([self.exptime*50000/peak,10])])
                logger.info('new exptime: {:.3f}'.format(self.exptime))
                if self.exptime == oldexptime : break
                exp=aposong.expose(self.exptime,filt=self.filt,bin=self.bin,display=self.disp)
                hdu=exp.hdu
                niter+=1

            x=objs[0]['x']
            y=objs[0]['y']

        if self.exptime < 1 : self.expavg = round(1/self.exptime)
        self.navg = np.min([3,np.max([1,int(5/self.exptime)])])

        # take full frame acquisition image to save and offset to brightest object
        exp=aposong.expose(self.exptime,display=self.disp,filt=self.filt,bin=self.bin,name='guide/acquire')
        mad=np.nanmedian(np.abs(hdu.data-np.nanmedian(hdu.data)))
        objs=stars.find(hdu.data,thresh=thresh*mad,fwhm=fwhm/self.pixscale,brightest=1)
        x=objs[0]['x']
        y=objs[0]['y']
        logger.info('offsetting {:.1f} {:.1f}'.format(x-self.x0,y-self.y0))
        aposong.offsetxy((x-self.x0),(y-self.y0),scale=self.pixscale)
        time.sleep(self.settle)
        return self.get_offset(name='guide/acquire')

    def get_offset(self,data=None,name=None) :
        """ Take exposure and get offset from star to desired center
        """
        if data is None :
            exp=aposong.expose(self.exptime,avg=self.expavg,display=self.disp,filt=self.filt,bin=self.bin,
                               box=self.box,name=name,max=self.display_max)
            data = exp.hdu.data
        else :
            self.disp.tv(data,max=self.display_max)
        x,y=self.guess
        center=centroid.rasym_centroid(data,x,y,self.rad,mask=self.mask,skyrad=self.skyrad,sat=self.sat)
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

        if self.disp is not None and not bad :
            self.disp.tvclear()
            self.disp.tvcirc(center.x,center.y,rad=self.rad)
            self.disp.tvcirc(center.x,center.y,rad=2,color='r')
            self.disp.tvcirc(self.x0-self.box.xmin,self.y0-self.box.ymin,rad=self.maskrad,color='g')
            self.guess = [center.x,center.y]

        self.center = center
        x0=self.x0-self.box.xmin
        y0=self.y0-self.box.ymin
        logger.debug('  instantaneous offset: {:d} {:.1f} {:.1f} {:.1f}'.format(self.nseq,self.center.x-x0,self.center.y-y0,self.center.tot))
        return self.center[0], self.center[1]

    def guide(self) :
        """ Take image and accumulate offset, make correction if accumulated
        """
        try : x,y = self.get_offset(name='guide/guide')
        except :
            logger.error('  ERROR in get_offset')
            self.center.x = -1
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
                try : self.ingest_offset(dx,dy)
                except : print('error ingesting offset')
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
                try: 
                    aposong.offsetxy(xoff,yoff,scale=self.pixscale)
                except: 
                    logger.error('Error in offsetxy')
                time.sleep(self.settle)
                try: self.ingest_correction(x,y,dx,dy,xoff,yoff)
                except : print('error ingesting correction')
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

    def postgres_write(self,guiding,acquired) :
        """ Write guide status to database
        """
        tab = Table()
        tab['guiders_id'] = [1]
        tab['dateobs'] = [Time.now().fits]
        tab['acquired'] = '1' if acquired else '0'
        tab['guiding'] = '1' if guiding else '0'
        tab['guide_target_x'] = self.x0
        tab['guide_target_y'] = self.y0
        tab['exp_time'] = self.exptime
        tab['exp_avg'] = self.expavg
        tab['ins_at'] = [Time.now().fits]

        d=database.DBSession(host='song1m_db.apo.nmsu.edu',database='db_apo',user='song')
        d.ingest('public.guiders',tab,onconflict='update',constraintname='guiders_id')
        d.close()
        return tab
    
def loop() :
    """ Main guider loop, accepting commands from command line or socket
    """

    print('opening socket...')
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.setblocking(False)  # Make the socket non-blocking
    server_socket.bind(('localhost', 5000))
    server_socket.listen()
    conn = None

    inputs=[sys.stdin,server_socket]
    outputs=[]

    guiding = False
    acquired = False
    t0=Time.now()
    while True :
      try:
        readable, writable, exceptional = select.select(inputs, outputs, inputs, 0.25)

        for s in readable : 
            txt=''
            if s is sys.stdin :
                txt = sys.stdin.readline().strip('\n')
                logger.info('stdin: {:s}'.format(txt))
            if s is server_socket :
                conn, addr = server_socket.accept()
                inputs.append(conn)
                logger.info('adding client')
                continue
            if s is conn :
                txt = conn.recv(1024).decode()
                logger.info('socket: {:s}'.format(txt))
            try : cmd = txt.split()[0]
            except : cmd=''
            if cmd == 'start' :
                kwargs=parse_string_to_kwargs(txt)
                g = Guider(pixscale=aposong.pixscale(),**kwargs)
                g.postgres_write(guiding,acquired)
                if g.findhole : g.get_hole_position()
                if g.x0 < 0 :
                    if s is conn : conn.send(b'failed acquire')
                    else : logger.info('failed acquire')
                    try : plt.close(g.disp.fig)
                    except: pass
                    del g
                    continue
                try : 
                    cent=g.acquire()
                    acquired = True
                    guiding = True
                    if s is conn : conn.send(f'acquired {cent[0]} {cent[1]}, exptime: {g.exptime}, expavg: {g.expavg}, navg: {g.navg}'.encode())
                    else : logger.info(f'acquired {cent[0]} {cent[1]}, exptime: {g.exptime}, expavg: {g.expavg}, navg: {g.navg}')
                except :
                    if s is conn : conn.send(b'failed acquire')
                    else : logger.info('failed acquire')
                try: g.postgres_write(guiding,acquired)
                except: logger.info('failed postres_write')

            elif cmd == 'stop' or cmd == 'pause' :
                guiding = False
                if cmd == 'stop' : 
                    acquired = False
                    if s is conn : conn.send(b'guiding stopped')
                    else : logger.info('guiding stopped')
                else :
                    if s is conn : conn.send(b'guiding paused')
                    else : logger.info('guiding paused')
                try :
                    try : g.postgres_write(guiding,acquired)
                    except: logger.info('failed postres_write')
                    if cmd == 'stop' :
                        plt.close(g.disp.fig)
                        del g
                except : pass
            elif cmd == 'resume' :
                if acquired : 
                    guiding = True
                    if s is conn : conn.send(b'guiding resumed')
                    else : logger.info('guiding resumed')
                else :
                    if s is conn : conn.send(b'not resumed, no guide star acquired')
                    else : logger.info('not resumed, no guide star acquired')
                try: g.postgres_write(guiding,acquired)
                except: logger.info('failed postres_write')
            elif cmd == 'status' :
                if s is conn : conn.send(b'guider alive')
                else : logger.info('guider alive')
            elif len(txt) > 0 :
                try: g.postgres_write(guiding,acquired)
                except: logger.info('failed postres_write')
                if s is conn : 
                    conn.send('unknown command: {:s} {:d}'.format(txt, len(txt)).encode())
                else :
                    logger.info('unknown command: {:s} {:d}'.format(txt, len(txt)))

            if conn is not None :
                inputs.remove(conn)
                conn.close()
                conn = None

        if guiding :
              g.guide()

      except KeyboardInterrupt :
        server_socket.close()
        break

if __name__ == "__main__" :
    loop()
