from astropy.coordinates import SkyCoord,EarthLocation, Angle
from astropy.time import Time, TimeDelta
from astropy.table import Table, hstack, vstack
from astropy.io import fits
import astropy.units as u
from astroplan import Observer, time_grid_from_range
import matplotlib
import matplotlib.pyplot as plt
from pyvista import imred,tv,spectra
from holtztools import html, plots

import copy
import glob
import io
import os
import pdb
import numpy as np
import time
import threading
import subprocess
import requests
import datetime

import logging
import yaml
import logging.config
 
import focus as dofocus
import reduce
import mail

try:
    with open('logging.yml', 'rt') as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)
    for handler in logging.root.manager.loggerDict['aposong'].handlers :
        if handler.name == 'daily' :
            handler.atTime = datetime.time(hour=14)
except:
    print("can't open logging.yml")

logger=logging.getLogger('aposong')

import aposong 
import cal
import database

class Target() :
    def __init__(self,name,ra,dec,mag=99.999,epoch=2000.) :
        self.name=name
        self.ra=ra
        self.dec=dec
        self.epoch=epoch
        self.mag=mag

    def table(self) :
        tab =Table()
        tab['targname'] = [self.name]
        tab['ra'] = [self.ra]
        tab['dec'] = [self.dec]
        tab['epoch'] = [self.epoch]
        tab['mag'] = [self.mag]
        return tab

    def acquire(self,display=None) :
        aposong.iodine_out()
        aposong.guide('stop')
        aposong.slew(self.ra,self.dec)

        ret = aposong.guide('start')
        if ret : time.sleep(30)
        return ret

class Schedule() :
    def __init__(self,name,min_airmass=1.005,max_airmass=1.8,nvisits=1,dt_visit=1.,nsequence=1) :
        self.name=name
        self.min_airmass=min_airmass
        self.max_airmass=max_airmass
        self.nvisits=nvisits
        self.dt_visit=dt_visit
        self.nsequence=nsequence

    def table(self) :
        tab =Table()
        tab['schedulename' ] = [self.name]
        tab['min_airmass'] = [self.min_airmass]
        tab['max_airmass'] = [self.max_airmass]
        tab['nvisits'] = [self.nvisits]
        tab['dt_visit'] = [self.dt_visit]
        tab['nsequence'] = [self.nsequence]
        return tab

class Sequence() :
    def __init__(self,name,filt=['U','B','V','R','I'],n_exp=[1,1,1,1,1],t_exp=[1,1,1,1,1],camera=[0,0,0,0,0],bin=[1,1,1,1,1]) :
        self.name=name
        self.n_exp=n_exp
        self.t_exp=t_exp
        self.filt=filt
        self.camera=camera
        self.bin=bin

    def length(self,acquisition=60,filtmove=5,read=15) :
        tot = acquisition
        for filt,nexp,texp in zip(self.filt,self.n_exp,self.t_exp) :
            tot+= filtmove+nexp*(texp+read)
        return tot

    def observe(self,name,display=None,fact=1,nfact=1,header=None,req_no=-1) :
        names = []
        # start back at foc0
        aposong.foc(foc0)
        foc1 = copy.copy(foc0)
        no_exp = 0
        for filt,nexp,texp,cam,bin in zip(self.filt,self.n_exp,self.t_exp,self.camera,self.bin) :
            for iexp in range(nexp*nfact) :
                # focus adjustment for altitude
                alt = np.max([np.min([75,aposong.T.Altitude]),35])
                # make relative focus adjustment to last position, in case we are using iodine with focus offset
                adjustedfoc = int(foc0 + 200*(70-alt)/30)
                logger.info('Focus adjustment: alt {:.0f} foc0: {:.0f} adjusted_foc: {:.0f}, foc1: {:.0f}'.format(alt,foc0,adjustedfoc,foc1))
                aposong.foc(adjustedfoc-foc1,relative=True)
                foc1=copy.copy(adjustedfoc)
                if filt == 'thar' :
                    logger.info('Expose camera: {:d} exptime: {:.2f} filt: {:s} bin: {:d}, '.format(cam,texp,filt,bin))
                    aposong.guide('pause')
                    time.sleep(10)
                    calnames=cal.cals(thar=nexp,flats=0,thar_exptime=texp,display=display)
                    names.extend(calnames)
                    aposong.guide('resume')
                    time.sleep(30)
                else :
                    logger.info('Expose camera: {:d} exptime: {:.2f} filt: {:s} bin: {:d}, '.format(cam,texp*fact,filt,bin))
                    exp=aposong.expose(texp*fact,filt,name=name,display=display,cam=cam,bin=bin,targ=self.name,header=header,imagetyp='STAR')
                    names.append(exp.name)
                    no_exp+=1
                    load_song_status(req_no,'exec',no_exp=no_exp)
                if (Time.now()-nautical_morn).to(u.hour) > 0*u.hour or not aposong.issafe(): 
                    logger.info('breaking sequence for twilight or not safe')
                    break
        aposong.iodine_out()
        return names

    def table(self) :
        tab =Table()
        tab['sequencename' ] = [self.name]
        tab['filter'] = [self.filt]
        tab['n_exp'] = [self.n_exp]
        tab['t_exp'] = [self.t_exp]
        tab['camera'] = [self.camera]
        tab['bin'] = [self.bin]
        return tab

class Observation() :
    def __init__(self,target,schedule,sequence,mjd,files) :
        self.target = target    
        self.schedule = schedule
        self.sequence = sequence
        self.mjd = mjd
        self.files = files

def hamax(dec,airmax,lat) :
    """ Determine maximum hour angle given maximum airmass, declination, and latitude
    """
    hamax=(1/airmax - np.sin(lat)*np.sin(dec) ) / np.cos(lat)*np.cos(dec)
    hamax=np.arccos(hamax)
    return hamax

def secz(ha,dec,lat) :
    """ Determine airmass (secz) given hour angle, declination, and latitude
    """
    secz = (np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.cos(ha))**-1
    return secz

def getrequests() :
    d=database.DBSession()
    tab=d.query(sql='SELECT request_pk,priority, targ.*,sched.*,seq.* FROM robotic.request as req  \
         JOIN robotic.target as targ on req.targname = targ.targname \
         JOIN robotic.schedule as sched on req.schedulename = sched.schedulename \
         JOIN robotic.sequence as seq on req.sequencename = seq.sequencename' ,fmt='table')
    d.close()
    return tab

def getsong(t=None,site='APO',verbose=True,max_airmass=2,dt_focus=10) :
    """ Try to get observing request from SONG database
    """
    apo=EarthLocation.of_site(site)
    if t is None :
        t = Time(Time.now(),location=apo)
    d=database.DBSession(host='song1m_db.apo.nmsu.edu',database='db_song',user='song')
    song_requests=d.query(sql="SELECT * FROM public.obs_request_4 as req JOIN public.obs_request_status_4 as stat ON req.req_no = stat.req_no WHERE stat.status = 'wait'" ,
                         fmt='table')
    d.close()
    if len(song_requests) == 0 : return None, None

    song_requests.rename_column('req_prio','priority')
    song_requests.rename_column('right_ascension','ra')
    song_requests.rename_column('declination','dec')
    song_requests.rename_column('object_name','targname')
    priorities = sorted(set(song_requests['priority']),reverse=True)
    for priority in priorities :
      gd=np.where(song_requests['priority'] == priority) 
      for request in song_requests[gd] :
          ra=Angle(request['ra'], unit=u.hour).to_string(sep=':')
          dec=Angle(request['dec'], unit=u.degree).to_string(sep=':')
          c=SkyCoord("{:s} {:s}".format(ra,dec),unit=(u.hourangle,u.deg))
          ha = t.sidereal_time('mean') - c.ra
          ha.wrap_at(12*u.hourangle,inplace=True)
          am = secz(ha,c.dec,apo.lat)
          while True :
              # if length of sequence will bring us above airmass limit, trim number of exposures
              seq = Sequence(request['targname'],filt=['none'],n_exp=[request['no_exp']],t_exp=[request['exp_time']],camera=3,bin=2)
              length=seq.length()
              haend=ha+(length/3600.)*u.hourangle
              ham=hamax(c.dec,max_airmass,apo.lat).to(u.hourangle)
              dt=ham-haend
              if dt > 0 or request['no_exp'] < 1: break
              request['no_exp'] -= 1
          if request['no_exp'] < 1 : 
              if verbose: logger.info('{:d}, {:s} not enough time for an exposure: {:.2f}, hamax-haend: {:.2f}'.format(request['req_no'],request['targname'],am,dt))
              continue

          if am < 1.0 or am > max_airmass or dt<0 : 
              if verbose: logger.info('{:d}, {:s} out of airmass range {:.2f}, hamax-haend: {:.2f}'.format(request['req_no'],request['targname'],am,dt))
              continue
          if t < Time(request['start_window']) or t > Time(request['stop_window']) :
              if verbose: logger.info('{:d}, {:s} out of window {:s} {:s}'.format(request['req_no'],request['targname'],request['start_window'],request['stop_window']))
              if t > Time(request['stop_window']) :
                  # if we've passed observing window set status to abort  
                  load_song_status(request['req_no'],'abort')
              continue
          else :
              n_exp = request['no_exp'] if request['no_exp']*request['exp_time'] < dt_focus else int(dt_focus/request['exp_time'])
              logger.info('n_exp: {:d} no_exp: {:d}'.format(n_exp,request['no_exp']))
              if request['obs_mode'] == 'iodine' :
                  request=Table(request)
                  request['filter'] = [['iodine']]
                  request['bin'] = [[2]]
                  request['camera'] = [[3]]
                  request['n_exp'] = [[request[0]['no_exp']]]
                  request['t_exp'] = [[request[0]['exp_time']]]
              elif request['obs_mode'] == 'thar' or request['obs_mode'] == 'none-iodine' :
                  request=Table(request)
                  request['filter'] = [['thar','none','thar']]
                  request['bin'] = [[2,2,2]]
                  request['camera'] = [[3,3,3]]
                  request['n_exp'] = [[request[0]['no_thar_exp'],request[0]['no_exp'],request[0]['no_thar_exp']]]
                  request['t_exp'] = [[120,request[0]['exp_time'],120]]
              header={}
              header['OBSERVER'] = request['observer'][0]
              header['OBS-MODE'] = request['obs_mode'][0]
              header['PROJECT'] = request['project_name'][0]
              header['PROJ-ID'] = request['project_id'][0]
              header['REQ_NO'] = request['req_no'][0]
              pk=loadrequest('song_{:d}'.format(request['req_no'][0]),request['targname'][0],'song','song',priority)
              loadtarg(request['targname'][0],request['ra'][0],request['dec'][0],epoch=request['epoch'][0],mag=request['magnitude'][0]) 
              request['request_pk'] = [pk]
              return request[0], header
    return None, None

def getlocal(t=None, requests=None, site='APO', criterion='setting',mindec=-90,maxdec=90,skip=None,verbose=True) :
    """ 
    Get best request from local database table of requests and time
    """
    apo=EarthLocation.of_site(site)
    if t is None :
        t = Time(Time.now(),location=apo)
    if requests is None :
        requests = getrequests()

    lst=t.sidereal_time('mean')
    logger.info('Looking for request, lst: {:s}'.format(lst))
    if criterion == 'setting' : 
        tmin=12*u.hourangle
    elif criterion == 'longest' : 
        tmin=0.*u.hourangle
    elif criterion == 'best' : 
        tmin=12.*u.hourangle
    else : 
        logger.error('unknown criterion: ', criterion)
        return
    best=None
    am_best=-1
    dt_best=-1
    d=database.DBSession()
    priorities = sorted(set(requests['priority']),reverse=True)
    for priority in priorities :
      gd=np.where(requests['priority'] == priority) 
      for request in requests[gd] :
        if request['targname'] is None : continue
        if skip is not None and request['targname'] == skip : continue
        c=SkyCoord("{:s} {:s}".format(request['ra'],request['dec']),unit=(u.hourangle,u.deg))
        ha = t.sidereal_time('mean') - c.ra
        ha.wrap_at(12*u.hourangle,inplace=True)
        am = secz(ha,c.dec,apo.lat)

        seq = Sequence(request['targname'],filt=request['filter'],n_exp=request['n_exp'],
                       t_exp=request['t_exp'],camera=request['camera'],bin=request['bin'])
        length=seq.length()
        haend=ha+(length/3600.)*u.hourangle
        ham=hamax(c.dec,request['max_airmass'],apo.lat).to(u.hourangle)
        dt=ham-haend
        hamid=ha+(length/3600./2.)*u.hourangle
        if am < request['min_airmass'] or am > request['max_airmass'] or dt<0 : 
            if verbose: logger.info('{:s} out of airmass range {:.2f}'.format(request['targname'],am))
            continue
        #if c.dec<mindec*u.deg or c.dec>maxdec*u.deg : 
        #    logger.info('{:s} out of declination range'.format(request['targname']))
        #    continue
        obslist=d.query(sql='SELECT * from robotic.observed WHERE request_pk = {:d}'.format(request['request_pk']),fmt='list')
        maxfiles=0
        for o in obslist[1:] :
            maxfiles=np.max([maxfiles,len(o[3])])

        obs=Table(names=('observed_pk','request_pk', 'mjd', 'files'), dtype=('i4','i4','f4','{:d}S24'.format(maxfiles)))
        for o in obslist[1:] :
            while len(o[3]) < maxfiles : o[3].append('')
            try:
                for i,oo in enumerate(o[3]) :
                    if oo == None : o[3][i] = ''
                obs.add_row([o[0],o[1],o[2],o[3]])
            except:
                pdb.set_trace()

        if len(obs)>0:
            if request['nvisits'] > 0 and len(obs) > request['nvisits'] : 
                if verbose: logger.info('{:s} observed enough visits {:d}'.format(request['targname'],len(obs)))
                continue
            dt_visit = t.mjd - np.max(obs['mjd'])
            if dt_visit < request['dt_visit'] :
                if verbose: logger.info('{:s} observed too recently'.format(request['targname']))
                continue

        if (criterion == 'setting' and dt< tmin) :
            best = request
            tmin=ham-haend
            am_best=am
            dt_best=ham-haend
        elif (criterion=='longest' and dt>tmin)  :
            best = request
            tmin=ham-haend
            am_best=am
            dt_best=ham-haend
        elif (criterion == 'best' and abs(hamid)<tmin) :
            best = request
            tmin = abs(hamid)
            am_best=am
            dt_best=abs(hamid)
      # if we found one at this priority, break and continue
      if best is not None : break

    d.close()
    if best is None:
        logger.info('No good request available')
        header = None
    else :
        logger.info('request selected: {:s} {:s}'.format(best['targname'],best['sequencename']))
        header={}
        header['PROJECT'] = 'APO'
        header['REQ_PK'] = best['request_pk']
        header['TARGNAME'] = best['targname']
        header['SCHNAME'] = best['schedulename']
        header['SEQNAME'] = best['sequencename']
    return best, header

def observe_object(request,display=None,acquire=True,fact=1,nfact=1,header=None,req_no=-1) :
    """ Given request, do the observation and record
    """

    logger.info('Observing {:s}, acquire: {}'.format(request['targname'],acquire))
    try :
        # acquire target
        targ = Target(request['targname'],request['ra'],request['dec'],epoch=request['epoch'])
        if acquire : 
            success = targ.acquire(display=display)
            if not success : return False
    except KeyboardInterrupt :
        raise KeyboardInterrupt('CTRL-C')
    except:
        logger.exception('  failed acquire')
        t=Time.now()
        return False

    try :
        # observe requested sequence
        seq = Sequence(request['targname'],filt=request['filter'],bin=request['bin'],
                       n_exp=request['n_exp'],t_exp=request['t_exp'],camera=request['camera'])
        t=Time.now()
        names=seq.observe(targ.name,display=display,fact=fact,nfact=nfact,header=header,req_no=req_no)
        return load_object(request,t.mjd,names)
    except KeyboardInterrupt :
        raise KeyboardInterrupt('CTRL-C')
    except:
        logger.exception('  failed observe')
        return False

def load_song_status(req_no,status,no_exp=None) :
    """ Load status into song obs_request_4_status table
    """
    if req_no < 0 : return
    try :
        tab=Table()
        tab['req_no'] = [req_no]
        tab['ins_at'] = [datetime.datetime.now(datetime.timezone.utc)]
        tab['status'] = [status]
        if no_exp is not None : tab['no_exp'] = [no_exp]
        d=database.DBSession(host='song1m_db.apo.nmsu.edu',database='db_apo',user='song')
        d.update('public.obs_request_status_4',tab)
        d.close()
    except :
        logger.exception('  failed loading song_status')

def load_status(status) :
    """ Load status into database
    """
    try :
        tab=Table()
        tab['pk'] = [1]
        tab['status'] = [status]
        d=database.DBSession()
        d.update('robotic.status',tab)
        d.close()
    except :
        logger.exception('  failed loading status')

def load_object(request,mjd,names) :
    """ Load observation into database
    """
    # load observation into observed table
    logger.info('Loading observation {:s}'.format(request['targname']))
    try :
        obs=Table()
        obs['request_pk'] = [request['request_pk']]
        obs['mjd'] = [mjd]
        obs['files'] = [names]
        d=database.DBSession()
        d.ingest('robotic.observed',obs,onconflict='update')
        d.close()
    except:
        logger.exception('  failed loading')
        return False

    return True

def observe(focstart=32400,dt_focus=[0.5,1.0,1.0,2.0],display=None,dt_sunset=0,dt_nautical=-0.2,obs='apo',tz='US/Mountain',
        criterion='best',maxdec=None,cals=True,ccdtemp=-10, initfoc=True, fact=1, nfact=1, usesong=True) :
  """ Start full observing night sequence 

  Parameters
  ==========
  focstart : integer, default=32400
         initial focus guess
  dt_focus : float, default=[1.5]
         minimum time to wait after focus run before triggering another (will wait for sequence to complete). 
         if list, increment list index each time focus run is done, e.g. to achieve more frequent focus at beginning of night
  display : pyvista TV object
         if specified, display images as they are taken
  dt_sunset : float, default=0
         time relative to sunset to open dome (only if opening conditions are met)
  dt_nautical : float, default=-0.2
         time relative to nautical twilight to start observations
  obs : str, default='apo'
         observatory name, for getting site coordinates
  tz : str, default='US/Mountain'
         time zone
  criterion : str, default='best'
         criterion for choosing object to observe, 'setting', 'best' or 'longest'
  maxdec : float, default=None
         if given, maximum declination
  cals : bool, default=Ture
         if True take cals at end of night
  ccdtemp : float, default=-10
         set temperature for CCDs (guide and spectrograph)
  initfoc : boot, default=True
         True to take initial focus run, so can set False if restarting after focus has been done
  fact : float, default=1
         factor to increase all database exposure times by
  nfact : int, default=1
         factor to increase all database number of exposures b
  """
  global nautical_morn
  global foc0
  global nightlogger

  if not isinstance(dt_focus,list) : dt_focus=[dt_focus]

  # change if you want to try to run across multiple nights!
  night = 0
  while night < 1 :
    # new night
    night+=1
    guideok = True
    domeok = True
    telescopeok = True
    if not aposong.isguideok(guideok) :
        print('guider not responding, restart it!')
        return
    if not aposong.isdomeok(domeok) :
        print('dome not responding, restart it!')
        return
    if not aposong.istelescopeok(telescopeok) :
        print('telescope not responding, restart it!')
        return

    nightlogger=logging.getLogger('night_logger')
    nightlogger_string = io.StringIO()
    handler = logging.StreamHandler(nightlogger_string)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    nightlogger.addHandler(handler)

    nightlogger.info('Starting night!')
    logger.info('Starting night!')

    site=Observer.at_site(obs,timezone=tz)
    sunset =site.sun_set_time(Time.now(),which='nearest')
    sunrise =site.sun_rise_time(Time.now(),which='next')
    nautical = site.twilight_evening_nautical(Time.now(),which='nearest')
    nautical_morn = site.twilight_morning_nautical(Time.now(),which='next')

    # open MJD file
    f=open('{:d}'.format(int(nautical.mjd)),'w')
    f.close()

    # set CCD temperatures
    aposong.settemp(ccdtemp,cam=0)
    aposong.settemp(ccdtemp,cam=3)

    # open dome when safe after desired time relative to sunset
    opentime = sunset+dt_sunset*u.hour
    load_status('closed')
    while (Time.now()-opentime)<0 :
        # wait until sunset + dt_sunset hours 
        logger.info('waiting for sunset+dt_sunset: {:.3f} '.format((opentime-Time.now()).to(u.hour).value,' hours'))
        time.sleep(60)

    while not aposong.issafe() and (Time.now()-nautical_morn).to(u.hour) < 0*u.hour : 
        # wait until safe to open based on Safety
        logger.info('waiting for issafe()')
        time.sleep(30)

    # if we haven't opened by morning twilight, we're done for the night!
    if (Time.now()-(nautical_morn+dt_nautical*u.hour)).to(u.hour) > 0*u.hour : continue

    # open if not already open (e.g., from previous robotic.observe invocation same night)
    if aposong.D.ShutterStatus != 0 :
        aposong.domeopen()
        load_status('open')
    logger.info('open at: {:s}'.format(Time.now().to_string()))
    nightlogger.info('open at: {:s}'.format(Time.now().to_string()))

    # wait for sunset to open louvers
    while (Time.now()-sunset).to(u.hour) < 0*u.hour :
        logger.info('waiting for sunset for louvers: {:.3f} '.format((sunset-Time.now()).to(u.hour).value,' hours'))
        time.sleep(60)
    aposong.louvers(True)

    # cals
    #if cals :
    #    cal.cals()

    # wait for nautical twilight
    while (Time.now()-(nautical+dt_nautical*u.hour)).to(u.hour) < 0*u.hour :
        try :
            # if dome has been closed by 3.5m but can now be opened again, do it
            if aposong.issafe() and aposong.D.ShutterStatus != 0 :
                logger.info('reopening dome')
                aposong.domeopen()
                load_status('open')
            elif aposong.D.ShutterStatus != 0 :
                load_status('closed')
            logger.info('waiting for nautical twilight+dt_nautical: {:.3f}'.format(
                        (nautical+dt_nautical*u.hour-Time.now()).to(u.hour).value,' hours'))
            time.sleep(60)
        except KeyboardInterrupt :
            return

    # in case 3.5m has closed while we were waiting....
    if aposong.D.ShutterStatus != 0 :
        aposong.domeopen()
        load_status('open')

    # fans off for observing?
    #aposong.fans_off()

    # focus star on meridian 
    if initfoc : 
        load_status('focus')
        foc0=focus(foc0=focstart,delta=75,n=15,display=display,iodine=False)
    else :
        foc0=focstart
    foctime=Time.now()

    oldtarg=''
    skiptarg=None
    closed = False
    nfocus = 0

    # loop for objects until morning nautical twilight
    while (Time.now()-nautical_morn).to(u.hour) < 0*u.hour : 
      load_status('ready')
      try :
        tnow=Time.now()
        logger.info('nautical twilight in : {:.3f}'.format((nautical_morn-tnow).to(u.hour).value))

        # close if not safe, open if it has become safe
        if not aposong.issafe() : 
            if not closed :
                logger.info('closing: {:s}'.format(Time.now().to_string()))
                nightlogger.info('closing: {:s}'.format(Time.now().to_string()))
                aposong.domeclose()
                load_status('closed')
                closed = True
            oldtarg=''
            time.sleep(90)
            continue
        elif aposong.D.ShutterStatus != 0 :
            logger.info('opening: {:s}'.format(Time.now().to_string()))
            nightlogger.info('opening: {:s}'.format(Time.now().to_string()))
            aposong.domeopen()
            load_status('open')
            closed = False

        # check if guider, dome, telescope are responding, send text if not
        guideok = aposong.isguideok(guideok,loggers=[logger,nightlogger],recipients=aposong.config['test_recipients'])
        domeok = aposong.isdomeok(domeok,loggers=[logger,nightlogger],recipients=aposong.config['test_recipients'])
        telescopeok = aposong.istelescopeok(telescopeok,loggers=[logger,nightlogger],recipients=aposong.config['test_recipients'])

        logger.info('tnow-foctime : {:.3f}'.format((tnow-foctime).to(u.hour).value))
        if (tnow-foctime).to(u.hour) > dt_focus[nfocus]*u.hour :
            # focus if it has been more than dt_focus since last focus
            nfocus = nfocus+1 if nfocus+1<len(dt_focus) else len(dt_focus)-1
            aposong.guide('stop')
            #foc0=focus(foc0=foc0,display=display,decs=[90,85,75,65,55,40],iodine=False)
            load_status('focus')
            foc0=focus(foc0=foc0,display=display,iodine=False)
            foctime=tnow
            oldtarg=''
        else :
            # observe best object!
            if usesong : best,header=getsong()
            else : best=None
            if best is None :
                best,header=getlocal(criterion=criterion,maxdec=maxdec,skip=skiptarg)
                if best is not None :
                    nightlogger.info('observe: LOCAL, {:s}'.format(best['targname']))
                req_no = -1
            else :
                req_no = best['req_no']
                nightlogger.info('observe: SONG {:d}, {:s}, {:s}'.format(req_no, best['targname'], best['project_name']))
            if best is None :
                time.sleep(300)
            else :
                nightlogger.info('observe: {:s}'.format(best['targname']))
                logger.info('observe: {:s}'.format(best['targname']))
                load_status('observing')
                load_song_status(req_no,'exec',no_exp=0)
                success = observe_object(best,display=display,acquire=(best['targname']!=oldtarg),
                                         fact=fact,nfact=nfact,req_no=req_no,header=header)
                nightlogger.info('success: {:d}'.format(success))
                logger.info('success: {:d}'.format(success))
                # if object failed, skip it for the next selection, but can try again after that
                if success : 
                    oldtarg=best['targname'] 
                    skiptarg=None
                    try : 
                        if req_no > 0 :
                            d=database.DBSession(host='song1m_db.apo.nmsu.edu',database='db_song',user='song')
                            best = d.query(sql="SELECT * FROM public.obs_request_4 as req JOIN public.obs_request_status_4 as stat ON req.req_no = stat.req_no WHERE req.req_no = {:d}".format(req_no), fmt='table')
                            d.close()
                            logger.info('no_exp: {:d} no_exp2: {:d}'.format(best[0]['no_exp'],best[0]['no_exp2']))
                            load_song_status(best[0]['req_no'],'done')
                    except KeyError : pass
                else : 
                    skiptarg=best['targname']
                    try: 
                        if req_no > 0 :
                            d=database.DBSession(host='song1m_db.apo.nmsu.edu',database='db_song',user='song')
                            best = d.query(sql="SELECT * FROM public.obs_request_4 as req JOIN public.obs_request_status_4 as stat ON req.req_no = stat.req_no WHERE req.req_no = {:d}".format(req_no), fmt='table')
                            d.close()
                            load_song_status(best[0]['req_no'],'abort')
                    except KeyError : pass
      except KeyboardInterrupt :
          return
 
      except:
          logger.exception('observe error!')

    # close
    try: aposong.guide('stop')
    except: pass
    logger.info('closing for night: {:s}'.format(Time.now().to_string()))
    nightlogger.info('closing for night: {:s}'.format(Time.now().to_string()))
    aposong.domeclose()
    load_status('closed')

    # morning cals
    if cals :
        header={}
        header['OBS-MODE'] = 'cal'
        header['PROJECT'] = 'No-Proj'
        header['PROJ-ID'] = 'No-Proj'
        cal.cals(header=header,flats=3,iodineflats=3)

    # night log
    mkhtml()

    logger.info('completed observing loop!')
    nightlogger.info('completed observing loop!')
    y,m,d,hr,mi,se = Time.now().ymdhms
    ut = 'UT{:d}{:02d}{:02d}'.format(y-2000,m,d)
    message='<A HREF=http://astronomy.nmsu.edu/1m/data/'+ut+'/'+ut+'.html>'+ut+' astronomy.nmsu.edu web page</A><P>'
    message+=nightlogger_string.getvalue()
    attachments = ['/data/1m/logs/daily.log','/data/1m/'+ut+'/focus.png','/data/1m/'+ut+'/throughput.png']
    mail.send(aposong.config['mail_recipients'],
              subject='APO SONG observing {:d} completed successfully'.format(int(Time.now().mjd)),
              message=message.replace('\n','<br>'), snapshot=True,attachment=attachments,html=True)
    subprocess.run('copy {:s}'.format(ut),shell=True)

    try :os.remove('{:d}'.format(int(nautical.mjd)))
    except : print('no MJD file to remove?')

    # restart display
    try :
        if display is not None :
            matplotlib.use('TkAgg') 
            plt.close('all')
            print('you will need to restart display window with disp=disp_init() if you use display=disp in robotic.observe()')
    except NameError : pass

def focus(foc0=28800,delta=75,n=9,decs=[52],iodine=True,display=None) :
    """ Do focus run for object on meridian
    """    

    for dec in decs :
        t=Time(Time.now(), location=EarthLocation.of_site('APO'))
        lst=t.sidereal_time('mean').value
        aposong.usno(ra=lst,dec=dec,magmin=9,magmax=10,rad=5*u.degree)

        # focus run with iodine cell (if desired)
        if iodine :
            foc0_0=copy.copy(foc0)
            aposong.iodine_in()
            f,focvals,best=aposong.focrun(foc0-4625,delta,n,exptime=3,filt=None,bin=1,thresh=100,display=display,max=5000)
            alt = aposong.T.Altitude
            logger.info('iodine focus: {:d}  besthf: {:.2f} foc0: {:d}  alt: {:.1f}'.format(f,best,foc0,alt))
            nightlogger.info('iodine focus: {:d}  besthf: {:.2f} foc0: {:d}  alt: {:.1f}'.format(f,best,foc0,alt))
            while f<focvals[2] or f>focvals[-3] :
                # redo if best focus is near end of run
                if (Time.now()-nautical_morn).to(u.hour) > 0*u.hour or not aposong.issafe(): 
                    break
                foc0=copy.copy(f)
                f,focvals,best=aposong.focrun(foc0-delta,delta,n+1,exptime=3,filt=None,bin=1,thresh=100,display=display,max=5000)
                alt = aposong.T.Altitude
                logger.info('iodine focus: {:d}  foc0: {:d}  alt: {:.1f}'.format(f,best,foc0,alt))
                nightlogger.info('iodine focus: {:d}  foc0: {:d}  alt: {:.1f}'.format(f,best,foc0,alt))
            aposong.iodine_out()
            foc0=copy.copy(foc0_0)

        # normal focus run
        f,focvals,best=aposong.focrun(foc0,delta,n,exptime=3,filt='open',bin=1,thresh=100,display=display,max=5000)
        alt = aposong.T.Altitude
        logger.info('focus: {:d}  besthf: {:.2f} foc0: {:d}  alt: {:.1f}'.format(f,best,foc0,alt))
        nightlogger.info('focus: {:d}  besthf: {:.2f} foc0: {:d}  alt: {:.1f}'.format(f,best,foc0,alt))
        while f<focvals[2] or f>focvals[-3] :
            # redo if best focus is near end of run
            if (Time.now()-nautical_morn).to(u.hour) > 0*u.hour or not aposong.issafe(): 
                break
            foc0=copy.copy(f)
            f,focvals,best=aposong.focrun(foc0,delta,n,exptime=3,filt=None,bin=1,thresh=100,display=display,max=5000)
            alt = aposong.T.Altitude
            logger.info('focus: {:d}  besthf: {:.2f} foc0: {:d}  alt: {:.1f}'.format(f,best,foc0,alt))
            nightlogger.info('focus: {:d}  besthf: {:.2f} foc0: {:d}  alt: {:.1f}'.format(f,best,foc0,alt))

    return f
       
def test() :
    t = Time.now()
    while True :
        best,am,dt = getlocal(t)
        print(t,am,dt)
        print(best)
        tab.remove_row
        t += 5400*u.s

def loadtargs(file,schedule='rv',sequence='UBVRI',insert=False) :
    """ Load database tables from old 1m input file
    """
    fp=open(file,'r')
    lines=fp.readlines()
    requests=[]
    for iline,line in enumerate(lines[2::2]) :
        name,ra,dec,epoch,amin,amax,y,m,d,m,d,umax,utmin,utmax,hmin,hmax,pmax,dmin,nv,dt,foc,df,g,gloc,skip=line.split()
        targ=Target(name,ra,dec)
        targs.append(targ.table())
        c=SkyCoord("{:s} {:s}".format(ra,dec),unit=(u.hourangle,u.deg))
        priority=20
        if iline == 0 :
            rtab = Table( [[name], [schedule],[sequence]], names = ('targname','schedulename','sequencename','priority') )
        else :
            rtab.add_row([name,'rv','UBVRI',priority])
    tab=vstack(targs)
    if insert :
        tab.rename_column('name','targname')
        d=database.DBSession()
        d.ingest('robotic.target',tab,onconflict='update')
        d.ingest('robotic.request',rtab,onconflict='update')
        d.close()

def loadtarg(targname,ra,dec,epoch=2000,mag=99) :
    """ Load a single target into robotic.target
    """
    tab=Table()
    tab['targname'] = [targname]
    tab['ra'] = [ra]
    tab['dec'] = [dec]
    tab['epoch'] = [epoch]
    tab['mag'] = [mag]
    d=database.DBSession()
    d.ingest('robotic.target',tab,onconflict='update')
    d.close()

def loadsched(name,min_airmass=1.0,max_airmass=2,nvisits=1,dt_visit=0) :
    """ Load a single schedule into robotic.schedule
    """
    d=database.DBSession()
    tab=Table()
    d=database.DBSession()
    schedule=Schedule(name,min_airmass=float(min_airmass),max_airmass=float(max_airmass),nvisits=int(nvisits),dt_visit=float(dt_visit))
    d.ingest('robotic.schedule',schedule.table(),onconflict='update')
    d.close()

def loadseq(name,t_exp=[1],n_exp=[1],filt=['V'],camera=[0],bin=[1]) :
    """ Load a single sequence into robotic.sequence
    """
    for var in t_exp, n_exp, filt, camera, bin :
        if len(var) < 5 :
            for i in range(len(var),5) :
                if isinstance(var[0],str) :
                    var.extend([' '])
                else :
                    var.extend([0])
    d=database.DBSession()
    sequence=Sequence(name,filt=filt,t_exp=t_exp,n_exp=n_exp,camera=camera,bin=bin)
    d.ingest('robotic.sequence',sequence.table(),onconflict='update')
    d.close()
    return sequence.table()

def loadrequest(name,targname,seqname,schedname,priority) :
    d=database.DBSession()
    tab=Table()
    tab['targname'] = [targname]
    tab['sequencename'] = [seqname]
    tab['schedulename'] = [schedname]
    tab['priority'] = [priority]
    d.ingest('robotic.request',tab,onconflict='update')
    pk=d.query(sql="select currval('robotic.request_request_pk_seq')")[0][0]
    d.close()
    return pk

def createrequest(targname,ra,dec,epoch=2000,filter=['none'],bin=2,n_exp=[1],t_exp=[60],camera=3) :
    request={}
    request['targname'] = targname
    request['ra'] = ra
    request['dec'] = dec
    request['epoch'] = epoch
    request['filter'] = filter if isinstance(filter,list) else [filter]
    request['n_exp'] = n_exp if isinstance(n_exp,list) else [n_exp]
    request['t_exp'] = t_exp if isinstance(t_exp,list) else [t_exp]
    request['bin'] = bin if isinstance(bin,list) else [bin]
    request['camera'] = camera if isinstance(camera,list) else [camera]
    return request

def mkhtml(mjd=None) :
    """ Make HTML pages for a night of observing

    Parameters
    ----------
    mjd : int, default=None
          make pages for specified MJD, now if mjd=None
    """
    if mjd is None : 
        mjd = int(Time.now().mjd)
    mkmovie(mjd)
    mkfocusplots(mjd,clobber=False)
    mklog(mjd,clobber=False)
    y,m,d,hr,mi,se = Time(mjd,format='mjd').ymdhms
    ut = 'UT{:d}{:02d}{:02d}'.format(y-2000,m,d)
    dofocus.mksum(mjd,hard='/data/1m/'+ut+'/focus.png')
    out=reduce.throughput_all(mjd=mjd,hard='/data/1m/'+ut+'/throughput.png')

def mkmovie(mjd,root='/data/1m/',clobber=False) :
    """ Make guider movies from guide images in guide subdirectory for specified MJD
    """
    y,m,d,hr,mi,se = Time(mjd,format='mjd').ymdhms
    ut = 'UT{:d}{:02d}{:02d}'.format(y-2000,m,d)

    dir=root+ut+'/guide'
    red=imred.Reducer(dir=dir)
    files=sorted(glob.glob(dir+'/acquire*.fits'))
    ims=[]
    seqs=[]
    for file in files :
        im=int(file.split('.')[-2])
        ims.append(im)
    for i,im in enumerate(ims[1:]) :
        if ims[i+1]-ims[i] > 1 :
            seqs.append([ims[i]+1,ims[i+1]])

    print(ims)
    print(seqs)
    grid=[]
    row=[]
    for i,seq in enumerate(seqs):
        try: plt.close('all')
        except: pass
        t=tv.TV()
        out='{:s}/{:d}.mp4'.format(dir,seq[0])
        if clobber or not os.path.isfile(out) :
            #red.movie(range(seq[0],seq[1]),display=t,max=10000,out=out)
            try : red.movie(range(seq[0],seq[1]),display=t,max=10000,out=out)
            except: pdb.set_trace()

        if i>0 and i%5 == 0 :
            grid.append(row)
            row=[]
        row.append('guide/{:d}.mp4'.format(seq[0]))
    for ii in range(i%5+1,5) :
        row.append('')
    grid.append(row)
    html.htmltab(grid,file=root+ut+'/guide.html',video=True)

def mkfocusplots(mjd,display=None,root='/data/1m/',clobber=False) :
    """ Make focus plot from focus sequences from database for specified MJD
    """
    d=database.DBSession()
    out=d.query('obs.focus',fmt='list')
    d.close()
    i =np.where(np.array(out[0]) == 'mjd')[0][0]
    ifile =np.where(np.array(out[0]) == 'files')[0][0]
    mjds=[]
    files=[]
    for o in out[1:] : 
        mjds.append(o[i])
        files.append(o[ifile])

    gd=np.where(np.array(mjds).astype(int) == mjd)[0]
    sort=np.argsort(np.array(mjds)[gd])

    grid=[]
    row=[]
    i=0
    for seq in gd[sort] :
        print(files[seq])
        try : 
            if clobber or not os.path.isfile(root+files[seq][0].replace('.fits','.png')) :
                dofocus.focus(files[seq],display=display,root=root,plot=True,hard=True,
                              pixscale=aposong.pixscale(0))
            if i>0 and i%5 == 0 :
                grid.append(row)
                row=[]
            row.append(os.path.basename(files[seq][0].replace('.fits','.png')))
            i+=1
        except : 
            continue
    i-=1
    for ii in range(i%5+1,5) :
        row.append('')
    grid.append(row)
    dir=files[seq][0].split('/')[0]
    html.htmltab(grid,file=root+dir+'/focus.html',size=250)

def mklog(mjd,root='/data/1m/',pause=False,clobber=False) :
    """ Makes master log page for specified MJD with observed table, exposure table, and links
    """
    y,m,d,hr,mi,se = Time(mjd,format='mjd').ymdhms
    ut = 'UT{:d}{:02d}{:02d}'.format(y-2000,m,d)

    d=database.DBSession()
    out=d.query(sql='select * from obs.exposure where mjd>{:d} and mjd<{:d} and camera=3'.format(mjd,mjd+1),fmt='table')
    guideout=d.query(sql='select * from obs.exposure where mjd>{:d} and mjd<{:d} and camera=0'.format(mjd,mjd+1),fmt='table')
    obslist=d.query(sql='''
          select * from robotic.observed as obs 
          join robotic.request as req on obs.request_pk = req.request_pk 
          where mjd>{:d} and mjd<{:d}'''.format(mjd,mjd+1),fmt='list')
    d.close()
    for o in [out,guideout] :
        o['mjd'].format='{:.3f}'
        o['exptime'].format='{:.2f}'
        o['alt'].format='{:.2f}'
        #o['ra'].format='{:.6f}'
        #o['dec'].format='{:.6f}'
        o['ccdtemp'].format='{:.2f}'
    out['airmass'] = 1./np.cos((90-out['alt'])*np.pi/180.)
    out['airmass'].format='{:.2f}'

    maxfiles=50
    #for o in obslist[1:] :
    #    maxfiles=np.max([maxfiles,len(o[3])])

    obs=Table()
    obs=Table(names=('request_pk', 'mjd', 'targname', 'schedulename','sequencename','priority','files'), dtype=('i4','f4', 'S', 'S','S','i4','{:d}S24'.format(maxfiles)))
    for o in obslist[1:] :
        while len(o[3]) < maxfiles : o[3].append('')
        try: obs.add_row([o[1],o[2],o[5],o[6],o[7],o[8],o[3]])
        except: pass

    j = np.argsort(obs['mjd'])
    obs = obs[j]
    obs['request'] = ' '*80
    red=imred.Reducer('SONG')

    # do ThAr first
    wavs=[]
    for req,o in enumerate(obs) :
        for i,f in enumerate(o['files']) : 
            if f.find(b'thar') >=0 :
                imec=reduce.specreduce(f.decode(),red=red,clobber=clobber,write=True)
                file=imec.header['FILE'].split('.')
                outfile='{:s}/{:s}_wav.{:s}.fits'.format(os.path.dirname(f.decode()).replace('1m/','1m/reduced/'),file[0],file[-2])
                wavs.append(spectra.WaveCal(outfile))

    for req,o in enumerate(obs) :
        fig,ax=plots.multi(1,3,figsize=(8,4),hspace=0.001)
        for i,f in enumerate(o['files']) : 
          if f != b'' and f.find( b'thar') <0 :
            try :
                im=red.rd(f.decode())
                imec=reduce.specreduce(f.decode(),red=red,clobber=clobber,write=True,wav=wavs)
                if f.find(b'thar')  >= 0 : continue
                t,sn,t_orders,sn_orders=reduce.throughput(imec,ax,os.path.basename(f.decode()),red=red,orders=[34])
                o['files'][i] = os.path.basename(o['files'][i])
                ind=np.where(np.char.find(out['file'],os.path.basename(f.decode()))>=0)[0]
                reduced=Table()
                reduced['request_pk'] = [o['request_pk']]
                reduced['exp_pk'] = [out[ind]['exp_pk'][0]]
                reduced['throughput'] = [t]
                reduced['sn'] = [sn]
                reduced['throughput_orders'] = [t_orders]
                reduced['sn_orders'] = [sn_orders]
                reduced['traceshift'] = [imec.header['T_SHIFT']]
                d=database.DBSession()
                d.ingest('obs.reduced',reduced,onconflict='update')
                d.close()
            except :
                print('error with ',f)
                pdb.set_trace()
                pass
        for i in range(3) :
            lim = ax[i].get_ylim()
            ax[i].set_ylim(0,lim[1])
        ax[0].legend(fontsize='xx-small')
        fig.savefig('{:s}/{:s}/{:s}_{:d}.png'.format(root,ut,o['targname'].replace(' ','_'),req))
        o['request'] = '<A HREF={:s}_{:d}.png> {:s} </A>'.format(o['targname'].replace(' ','_'),req,o['targname'])
        if pause : pdb.set_trace()
        plt.close()

    for o in [out,guideout] :
        gd = np.where(o['file'] != '')[0]
        o=o[gd]
        j = np.argsort(o['dateobs'])
        o=o[j]
    # not sure why following but seems to be necessary
    j=np.argsort(out['dateobs'])
    out=out[j]

    fp = open('{:s}/{:s}/{:s}.html'.format(root,ut,ut),'w')

    fp.write('<HTML><BODY>\n')
    fp.write('<h2>{:s}  MJD: {:d}</h2><br>'.format(ut,mjd))
    fp.write('<A HREF=guide.html>Guider movies</A><BR>\n')
    fp.write('<A HREF=focus.html>Focus curves</A><BR>\n')
    fp.write('<A HREF=../reduced/{:s}>Reduced data</A><BR>\n'.format(ut))
    fp.write('<A HREF=focus.png><IMG SRC=focus.png WIDTH=30%></A>\n')
    fp.write('<A HREF=throughput.png><IMG SRC=throughput.png WIDTH=30%></A>\n')
    fp.write('<A HREF=https://irsc.apo.nmsu.edu/graphArchive/{:d}graph.png><IMG SRC=https://irsc.apo.nmsu.edu/graphArchive/{:d}graph.png WIDTH=30%></A>'.format(mjd,mjd))

    fp.write('<p>Observed robotic requests: <BR>\n')
    obs['files'] = obs['files'].astype('<U')
    tab = html.tab(obs['request','mjd','schedulename','sequencename','priority','files'],file=fp)

    fp.write('<p>Spectrograph exposures: <BR>\n')
    tab = html.tab(out['file','dateobs','ra','dec','exptime','airmass','camera','filter','focus','ccdtemp'],file=fp)

    fp.write('<p>Guider exposures: <BR>\n')
    tab = html.tab(guideout['file','dateobs','ra','dec','exptime','camera','filter','focus','ccdtemp'],file=fp)
    fp.write('</BODY></HTML>\n')
    fp.close()

    try : os.symlink('{:s}.html'.format(ut),'{:s}/{:s}/0_{:s}.html'.format(root,ut,ut))
    except: pass

    return out
