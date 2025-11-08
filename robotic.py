from astropy.coordinates import SkyCoord,EarthLocation
from astropy.time import Time, TimeDelta
from astropy.table import Table, hstack, vstack
from astropy.io import fits
import astropy.units as u
from astroplan import Observer, time_grid_from_range
import matplotlib
import matplotlib.pyplot as plt
from pyvista import imred,tv
from holtztools import html, plots

import copy
import glob
import os
import pdb
import numpy as np
import time
import threading

import logging
import yaml
import logging.config
 
import focus as dofocus
import reduce

try:
    with open('logging.yml', 'rt') as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)
except:
    print("can't open logging.yml")

logger=logging.getLogger(__name__)

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
        aposong.newguider('stop')
        aposong.slew(self.ra,self.dec)
        aposong.newguider('start')
        time.sleep(30)

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

    def observe(self,name,display=None,fact=1,nfact=1) :
        names = []
        for filt,nexp,texp,cam,bin in zip(self.filt,self.n_exp,self.t_exp,self.camera,self.bin) :
            for iexp in range(nexp*nfact) :
                logger.info('Expose camera: {:d} exptime: {:.2f} bin: {:d}, '.format(cam,texp,bin))
                if filt == 'thar' :
                    aposong.newguider('pause')
                    time.sleep(10)
                    cal.cals(thar=nexp,flats=0,thar_exptime=texp)
                    names.append('thar')
                    aposong.newguider('resume')
                    time.sleep(30)
                else :
                    exp=aposong.expose(texp*fact,filt,name=name,display=display,cam=cam,bin=bin,targ=self.name)
                    names.append(exp.name)
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

def getbest(t=None, requests=None, site='APO', criterion='setting',mindec=-90,maxdec=90,skip=None,verbose=True) :
    """ 
    Get best request given table of requests and time
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
        obs=d.query(sql='SELECT * from robotic.observed WHERE request_pk = {:d}'.format(request['request_pk']))
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
    else :
        logger.info('request selected: {:s} {:s}'.format(best['targname'],best['sequencename']))
    return best

def observe_object(request,display=None,acquire=True,fact=1,nfact=1) :
    """ Given request, do the observation and record
    """

    logger.info('Observing {:s}, acquire: {}'.format(request['targname'],acquire))
    try :
        # acquire target
        targ = Target(request['targname'],request['ra'],request['dec'],epoch=request['epoch'])
        if acquire : 
            targ.acquire(display=display)
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
        names=seq.observe(targ.name,display=display,fact=fact,nfact=nfact)
        return load_object(request,t.mjd,names)
    except KeyboardInterrupt :
        raise KeyboardInterrupt('CTRL-C')
    except:
        logger.exception('  failed observe')
        return False

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

def obsopen(opentime) :
    """ Open observatory at/after requested time and when safe 
    """
    while (Time.now()-opentime)<0 :
        # wait until sunset + dt_sunset hours 
        logger.info('waiting for sunset: {:.3f} '.format((opentime-Time.now()).to(u.hour).value,' hours'))
        time.sleep(60)

    while not aposong.issafe() :
        # wait until safe to open based on Safety
        time.sleep(30)

    # open if not already open (e.g., from previous observe invocation same night)
    if aposong.D.ShutterStatus != 0 :
        aposong.domeopen()
    logger.info('open at: {:s}'.format(Time.now().to_string()))

def observe(foc0=32400,dt_focus=1.5,display=None,dt_sunset=0,dt_nautical=-0.2,obs='apo',tz='US/Mountain',
            criterion='best',maxdec=None,cals=True,ccdtemp=-10, initfoc=True, fact=1, nfact=1, override=0) :
  """ Start full observing night sequence 

  Parameters
  ==========
  foc0 : integer, default=32400
         initial focus guess
  dt_focus : float, default=1.5
         minimum time to wait after focus run before triggering another (will wait for sequence to complete)
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
  override : int, default=0
         DANGEROUS: if given, set override for specified number of seconds to bypass need for 2.5m or 3.5m to be open
  """
  while True :
    site=Observer.at_site(obs,timezone=tz)
    sunset =site.sun_set_time(Time.now(),which='nearest')
    sunrise =site.sun_rise_time(Time.now(),which='next')
    nautical = site.twilight_evening_nautical(Time.now(),which='nearest')
    nautical_morn = site.twilight_morning_nautical(Time.now(),which='next')

    # open dome when safe after desired time relative to sunset
    if override > 0 :
        if override > 1800 :
            input('you are requesting an override of >1800s, confirm with keystrok: ')
        aposong.S.Action('override',override)
        override=0

    aposong.settemp(ccdtemp,cam=0)
    aposong.settemp(ccdtemp,cam=3)
    obsopen(sunset+dt_sunset*u.hour)

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
            logger.info('waiting for nautical twilight: {:.3f}'.format(
                        (nautical+dt_nautical*u.hour-Time.now()).to(u.hour).value,' hours'))
            time.sleep(60)
        except KeyboardInterrupt :
            return

    # in case 3.5m has closed while we were waiting....
    if aposong.D.ShutterStatus != 0 :
        aposong.domeopen()

    # fans off for observing?
    #aposong.fans_off()

    # focus star on meridian 
    if initfoc : foc=focus(foc0=foc0,delta=75,n=15,display=display)
    foctime=Time.now()

    oldtarg=''
    skiptarg=None
    while (Time.now()-nautical_morn).to(u.hour) < 0*u.hour : 
      try :
        tnow=Time.now()
        logger.info('nautical twilight in : {:.3f}'.format((nautical_morn-tnow).to(u.hour).value))
        if not aposong.issafe() : 
            aposong.domeclose()
            time.sleep(90)
            continue
        elif aposong.D.ShutterStatus != 0 :
            aposong.domeopen()

        logger.info('tnow-foctime : {:.3f}'.format((tnow-foctime).to(u.hour).value))
        if (tnow-foctime).to(u.hour) > dt_focus*u.hour :
            aposong.newguider('stop')
            foc0=focus(foc0=foc0,display=display)
            foctime=tnow
            oldtarg=''
        else :
            best=getbest(criterion=criterion,maxdec=maxdec,skip=skiptarg)
            if best is None :
                time.sleep(60)
            else :
                success = observe_object(best,display=display,acquire=(best['targname']!=oldtarg),fact=fact,nfact=nfact)
                # if object failed, skip it for the next selection, but can try again after that
                if success : 
                    oldtarg=best['targname'] 
                    skiptarg=None
                else : 
                    skiptarg=best['targname']
      except KeyboardInterrupt :
          return
 
      except:
          logger.exception('observe error!')

    logger.info('closing: {:s}'.format(Time.now().to_string()))
    aposong.domeclose()

    if cals :
        cal.cals()

    matplotlib.use('Agg') 
    #mkhtml()
    html_process=threading.Thread(target=mkhtml)
    html_lock=threading.Lock() 
    html_process.start()
    html_process.join()

    input("Hit a key to start next night: ")

def focus(foc0=28800,delta=75,n=9,display=None) :
    """ Do focus run for object on meridian
    """    
    t=Time(Time.now(), location=EarthLocation.of_site('APO'))
    lst=t.sidereal_time('mean').value
    aposong.usno(ra=lst,dec=50.,magmin=9,magmax=10,rad=5*u.degree)

    foc0_0=copy.copy(foc0)
    aposong.iodine_in()
    f=aposong.focrun(foc0-4625,delta,n,exptime=3,filt=None,bin=1,thresh=100,display=display,max=5000)
    while f<foc0-2.1*delta or f>foc0+2.1*delta :
        logger.info('focus: {:d}  foc0: {:d}'.format(f,foc0))
        foc0=copy.copy(f)
        f=aposong.focrun(foc0-delta,delta,n+1,exptime=3,filt=None,bin=1,thresh=100,display=display,max=5000)
    aposong.iodine_out()
    foc0=copy.copy(foc0_0)

    f=aposong.focrun(foc0-delta,delta,n+1,exptime=3,filt='open',bin=1,thresh=100,display=display,max=5000)
    while f<foc0-2.1*delta or f>foc0+2.1*delta :
        logger.info('focus: {:d}  foc0: {:d}'.format(f,foc0))
        foc0=copy.copy(f)
        f=aposong.focrun(foc0,delta,n,exptime=3,filt=None,bin=1,thresh=100,display=display,max=5000)
    return f
       
def test() :
    t = Time.now()
    while True :
        best,am,dt = getbest(t)
        print(t,am,dt)
        print(best)
        tab.remove_row
        t += 5400*u.s

def loadtarg(file,schedule='rv',sequence='UBVRI',insert=False) :
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

def loadsched(name,min_airmass=1.0,max_airmass=2,nvisits=1,dt_visit=0) :
    d=database.DBSession()
    schedule=Schedule(name,min_airmass=float(min_airmass),max_airmass=float(max_airmass),nvisits=int(nvisits),dt_visit=float(dt_visit))
    d.ingest('robotic.schedule',schedule.table(),onconflict='update')
    d.close()

def loadseq(name,t_exp=[1],n_exp=[1],filt=['V'],camera=[0],bin=1) :
    d=database.DBSession()
    sequence=Sequence(name,filt=filt,t_exp=t_exp,n_exp=n_exp,camera=camera,bin=bin)
    d.ingest('robotic.sequence',sequence.table(),onconflict='update')
    d.close()

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
    mkfocusplots(mjd)
    mklog(mjd)

def mkmovie(mjd,root='/data/1m/',clobber=False) :
    """ Make guider movies from guide images in guide subdirectory for specified MJD
    """
    y,m,d,hr,mi,se = Time(mjd,format='mjd').ymdhms
    ut = 'UT{:d}{:02d}{:02d}'.format(y-2000,m,d)

#    matplotlib.use('Agg')
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
            try : red.movie(range(seq[0],seq[1]),display=t,max=10000,out=out)
            except: continue
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
#    matplotlib.use('Agg')
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

def mklog(mjd,root='/data/1m/',pause=False) :
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

    maxfiles=0
    for o in obslist[1:] :
        maxfiles=np.max([maxfiles,len(o[3])])

    obs=Table()
    obs=Table(names=('mjd', 'targname', 'schedulename','sequencename','priority','files'), dtype=('f4', 'S', 'S','S','i4','{:d}S24'.format(maxfiles)))
    out['mjd'].format='{:.3f}'
    guideout['mjd'].format='{:.3f}'
    for o in obslist[1:] :
        while len(o[3]) < maxfiles : o[3].append('')
        try: obs.add_row([o[2],o[5],o[6],o[7],o[8],o[3]])
        except: pdb.set_trace()

    j = np.argsort(obs['mjd'])
    obs = obs[j]
    obs['request'] = ' '*80
    red=imred.Reducer('SONG')
    for req,o in enumerate(obs) :
        fig,ax=plots.multi(1,3,figsize=(8,4),hspace=0.001)
        for i,f in enumerate(o['files']) : 
            try :
                t,sn=reduce.plot(ax,f.decode(),os.path.basename(f.decode()),red=red,write=True)
                o['files'][i] = os.path.basename(o['files'][i])
                ind=np.where(np.char.find(out['file'],os.path.basename(f.decode()))>=0)[0]
                reduced=Table()
                reduced['exp_pk'] = [out[ind]['exp_pk'][0]]
                reduced['throughput'] = [t]
                reduced['sn'] = [sn]
                d=database.DBSession()
                d.ingest('obs.reduced',reduced,onconflict='update')
                d.close()
            except :
                print('error with ',f)
                pass
        for i in range(3) :
            lim = ax[i].get_ylim()
            ax[i].set_ylim(0,lim[1])
        ax[0].legend(fontsize='xx-small')
        fig.savefig('{:s}/{:s}/{:s}_{:d}.png'.format(root,ut,o['targname'],req))
        o['request'] = '<A HREF={:s}_{:d}.png> {:s} </A>'.format(o['targname'],req,o['targname'])
        if pause : pdb.set_trace()
        plt.close()

    for o in [out,guideout] :
        gd = np.where(o['file'] != '')[0]
        o=o[gd]
        j = np.argsort(o['dateobs'])
        o=o[j]
        o['exptime'].format='{:.2f}'

    fp = open('{:s}/{:s}/{:s}.html'.format(root,ut,ut),'w')

    fp.write('<HTML><BODY>\n')
    fp.write('<h2>{:s}  MJD: {:d}</h2><br>'.format(ut,mjd))
    fp.write('<A HREF=guide.html>Guider movies</A><BR>\n')
    fp.write('<A HREF=focus.html>Focus curves</A><BR>\n')
    fp.write('<A HREF=../reduced/{:s}>Reduced data</A><BR>\n'.format(ut))

    fp.write('<p>Observed robotic requests: <BR>\n')
    obs['files'] = obs['files'].astype('<U')
    tab = html.tab(obs['request','mjd','schedulename','sequencename','priority','files'],file=fp)

    fp.write('<p>Spectrograph exposures: <BR>\n')
    tab = html.tab(out['file','dateobs','ra','dec','exptime','camera','filter','focus','ccdtemp'],file=fp)

    fp.write('<p>Guider exposures: <BR>\n')
    tab = html.tab(guideout['file','dateobs','ra','dec','exptime','camera','filter','focus','ccdtemp'],file=fp)
    fp.write('</BODY></HTML>\n')
    fp.close()

    try : os.symlink('{:s}.html'.format(ut),'{:s}/{:s}/0_{:s}.html'.format(root,ut,ut))
    except: pass

    return out
