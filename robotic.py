from astropy.coordinates import SkyCoord,EarthLocation
from astropy.time import Time
from astropy.table import Table, hstack, vstack
import astropy.units as u
from astroplan import Observer, time_grid_from_range
import matplotlib
import matplotlib.pyplot as plt
from pyvista import imred,tv
from holtztools import html

import copy
import glob
import os
import pdb
import numpy as np
import time

import logging
import yaml
import logging.config
 
import focus as dofocus

try:
    with open('logging.yml', 'rt') as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)
except:
    print("can't open logging.yml")

logger=logging.getLogger(__name__)

import aposong 
import database
import APOSafety

class Target() :
    def __init__(self,name,ra,dec,epoch=2000.) :
        self.name=name
        self.ra=ra
        self.dec=dec
        self.epoch=epoch

    def table(self) :
        tab =Table()
        tab['targname'] = [self.name]
        tab['ra'] = [self.ra]
        tab['dec'] = [self.dec]
        tab['epoch'] = [self.epoch]

    def acquire(self,display=None,prop=0.0) :
        aposong.guide(False)
        aposong.slew(self.ra,self.dec)
        aposong.guide(True,exptime=0.5,name='guide/{:s}'.format(self.name),display=None,vmax=30000,prop=prop)

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
    def __init__(self,name,filt=['U','B','V','R','I'],n_exp=[1,1,1,1,1],t_exp=[1,1,1,1,1],camera=[0,0,0,0,0]) :
        self.name=name
        self.n_exp=n_exp
        self.t_exp=t_exp
        self.filt=filt
        self.camera=camera

    def length(self,acquisition=60,filtmove=5,read=15) :
        tot = acquisition
        for filt,nexp,texp in zip(self.filt,self.n_exp,self.t_exp) :
            tot+= filtmove+nexp*(texp+read)
        return tot

    def observe(self,name,display=None) :
        names = []
        for filt,nexp,texp,cam in zip(self.filt,self.n_exp,self.t_exp,self.camera) :
            for iexp in range(nexp) :
                logger.info('Expose camera: {:d} exptime: {:.2f}, '.format(cam,texp))
                exp=aposong.expose(texp,filt,name=name,display=display,cam=cam)
                names.append(exp.name)
        return names

    def table(self) :
        tab =Table()
        tab['sequencename' ] = [self.name]
        tab['filter'] = [self.filt]
        tab['n_exp'] = [self.n_exp]
        tab['t_exp'] = [self.t_exp]
        tab['camera'] = [self.camera]
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

def getbest(t=None, requests=None, site='APO', criterion='setting',mindec=-90,maxdec=90,skip=None) :
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
                       t_exp=request['t_exp'],camera=request['camera'])
        length=seq.length()
        haend=ha+(length/3600.)*u.hourangle
        ham=hamax(c.dec,request['max_airmass'],apo.lat).to(u.hourangle)
        dt=ham-haend
        hamid=ha+(length/3600./2.)*u.hourangle
        if am < request['min_airmass'] or am > request['max_airmass'] or dt<0 : 
            logger.info('{:s} out of airmass range {:.2f}'.format(request['targname'],am))
            continue
        #if c.dec<mindec*u.deg or c.dec>maxdec*u.deg : 
        #    logger.info('{:s} out of declination range'.format(request['targname']))
        #    continue
        obs=d.query(sql='SELECT * from robotic.observed WHERE request_pk = {:d}'.format(request['request_pk']))
        if len(obs)>0:
            if request['nvisits'] > 0 and len(obs) > request['nvisits'] : 
                logger.info('{:s} observed enough visits {:d}'.format(request['targname'],len(obs)))
                continue
            dt_visit = t.mjd - np.max(obs['mjd'])
            if dt_visit < request['dt_visit'] :
                logger.info('{:s} observed too recently'.format(request['targname']))
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
        logger.info('request selected: {:s}'.format(best['targname']))
    return best

def observe_object(request,display=None,acquire=True) :
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
        seq = Sequence(request['targname'],filt=request['filter'],
                       n_exp=request['n_exp'],t_exp=request['t_exp'],camera=request['camera'])
        t=Time.now()
        names=seq.observe(targ.name,display=display)
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

def obsopen(opentime,safety) :
    """ Open observatory at/after requested time and when safe 
    """
    while (Time.now()-opentime)<0 :
        # wait until sunset + dt_sunset hours 
        logger.info('waiting for sunset: {:.3f} '.format((opentime-Time.now()).to(u.hour).value,' hours'))
        time.sleep(60)

    while not safety.issafe() :
        # wait until safe to open based on Safety
        time.sleep(30)

    # open if not already open (e.g., from previous observe invocation same night)
    if aposong.D.ShutterStatus != 0 :
        aposong.domeopen()
    logger.info('open at: {:s}'.format(Time.now().to_string()))

def observe(foc0=28800,dt_focus=1.5,display=None,dt_sunset=0,dt_nautical=-0.4,obs='apo',tz='US/Mountain',criterion='best',maxdec=None) :
    """ Full observing night sequence 
    """
    site=Observer.at_site(obs,timezone=tz)
    sunset =site.sun_set_time(Time.now(),which='nearest')
    sunrise =site.sun_rise_time(Time.now(),which='next')
    nautical = site.twilight_evening_nautical(Time.now(),which='nearest')
    nautical_morn = site.twilight_morning_nautical(Time.now(),which='next')

    # setup up Safety object and open dome when safe after desired time relative to sunset
    safety = APOSafety.Safety()
    obsopen(sunset+dt_sunset*u.hour,safety)

    # wait for nautical twilight
    while (Time.now()-(nautical+dt_nautical*u.hour)).to(u.hour) < 0*u.hour :
        try :
            # if dome has been closed by 3.5m but can now be opened again, do it
            if safety.issafe() and aposong.D.ShutterStatus != 0 :
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

    # focus star on meridian 
    foc=focus(foc0=foc0,display=aposong.disp)
    foctime=Time.now()

    oldtarg=''
    skiptarg=None
    while (Time.now()-nautical_morn).to(u.hour) < 0*u.hour : 
      try :
        tnow=Time.now()
        logger.info('nautical twilight in : {:.3f}'.format((nautical_morn-tnow).to(u.hour).value))
        if not safety.issafe() : 
            aposong.domeclose()
            time.sleep(90)
            continue
        elif aposong.D.ShutterStatus != 0 :
            aposong.domeopen()

        logger.info('tnow-foctime : {:.3f}'.format((tnow-foctime).to(u.hour).value))
        if (tnow-foctime).to(u.hour) > dt_focus*u.hour :
            aposong.guide(False)
            foc0=focus(foc0=foc0,display=display)
            foctime=tnow
        else :
            best=getbest(criterion=criterion,maxdec=maxdec,skip=skiptarg)
            if best is None :
                time.sleep(60)
            else :
                success = observe_object(best,display=display,acquire=(best['targname']!=oldtarg))
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

def focus(foc0=28800,display=None) :
    """ Do focus run for object on meridian
    """    
    t=Time(Time.now(), location=EarthLocation.of_site('APO'))
    lst=t.sidereal_time('mean').value
    aposong.usno(ra=lst,dec=10.,magmin=9,magmax=10,rad=5*u.degree)

    f=aposong.focrun(foc0,75,9,2,None,bin=1,thresh=100,display=display,max=5000)
    while f<foc0-150 or f>foc0+150 :
        logger.info('focus: {:d}  foc0: {:d}'.format(f,foc0))
        foc0=copy.copy(f)
        f=aposong.focrun(foc0,75,9,2,None,bin=1,thresh=100,display=display,max=5000)
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
        if iline == 0 :
            rtab = Table( [[name], [schedule],[sequence]], names = ('targname','schedulename','sequencename') )
        else :
            rtab.add_row([name,'rv','UBVRI'])
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

def loadseq(name,t_exp=[1],n_exp=[1],filt=['V'],camera=[0]) :
    d=database.DBSession()
    sequence=Sequence(name,filt=filt,t_exp=t_exp,n_exp=n_exp,camera=camera)
    d.ingest('robotic.sequence',sequence.table(),onconflict='update')
    d.close()

def mkhtml(t=None) :
    if t is None : t = Time.now()
    y,m,d,hr,mi,se = t.ymdhms
    ut = 'UT{:d}{:02d}{:02d}'.format(y-2000,m,d)
    mkmovie(ut)
    mkplots(int(t.mjd))


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
        red.movie(range(seq[0],seq[1]),display=t,max=10000,out=out)
        if i>0 and i%5 == 0 :
            grid.append(row)
            row=[]
        row.append('guide/{:d}.gif'.format(seq[0]))
    for ii in range(i%5+1,5) :
        row.append('')
    grid.append(row)
    html.htmltab(grid,file=root+ut+'/guide.html',size=200)

def mkplots(mjd,display=None,root='/data/1m/') :

    matplotlib.use('Agg')
    d=database.DBSession()
    out=d.query('obs.focus',fmt='list')
    i =np.where(np.array(out[0]) == 'mjd')[0][0]
    ifile =np.where(np.array(out[0]) == 'files')[0][0]
    mjds=[]
    files=[]
    for o in out[1:] : 
        mjds.append(o[i])
        files.append(o[ifile])

    gd=np.where(np.array(mjds).astype(int) == mjd)[0]

    grid=[]
    row=[]
    for i,seq in enumerate(gd) :
        print(files[seq])
        try : 
            dofocus.focus(files[seq],pixscale=0.22,display=display,root=root,plot=True,hard=True)
            if i>0 and i%5 == 0 :
                grid.append(row)
                row=[]
            row.append(os.path.basename(files[seq][0].replace('.fits','.png')))
        except : continue
    for ii in range(i%5+1,5) :
        row.append('')
    grid.append(row)
    dir=files[seq][0].split('/')[0]
    html.htmltab(grid,file=root+dir+'/focus.html',size=250)

