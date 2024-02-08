from astropy.coordinates import SkyCoord,EarthLocation
from astropy.time import Time
from astropy.table import Table, hstack, vstack
import astropy.units as u

import pdb
import numpy as np

import aposong 
import database

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

    def acquire(self) :
        aposong.slew(self.ra,self.dec)

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
    def __init__(self,name,filt=['U','B','V','R','I'],n_exp=[1,1,1,1,1],t_exp=[1,1,1,1,1]) :
        self.name=name
        self.n_exp=n_exp
        self.t_exp=t_exp
        self.filt=filt

    def length(acquisition=60,filtmove=5,read=15) :
        tot = acquisition
        for filt,nexp,texp in zip(self.filt,self.n_exp,self.t_exp) :
            tot+= filtmove+nexp*(texp+read)
        return length

    def observe(self,name,display=None) :
        names = []
        for filt,nexp,texp in zip(self.filt,self.n_exp,self.t_exp) :
            for iexp in range(nexp) : 
                hdu,outname=aposong.expose(texp,filt,name=name,display=display)
                names.append(outname)
        return names

    def table(self) :
        tab =Table()
        tab['sequencename' ] = [self.name]
        tab['filter'] = [self.filt]
        tab['n_exp'] = [self.n_exp]
        tab['t_exp'] = [self.t_exp]
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
    tab=d.query(sql='SELECT request_pk,targ.*,sched.*,seq.* FROM robotic.request as req  \
         JOIN robotic.target as targ on req.targname = targ.targname \
         JOIN robotic.schedule as sched on req.schedulename = sched.schedulename \
         JOIN robotic.sequence as seq on req.sequencename = seq.sequencename' ,fmt='table')
    d.close()
    return tab

def getbest(t=None, requests=None, site='APO', criterion='setting') :
    """ 
    Get best request given table of requests and time
    """
    if t is None :
        t = Time.now()
    if requests is None :
        requests = getrequests()

    apo=EarthLocation.of_site(site)
    t.location=apo
    lst=t.sidereal_time('mean')
    print('lst: ', lst)
    if criterion == 'setting' : tmin=12*u.hourangle
    elif criterion == 'longest' : tmin=0.*u.hourangle
    else : 
        print('unknown criterion: ', criterion)
        return
    best=-1
    am_best=-1
    dt_best=-1
    d=database.DBSession()
    for iline,line in enumerate(requests) :
        if line['targname'] is None : continue
        c=SkyCoord("{:s} {:s}".format(line['ra'],line['dec']),unit=(u.hourangle,u.deg))
        ha = t.sidereal_time('mean') - c.ra
        ha.wrap_at(12*u.hourangle,inplace=True)
        am = secz(ha,c.dec,apo.lat)

        length = 5400.
        haend=ha+(length/3600.)*u.hourangle
        ham=hamax(c.dec,line['max_airmass'],apo.lat).to(u.hourangle)
        dt=ham-haend
        if am < line['min_airmass'] or am > line['max_airmass'] or dt<0 : continue
        obs=d.query(sql='SELECT * from robotic.observed WHERE request_pk = {:d}'.format(line['request_pk']))
        if len(obs)>0:
            print(obs)
            if line['nvisits'] > 0 and len(obs) > line['nvisits'] : 
                print(line['targname'], 'observed enough visits',len(obs))
                continue
            dt_visit = t.mjd - np.max(obs['mjd'])
            if dt_visit < line['dt_visit'] :
                print(line['targname'], 'observed too recently')
                continue

        if (criterion == 'setting' and dt< tmin) or (criterion=='longest' and dt>tmin) :
            best = line
            bestline = iline
            tmin=ham-haend
            am_best=am
            dt_best=ham-haend
    d.close()
    print(best)
    return best

def observe(request,display=None) :
    """ Given request, do the observation and record
    """
    targ = Target(request['targname'],request['ra'],request['dec'],epoch=request['epoch'])
    targ.acquire()
    t=Time.now()
    seq = Sequence(request['targname'],filt=request['filter'],n_exp=request['n_exp'],t_exp=request['t_exp'])
    names=seq.observe(targ.name,display=display)

    obs=Table()
    obs['request_pk'] = [request['request_pk']]
    obs['mjd'] = [t.mjd]
    obs['files'] = [names]
    d=database.DBSession()
    d.ingest('robotic.observed',obs,onconflict='update')
    d.close()
    return obs

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

def loadseq(name,t_exp=[1],n_exp=[1],filt=['V']) :
    d=database.DBSession()
    sequence=Sequence(name,filt=filt,t_exp=t_exp,n_exp=n_exp)
    d.ingest('robotic.sequence',sequence.table(),onconflict='update')
    d.close()
