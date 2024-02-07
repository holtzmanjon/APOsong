from astropy.coordinates import SkyCoord,EarthLocation
from astropy.time import Time
from astropy.table import Table, hstack, vstack
import astropy.units as u

import pdb
import numpy as np

import aposong 
import database

class Request() :
    def __init__(self,target,schedule,sequence) :
        self.target = target    
        self.schedule = schedule
        self.sequence = sequence

    def table(self) :
        tab =Table()
        tab['name'] = [self.target.name]
        tab['ra'] = [self.target.ra]
        tab['dec'] = [self.target.dec]
        tab['epoch'] = [self.target.epoch]
        tab['min_airmass'] = [self.schedule.min_airmass]
        tab['max_airmass'] = [self.schedule.max_airmass]
        tab['nvisits'] = [self.schedule.nvisits]
        tab['dt_visit'] = [self.schedule.dt_visit]
        tab['nsequence'] = [self.schedule.nsequence]
        tab['sequence'] = [self.sequence.name]
        return tab

    def observe(self) :
        aposong.slew(self.target.ra,self.target.dec)
        aposong.expose(self.sequence.expose,self.target.name)


class Target() :
    def __init__(self,name,ra,dec,epoch=2000.) :
        self.name=name
        self.ra=ra
        self.dec=dec
        self.epoch=epoch

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

    def observe(self,name) :
        for filt,nexp,texp in zip(self.filt,self.n_exp,self.t_exp) :
            for iexp in range(nexp) : 
                aposong.expose(texp,filt,write=name)

    def table(self) :
        tab =Table()
        tab['sequencename' ] = [self.name]
        tab['filter'] = [self.filt]
        tab['n_exp'] = [self.n_exp]
        tab['t_exp'] = [self.t_exp]
        return tab

#class Observation() :
#    def __init() :

def hamax(dec,airmax,lat) :
    hamax=(1/airmax - np.sin(lat)*np.sin(dec) ) / np.cos(lat)*np.cos(dec)
    hamax=np.arccos(hamax)
    return hamax

def secz(ha,dec,lat) :
    secz = (np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.cos(ha))**-1
    return secz

def getrequests() :
    d=database.DBSession()
    tab=d.query(sql='SELECT targ.*,sched.*,seq.* FROM robotic.request as req  \
         JOIN robotic.target as targ on req.targname = targ.targname \
         JOIN robotic.schedule as sched on req.schedulename = sched.schedulename \
         JOIN robotic.sequence as seq on req.sequencename = seq.sequencename',fmt='table')
    d.close()
    return tab

def getbest(tab, t, site='APO', criterion='setting') :
    apo=EarthLocation.of_site(site)
    t.location=apo
    lst=t.sidereal_time('mean')
    print('lst: ', lst)
    if criterion == 'setting' : tmin=12*u.hourangle
    elif criterion == 'longest' : tmin=0.*u.hourangle
    best=-1
    am_best=-1
    dt_best=-1
    for iline,line in enumerate(tab) :
        c=SkyCoord("{:s} {:s}".format(line['ra'],line['dec']),unit=(u.hourangle,u.deg))
        ha = t.sidereal_time('mean') - c.ra
        ha.wrap_at(12*u.hourangle,inplace=True)
        am = secz(ha,c.dec,apo.lat)

        length = 5400.
        haend=ha+(length/3600.)*u.hourangle
        ham=hamax(c.dec,line['max_airmass'],apo.lat).to(u.hourangle)
        dt=ham-haend
        if am < line['min_airmass'] or am > line['max_airmass'] or dt<0 : continue
        if line['observed'] : continue
        if (criterion == 'setting' and dt< tmin) or (criterion=='longest' and dt>tmin) :
            best = line
            bestline = iline
            tmin=ham-haend
            am_best=am
            dt_best=ham-haend

    return best,am_best,dt_best,bestline

def test() :
    t = Time.now()
    while True :
        best,am,dt = getbest(tab,t)
        print(t,am,dt)
        print(best)
        tab.remove_row
        t += 5400*u.s

def maketarg(file,insert=False) :
    fp=open(file,'r')
    lines=fp.readlines()
    requests=[]
    for iline,line in enumerate(lines[2::2]) :
        name,ra,dec,epoch,amin,amax,y,m,d,m,d,umax,utmin,utmax,hmin,hmax,pmax,dmin,nv,dt,foc,df,g,gloc,skip=line.split()
        targ=Target(name,ra,dec)
        schedule=Schedule('rv',min_airmass=float(amin),max_airmass=float(amax),nvisits=int(nv),dt_visit=float(dt))
        sequence=Sequence('UBVRI',n_exp=[2,2,2,2,2],filt=['U','B','V','R','I'],t_exp=[300,60,60,60,120])
        request=Request(targ,schedule,sequence)
        requests.append(request.table())
        c=SkyCoord("{:s} {:s}".format(ra,dec),unit=(u.hourangle,u.deg))
        print(request.table())
        if iline == 0 :
            rtab = Table( [[name], ['rv'],['UBVRI']], names = ('targname','schedulename','sequencename') )
        else :
            rtab.add_row([name,'rv','UBVRI'])
    tab=vstack(requests)
    if insert :
        tab.rename_column('name','targname')
        d=database.DBSession()
        d.ingest('robotic.target',tab['targname','ra','dec','epoch'],verbose=True,onconflict='update')
        pdb.set_trace()
        d.ingest('robotic.schedule',schedule.table())
        d.ingest('robotic.sequence',sequence.table())
        d.ingest('robotic.request',rtab)
        d.close()
    return tab
