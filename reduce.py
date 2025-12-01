import os
import pdb
import matplotlib
import yaml
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

import numpy as np
from pyvista import imred, spectra, image, utils, centroid
from pyvista.dataclass import Data
import database
from holtztools import plots
from astropy.table import Table
from astropy.time import Time
import robotic
import aposong

def specreduce(n, red=None, trace=None, wav=None, retrace=False, cr=True, scat=False, response=None, 
               write=False, display=None, clobber=False, twod=False) :
    """ Quick reduction

    Parameters
    ----------
    n : str or int
        number or name of file to reduce
    red : Reducer object, default=None
        pyvista imred.Reducer object for doing basic reduction
    trace : Trace object, default=None
        pyvista spectra.Trace object for extraction
    wav : WaveCal object, default=None
        pyvista spectra.WaveCal object for wavelength calibration
    retrace : bool, default=False
        set to True to retrace on object
    cr : bool, default=True
        set to True to do CR cleaning
    scat : bool, default=False
        set to True to do scattered light correction
    response :  pyvista Data object, default=None
        response function to load into output
    write : bool, default=False
        set to True to write extracted frame to disk
    display : pyvista TV object, default=None
        if given, display data during reduction, and pause for inspection
    clobber : bool, default=False
        if True, re-reduce objects even if output file exists already
    """

    try :
        with open('aposong.yml','r') as config_file :
            config = yaml.safe_load(config_file)
    except:
        print('no configuration file found')
    dataroot = config['dataroot']+'/'

    if red == None :
        red=imred.Reducer('SONG',dir=dataroot+'UT251003')
    if response == None :
        response=Data.read(dataroot+'cal/response.fits')
    if cr : crbox=[1,11]
    else : crbox=None
    if scat : doscat=red.scat
    else : doscat=None

    # get header for filename and exptime for dark (dark subtraction will still scale)
    im=red.rd(n)
    file=im.header['FILE'].split('.')
    if isinstance(n,str) :
        outfile='{:s}/{:s}_ec.{:s}.fits'.format(os.path.dirname(n).replace('1m/','1m/reduced/'),file[0],file[-2])
    else :
        outfile='{:s}/{:s}_ec.{:s}.fits'.format(red.dir.replace('1m/','1m/reduced/'),file[0],file[-2])

    # if alreadly done, read and return, unless clobber
    if os.path.exists(outfile) and not clobber :
        return Data.read(outfile)

    # set appropriate cals based on mjd
    if Time(im.header['DATE-OBS']).mjd < 60997 :
        darktimes=np.array([25,60,120,180,240,300])
        utdark='UT251010'
        utflat='UT251119'
        if trace == None : trace=spectra.Trace(dataroot+'cal/trace/UT251007_Trace_fiber2.fits')
        if wav == None : wav=spectra.WaveCal(dataroot+'cal/wavecal/UT251007_WaveCal_fiber2.fits')
    else :
        darktimes=np.array([30,60,120,180,240,300,600])
        utdark='UT251119'
        utflat='UT251119'
        if trace == None : trace=spectra.Trace(dataroot+'cal/trace/UT251119_Trace_fiber2.fits')
        if wav == None : wav=spectra.WaveCal(dataroot+'cal/wavecal/UT251119_WaveCal_fiber2.fits')

    # read dark and flat frames
    dtime=darktimes[np.argmin(abs(im.header['EXPTIME']-darktimes))]
    dark=Data.read(dataroot+'cal/darks/dark_{:d}_-10_{:s}.fits'.format(dtime,utdark))
    flat=Data.read(dataroot+'cal/pixflats/pixflat_flat_{:s}.fits'.format(utflat))

    # Reduce
    im=red.reduce(n,crbox=crbox,scat=doscat,dark=dark,flat=flat,display=display)
    if twod : return im

    # Extract
    if retrace : trace.retrace(im,display=display)
    traceshift=trace.find(im,lags=range(-10,11))
    imec=trace.extract(im,display=display)
    imec.header['T_SHIFT'] = traceshift[0]

    # Add wave and response, barycentric
    wav.add_wave(imec)
    imec.add_response(response.data)
    bv, bt = utils.getbc(imec.header)
    imec.header['BVC'] = bv.value
    imec.header['BJD-MID'] = bt.jd

    # write if requested
    if write :
        try : os.makedirs(os.path.dirname(outfile))
        except : pass
        imec.write(outfile)
    return imec

def throughput(spec,ax,name,mag=None,song=None,red=None,orders=[34]) :
    """ Reduce and add plot of extracted spectrum, throughput, and S/N for an image to existing Axes

    Parameters
    ==========
    ax : Axes[3]
        matplotlib Axes to plot into
    i : int or str
        number or name of file to reduce
    name : str
        name for plot legend
    mag : float, default=None
        magnitude of object to compute throughput, if not given, query database
    orders : list, default=[34]
        Order(s) to plot
    """
    #spec=specreduce(i,red=red,write=write,clobber=clobber)
    targ=spec.header['OBJECT'].split('.')[0].replace('_iodine','')
    if mag == None :
        d=database.DBSession() 
        targs=d.query('robotic.target')
        j=np.where(targs['targname'] == targ)[0]
        mag=targs['mag'][j][0]
        d.close()
    if targ == 'muHer' : song = 25000/120
    elif targ == 'HD185144' : song = 17500/300
    elif targ == 'gammaCep' or targ == 'gamCep' : song = 35000/120
    for order in orders :
        ax[0].plot(spec.wave[order],spec.data[order]*red.gain/spec.header['EXPTIME']*10**(0.4*mag),label=name)
        ax[1].plot(spec.wave[order],spec.data[order]*red.gain,label=name)
        ax[2].plot(spec.wave[order],spec.data[order]/spec.uncertainty.array[34],label=name)
    if song != None : ax[0].scatter([5564],[song*10**(0.4*mag)],s=50)
    ax[0].legend(fontsize='xx-small')
    ax[0].set_ylabel('photons/s for m=0')
    ax[1].set_ylabel('photons')
    ax[2].set_ylabel('S/N')
    ax[2].set_xlabel('wavelength')
    j1=np.argmin(np.abs(spec.wave-5570))
    j2=np.argmin(np.abs(spec.wave-5560))
    t=np.max(spec.data.flatten()[j1:j2])*red.gain[0]/spec.header['EXPTIME']*10**(0.4*mag)
    sn=np.max(spec.data.flatten()[j1:j2]/spec.uncertainty.array.flatten()[j1:j2])
    t_orders=np.nanmedian(spec.data[:,1400:3400],axis=1)*red.gain[0]/spec.header['EXPTIME']*10**(0.4*mag)
    sn_orders=np.nanmedian(spec.data[:,1400:3400]/spec.uncertainty.array[:,1400:3400],axis=1)

    return t,sn,t_orders,sn_orders

def throughput_all() :
    """ Make plot of throughput from database query
    """

    pdb.set_trace()
    d=database.DBSession()
    outlist=d.query(sql="select * from obs.reduced as red join obs.exposure as exp on red.exp_pk = exp.exp_pk")
    out=Table(dtype=outlist.dtype)
    outlist=d.query(sql="select * from obs.reduced as red join obs.exposure as exp on red.exp_pk = exp.exp_pk",fmt="list")
    for o in outlist[1:] :
        if len(o[4]) == 56 : o[4].append(0.)
        if len(o[5]) == 56 : o[5].append(0.)
        out.add_row(o)

    targs=d.query('robotic.target')
    d.close()
    i=np.where(out['filter'] == 'iodine')
    o=np.where(out['filter'] != 'iodine')
    mag=[]
    for f in out['file'] :
        targ = os.path.basename(f).split('.')[0]
        if targ != '':
            j=np.where(targs['targname'] == targ)[0][0]
            mag.append(targs[j]['mag'])
    mag=np.array(mag)
    fig,ax=plots.multi(1,2,figsize=(10,8),hspace=0.001,sharex=True)
    ax[0].scatter(out['mjd'][i],out['throughput'][i],c=out['alt'][i],marker='+',label='iodine')
    scat=ax[0].scatter(out['mjd'][o],out['throughput'][o],c=out['alt'][o],marker='s',label='no iodine')
    plt.colorbar(scat)
    ax[1].scatter(out['mjd'][i],out['throughput'][i],c=mag[i],marker='+',label='iodine')
    scat=ax[1].scatter(out['mjd'][o],out['throughput'][o],c=mag[o],marker='s',label='no iodine')
    plt.colorbar(scat)
    xlim=ax[0].get_xlim()
    for targ in ['muHer','HD185144','gammaCep'] :
        if targ == 'muHer' :
            song = 25000/120
            mag = 3.42
        elif targ == 'HD185144' :
            song = 17500/300
            mag=4.67
        elif targ == 'gammaCep' or targ == 'gamCep' :
            mag=3.22
            song = 35000/120
        ax[0].plot(xlim,[song*10**(0.4*mag),song*10**(0.4*mag)],label='Tenerife {:s}'.format(targ))

    ax[0].set_ylim(0,6000)
    ax[1].set_ylim(0,6000)
    ax[0].legend()
    ax[1].set_xlabel('MJD')
    ax[0].set_ylabel('cnts/s scaled to m=0')
    ax[0].set_title('Throughput at 5560-70 (color coded by altitude)')
    ax[1].set_title('Throughput at 5560-70 (color coded by mag)')
    fig.tight_layout()

def guider(i1,i2,red=None,sat=65000,title='') :
    """ Make plots of derived star positions from guider sequence
    """
    fig,ax=plots.multi(2,2)
    for i in range(i1,i2) :
        a=red.rd(i)
        c=centroid.maxtot(a,90,93,rad=6)
        s=centroid.rasym_centroid(a,c.x,c.y,sat=sat)
        if i==i1 : label='maxtot position'
        else : label=None
        ax[0,0].scatter(i,c.x,c='r',label=label)
        ax[1,0].scatter(i,c.y,c='r',label=label)
        ax[0,1].scatter(c.x,c.y,c='r',label=label)
        if i==i1 : label='asym position'
        else : label=None
        ax[0,0].scatter(i,s.x,c='b',label=label)
        ax[1,0].scatter(i,s.y,c='b',label=label)
        ax[0,1].scatter(s.x,s.y,c='b',label=label)

    ax[0,0].set_xlabel('image number')
    ax[0,0].set_ylabel('x position')
    ax[0,0].legend()
    ax[1,0].set_xlabel('image number')
    ax[1,0].set_ylabel('y position')
    ax[1,0].legend()
    ax[0,1].set_xlim(81.5,91.5)
    ax[0,1].set_ylim(81.5,91.5)
    ax[0,1].set_xlabel('x position')
    ax[0,1].set_ylabel('y position')
    ax[0,1].legend()
    ax[0,1].grid()
    ax[0,1].scatter([86.5],[86.5],s=200,c='g',marker='+')
    ax[1,1].imshow(a.data,vmin=0,vmax=sat)

    fig.suptitle('Guider positions, sat: {:d} {:s}'.format(sat,title))
    fig.tight_layout()

def getobs(targ) :

    d=database.DBSession()
    out=d.query(sql='select * from obs.reduced as red join obs.exposure as exp on red.exp_pk = exp.exp_pk')
    obslist=d.query(sql='select * from robotic.observed as obs join robotic.request as req on obs.request_pk = req.request_pk',fmt='list')
    d.close()

    maxfiles=0
    for o in obslist[1:] :
        maxfiles=np.max([maxfiles,len(o[3])])

    obs=Table(names=('request_pk', 'mjd', 'targname', 'schedulename','sequencename','priority','files'), dtype=('i4','f4', 'S', 'S','S','i4','{:d}S24'.format(maxfiles)))
    for o in obslist[1:] :
        while len(o[3]) < maxfiles : o[3].append('')
        try: 
            for i,oo in enumerate(o[3]) :
                if oo == None : o[3][i] = ''
            obs.add_row([o[1],o[2],o[5],o[6],o[7],o[8],o[3]])
        except: 
            pdb.set_trace()

    j=np.where(obs['targname']==targ)
    ind=[]
    for file in obs[j]['files'] :
        for f in file :
            i=np.where(out['file']  == f.decode().replace('/data/1m/',''))[0]
            if len(i) > 0 :
                ind.append(i[0])
    ind=np.array(ind)
    sort =np.argsort(out[ind]['sn'])[::-1]
    tab=Table(out[ind[sort]])['file','filter','throughput','sn']
    tab['sn'].format = '.2f'
    tab['throughput'].format = '.2f'
    return tab


def reduceall(mjdmin) :

    matplotlib.use('Agg')
    d=database.DBSession()
    obs=d.query(sql='select * from robotic.observed as obs join robotic.request as req on obs.request_pk = req.request_pk',fmt='list')
    d.close()

    mjds=[]
    for o in obs[1:] :
        mjds.append(o[2])
    mjds=sorted(set(np.array(mjds).astype(int)))

    for mjd in mjds :
        if mjd > mjdmin :
            robotic.mklog(mjd,clobber=True)
            #robotic.mkhtml(mjd)




