import os
import pdb
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.ion()

import numpy as np
from pyvista import imred, spectra, image, utils, centroid
from pyvista.dataclass import Data
import database
from holtztools import plots

def specreduce(n, red=None, trace=None, wav=None, retrace=False, cr=True, scat=False, response=None, write=False, display=None, clobber=False) :
    """ Quick reduction
    """

    if red == None :
        red=imred.Reducer('SONG',dir='/data/1m/UT251003')
    if trace == None :
        trace=spectra.Trace('SONG/UT251007_Trace_fiber2.fits')
    if wav == None :
        wav=spectra.WaveCal('SONG/UT251007_WaveCal_fiber2.fits')
    if response == None :
        response=Data.read('/data/1m/cal/response.fits')
    if cr : crbox=[1,11]
    else : crbox=None
    if scat : doscat=red.scat
    else : doscat=None

    # get closest dark for exptime (dark subtraction will still scale)
    darktimes=np.array([25,60,120,180,240,300])
    im=red.rd(n)
    file=im.header['FILE'].split('.')
    if isinstance(n,str) :
        out='{:s}/{:s}_ec.{:s}.fits'.format(os.path.dirname(n).replace('1m/','1m/reduced/'),file[0],file[-2])
    else :
        out='{:s}{:s}_ec.{:s}.fits'.format(red.dir.replace('1m/','1m/reduced'),file[0],file[-2])
    if os.path.exists(out) and not clobber :
        return Data.read(out)
    dtime=darktimes[np.argmin(abs(im.header['EXPTIME']-darktimes))]
    dark=Data.read('/data/1m/cal/darks/dark_{:d}_-10_UT251010.fits'.format(dtime))
    im=red.reduce(n,crbox=crbox,scat=doscat,dark=dark,display=display)
    if retrace : trace.retrace(im,display=display)
    imec=trace.extract(im,display=display)
    wav.add_wave(imec)
    imec.add_response(response.data)
    bv, bt = utils.getbc(imec.header)
    imec.header['BVC'] = bv.value
    imec.header['BJD-MID'] = bt.jd
    if write :
        try : os.makedirs(os.path.dirname(out))
        except : pass
        imec.write(out)
    return imec

def plot(ax,i,name,mag=None,song=None,red=None,orders=[34],write=False):
    spec=specreduce(i,red=red,write=write)
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
    return t,sn

def mjd(ut='UT251003',rn=None, write=False) :
    red=imred.Reducer('SONG',dir='/data/1m/'+ut)
    if rn is not None : red.rn = [rn]
    fig,ax=plots.multi(1,3,figsize=(8,4),hspace=0.001)
    if ut == 'UT251003' :
        plot(ax,61,'muHer_a',song=25000/120,red=red,write=write)
        plot(ax,62,'muHer_b',red=red,write=write)
        plot(ax,63,'muHer_i_a',red=red,write=write)
        plot(ax,64,'muHer_i_b',red=red,write=write)
        plot(ax,74,'HD185144_a',red=red,song=17500/300,write=write)
        plot(ax,75,'HD185144_b',red=red,write=write)
        plot(ax,76,'HD185144_i_a',red=red,write=write)
        plot(ax,77,'HD185144_i_b',red=red,write=write)
        plot(ax,78,'gammaCep_a',red=red,song=35000/120,write=write)
        plot(ax,79,'gammaCep_b',red=red,write=write)
        plot(ax,80,'gammaCep_i_a',red=red,write=write)
        plot(ax,81,'gammaCep_i_b',red=red,write=write)
    elif ut == 'UT251006' :
        plot(ax,113,'gammaCep_a',red=red,write=write)
        plot(ax,114,'gammaCep_b',red=red,write=write)
        plot(ax,115,'gammaCep_i_a',red=red,write=write)
        plot(ax,116,'gammaCep_i_b',red=red,write=write)
        plot(ax,117,'HD185144_a',red=red,write=write)
        plot(ax,118,'HD185144_b',red=red,write=write)
    elif ut == 'UT251011' :
        plot(ax,64,'HR7235_a',red=red,write=write)
        plot(ax,65,'HR7235_b',red=red,write=write)
        plot(ax,66,'HR7235_i_a',red=red,write=write)
        plot(ax,67,'HR7235_i_b',red=red,write=write)
        plot(ax,83,'HR6410_a',red=red,write=write)
        plot(ax,84,'HR6410_b',red=red,write=write)
        plot(ax,85,'HR6410_i_a',red=red,write=write)
        plot(ax,86,'HR6410_i_b',red=red,write=write)
        plot(ax,86,'HR6410_i_b',red=red,write=write)
        plot(ax,117,'HR264_a',red=red,write=write)
        plot(ax,118,'HR264_b',red=red,write=write)
        plot(ax,119,'HR264_i_a',red=red,write=write)
        plot(ax,120,'HR264_i_b',red=red,write=write)
        plot(ax,232,'HD10780_a',red=red,write=write)
        plot(ax,233,'HD10780_b',red=red,write=write)
        plot(ax,234,'HD10780_i_a',red=red,write=write)
        plot(ax,235,'HD10780_i_b',red=red,write=write)
        plot(ax,266,'gammaCep_a',red=red,write=write)
        plot(ax,267,'gammaCep_b',red=red,write=write)
        plot(ax,268,'gammaCep_i_a',red=red,write=write)
        plot(ax,269,'gammaCep_i_b',red=red,write=write)
    elif ut == 'UT251016' :
        plot(ax,224,'HR2845_a',red=red,write=write)
        plot(ax,225,'HR2845_b',red=red,write=write)
        plot(ax,226,'HR2845_i_a',red=red,write=write)
        plot(ax,227,'HR2845_i_b',red=red,write=write)
    elif ut == 'UT251017' :
        for targ,num in zip(['HR7235','HR6410','HR264'],[46,59,72]) :
            for i in range(4) :
              if i== 0 : ext = '_a'
              elif i== 1 : ext = '_b'
              elif i== 2 : ext = '_i_a'
              elif i== 3 : ext = '_i_b'
              plot(ax,num+i,targ+ext,red=red,write=write)
    elif ut == 'UT251025' :
        for targ,num in zip(['HD185144','HR264','HR1087','HR2845','HR3494'],[89,126,148,209,222]) :
            for i in range(4) :
              if i== 0 : ext = '_a'
              elif i== 1 : ext = '_b'
              elif i== 2 : ext = '_i_a'
              elif i== 3 : ext = '_i_b'
              plot(ax,num+i,targ+ext,red=red,write=write)
        for targ,num in zip(['HD23249','HD23249','HD23249','HD23249'],[161,173,185,197]) :
            for i in range(3) :
              if i== 0 : ext = '_a'
              elif i== 1 : ext = '_b'
              elif i== 2 : ext = '_c'
              plot(ax,num+i,targ+ext,red=red,write=write)

    d=database.DBSession() 
    targs=d.query('robotic.target')
    d.close()
    for targ in ['muHer','HD185144','gamCep'] :
        j=np.where(targs['targname'] == targ)[0]
        mag=targs['mag'][j][0]
        if targ == 'muHer' : song = 25000/120
        elif targ == 'HD185144' : song = 17500/300
        elif targ == 'gammaCep' or targ == 'gamCep' : song = 35000/120
        ax[0].scatter([5564],[song*10**(0.4*mag)],s=50)
    for i in range(3) :
        lim = ax[i].get_ylim()
        ax[i].set_ylim(0,lim[1])
    ax[0].legend(fontsize='xx-small',ncol=3)


    fig.suptitle(ut)
    fig.savefig(ut+'.png')
    #plt.close()
    return fig,ax

def throughput() :

    d=database.DBSession()
    out=d.query(sql='select * from obs.reduced as red join obs.exposure as exp on red.exp_pk = exp.exp_pk')
    d.close()
    i=np.where(out['filter'] == 'iodine')
    o=np.where(out['filter'] != 'iodine')
    plt.figure(figsize=(8,4))
    plt.scatter(out['mjd'][i],out['throughput'][i],c=out['alt'][i],marker='+',label='no iodine')
    plt.scatter(out['mjd'][o],out['throughput'][o],c=out['alt'][o],marker='s',label='iodine')
    plt.ylim(0,6000)
    plt.colorbar()
    xlim=plt.xlim()
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
        plt.plot(xlim,[song*10**(0.4*mag),song*10**(0.4*mag)],label='Tenerife {:s}'.format(targ))

    plt.legend()
    plt.xlabel('MJD')
    plt.ylabel('cnts/s scaled to m=0')
    plt.title('Throughput (color coded by altitude)')

def guider(i1,i2,red=None,sat=65000,title='') :

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

