import tkthread

import os
import pdb
from tkinter import *
from tkinter import ttk
import tkinter.font
import numpy as np
import yaml
from numpy.random import randint

from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
import astropy.units as u
import pwi4_client

import weather
import influx
import aposong
import eshel

class TelescopeWgt(ttk.Frame) :

    def __init__(self,container) :
        super().__init__(container)
        row=1
        ttk.Label(self,text="RA").grid(column=1,row=row,sticky=(W))
        self.ra = StringVar()
        ttk.Label(self, textvariable=self.ra).grid(column=2,row=row,sticky=(E),padx=15)

        ttk.Label(self,text="AZ").grid(column=3,row=row,sticky=(W))
        self.az = StringVar()
        ttk.Label(self, textvariable=self.az).grid(column=4,row=row,sticky=(E),padx=15) 

        ttk.Label(self,text="UT").grid(column=5,row=row,sticky=(W))
        self.ut = StringVar()
        ttk.Label(self, textvariable=self.ut).grid(column=6,row=row,sticky=(E),padx=10)


        row+=1
        ttk.Label(self,text="DEC").grid(column=1,row=row,sticky=(W))
        self.dec = StringVar()
        ttk.Label(self, textvariable=self.dec).grid(column=2,row=row,sticky=(E),padx=15)

        ttk.Label(self,text="ALT").grid(column=3,row=row,sticky=(W))
        self.alt = StringVar()
        ttk.Label(self, textvariable=self.alt).grid(column=4,row=row,sticky=(E),padx=15)

        ttk.Label(self,text="LST").grid(column=5,row=row,sticky=(W))
        self.lst = StringVar()
        ttk.Label(self, textvariable=self.lst).grid(column=6,row=row,sticky=(E),padx=10)

        row+=1
        ttk.Label(self,text="PA").grid(column=1,row=row,sticky=(W))
        self.pa = StringVar()
        ttk.Label(self, textvariable=self.pa).grid(column=2,row=row,sticky=(E),padx=15)

        ttk.Label(self,text="ROT").grid(column=3,row=row,sticky=(W))
        self.rot = StringVar()
        ttk.Label(self, textvariable=self.rot).grid(column=4,row=row,sticky=(E),padx=15)

        ttk.Label(self,text="MJD").grid(column=5,row=row,sticky=(W))
        self.mjd = StringVar()
        ttk.Label(self, textvariable=self.mjd).grid(column=6,row=row,sticky=(E),padx=10)

        row+=1
        ttk.Label(self,text="HA").grid(column=1,row=row,sticky=(W))
        self.ha = StringVar()
        ttk.Label(self, textvariable=self.ha).grid(column=2,row=row,sticky=(E),padx=15) 

        ttk.Label(self,text="FOCUS").grid(column=3,row=row,sticky=(W))
        self.focus = StringVar()
        ttk.Label(self, textvariable=self.focus).grid(column=4,row=row,sticky=(E),padx=15) 

        ttk.Label(self,text="PORT").grid(column=5,row=row,sticky=(W))
        self.port = StringVar()
        ttk.Label(self, textvariable=self.port).grid(column=6,row=row,sticky=(E),padx=15) 


class DomeWgt(ttk.Frame) :

    def __init__(self,container) :
        super().__init__(container) 

        row=1
        ttk.Label(self,text="SHUTTER").grid(column=1,row=row,sticky=(W))
        self.shutter = StringVar()
        ttk.Label(self, textvariable=self.shutter).grid(column=2,row=row,sticky=(E),padx=10)

        ttk.Label(self,text="DOME AZ").grid(column=3,row=row,sticky=(W))
        self.az = StringVar()
        ttk.Label(self, textvariable=self.az).grid(column=4,row=row,sticky=(E),padx=10)

        self.slewing = StringVar()
        ttk.Label(self, textvariable=self.slewing).grid(column=5,row=1,sticky=(W),padx=10)

        row+=1
        ttk.Label(self,text="MIRROR COVERS").grid(column=1,row=row,sticky=(W))
        self.coverstate = StringVar()
        ttk.Label(self, textvariable=self.coverstate).grid(column=2,row=row,sticky=(E),padx=10) 

        ttk.Label(self,text="35M/25M").grid(column=3,row=row,sticky=(W))
        self.stat35m = StringVar()
        self.statcolor = StringVar()
        ttk.Label(self, textvariable=self.stat35m, foreground=self.statcolor.get()).grid(column=4,row=row,sticky=(E),padx=10) 

     
class CameraWgt(ttk.Frame) :

    def __init__(self,container) :
        super().__init__(container) 

        ttk.Label(self,text="FILTER",width=8).grid(column=1,row=1,sticky=(W))
        self.filter = StringVar()
        ttk.Label(self, textvariable=self.filter).grid(column=2,row=1,sticky=(W,E),padx=10)

        ttk.Label(self,text="BINNING",width=8).grid(column=3,row=1,sticky=(W))
        self.binning = StringVar()
        ttk.Label(self, textvariable=self.binning).grid(column=4,row=1,sticky=(W,E),padx=10)

        ttk.Label(self,text="STATE",width=8).grid(column=5,row=1,sticky=(W))
        self.state = StringVar()
        ttk.Label(self, textvariable=self.state).grid(column=6,row=1,sticky=(W,E),padx=10)

        ttk.Label(self,text="TEMP",width=8).grid(column=1,row=2,sticky=(W))
        self.temperature = StringVar()
        ttk.Label(self, textvariable=self.temperature).grid(column=2,row=2,sticky=(W,E),padx=10)
        ttk.Label(self,text="COOLER POWER",width=8).grid(column=3,row=2,sticky=(W))
        self.cooler = StringVar()
        ttk.Label(self, textvariable=self.cooler).grid(column=4,row=2,sticky=(W,E),padx=10)

class IodineWgt(ttk.Frame) :

    def __init__(self,container) :
        super().__init__(container) 

        ttk.Label(self,text="IODINE STAGE",width=16).grid(column=1,row=1,sticky=(W))
        self.position = StringVar()
        ttk.Label(self, textvariable=self.position).grid(column=2,row=1,sticky=(W,E),padx=10)

        ttk.Label(self,text="IODINE TEMP",width=16).grid(column=3,row=1,sticky=(W))
        self.temp = StringVar()
        ttk.Label(self, textvariable=self.temp).grid(column=4,row=1,sticky=(W,E),padx=10)

        ttk.Label(self,text="IODINE VOLTAGE",width=16).grid(column=1,row=2,sticky=(W))
        self.voltage = StringVar()
        ttk.Label(self, textvariable=self.voltage).grid(column=2,row=2,sticky=(W,E),padx=10)

        ttk.Label(self,text="IODINE CURRENT",width=16).grid(column=3,row=2,sticky=(W))
        self.current = StringVar()
        ttk.Label(self, textvariable=self.current).grid(column=4,row=2,sticky=(W,E),padx=10)

class eShelWgt(ttk.Frame) :

    def __init__(self,container) :
        super().__init__(container) 

        ttk.Label(self,text="ESHEL CALIBRATION",width=24).grid(column=1,row=1,sticky=(W))

        ttk.Label(self,text="MIRROR",width=8).grid(column=1,row=2,sticky=(W))
        self.mirror = StringVar()
        ttk.Label(self, textvariable=self.mirror).grid(column=2,row=2,sticky=(W,E),padx=10)

        ttk.Label(self,text="QUARTZ",width=8).grid(column=3,row=2,sticky=(W))
        self.quartz = StringVar()
        ttk.Label(self, textvariable=self.quartz).grid(column=4,row=2,sticky=(W,E),padx=10)

        ttk.Label(self,text="LED",width=8).grid(column=5,row=2,sticky=(W))
        self.led = StringVar()
        ttk.Label(self, textvariable=self.led).grid(column=6,row=2,sticky=(W,E),padx=10)

        ttk.Label(self,text="ThAr",width=8).grid(column=7,row=2,sticky=(W))
        self.thar = StringVar()
        ttk.Label(self, textvariable=self.thar).grid(column=8,row=2,sticky=(W,E),padx=10)

def status() :
    """ Start status window and updater
    """
    root = Tk()
    default_font=tkinter.font.nametofont('TkDefaultFont')
    default_font.configure(size=16,weight=tkinter.font.BOLD,family='Arial')

    root.rowconfigure(0, weight=1)
    root.columnconfigure(0, weight=1,minsize=200)
    root.minsize(500,20)
    root.geometry('+0+0')
    root.resizable(False,False)

    mainframe = ttk.Frame(root, padding="3 3 12 12",width=100,height=10)
    mainframe.grid(column=0, row=0, stick=(N,W,E,S))
    mainframe.columnconfigure(1, minsize=100, weight=1)

    telframe=TelescopeWgt(mainframe)
    telframe.grid(column=1,row=1,stick=(W))
    for i in range(1,7) :
        telframe.columnconfigure(i, minsize=50, weight=1)

    line=ttk.Separator(mainframe,orient='horizontal')
    line.grid(column=1,row=2,stick=(E,W))

    domeframe=DomeWgt(mainframe)
    domeframe.columnconfigure(1, minsize=200, weight=1)
    domeframe.grid(column=1,row=3,stick=(W))

    line=ttk.Separator(mainframe,orient='horizontal')
    line.grid(column=1,row=4,stick=(E,W))

    camframe=CameraWgt(mainframe)
    camframe.grid(column=1,row=5,stick=(W))

    line=ttk.Separator(mainframe,orient='horizontal')
    line.grid(column=1,row=6,stick=(E,W))

    iodineframe=IodineWgt(mainframe)
    iodineframe.grid(column=1,row=7,stick=(W))

    line=ttk.Separator(mainframe,orient='horizontal')
    line.grid(column=1,row=8,stick=(E,W))

    eShelframe=eShelWgt(mainframe)
    eShelframe.grid(column=1,row=9,stick=(W))

    for child in mainframe.winfo_children(): 
        child.grid_configure(padx=5, pady=5)
    mainframe.focus()

    shutter=['Open','Closed','Opening','Closing','Error']
    camerastate=['Idle','Waiting','Exposing','Reading','Download','Error']
    coverstate=['NotPresent','Closed','Moving','Open','Unknown','Error']

    apo=EarthLocation.of_site('APO')
    #aposite=site.Site('APO')

    def update() :
        try :
            # weather status to influxDB
            wdict=weather.getapo()
            weather.influx_write(wdict)
        except : pass

        try :
            ccd_dict={}
            for i in range(4) :
                try :
                    icam = aposong.getcam(i)
                    ccd_dict[f'camera_{i}_temp'] = aposong.C[icam].CCDTemperature
                    ccd_dict[f'camera_{i}_power'] = aposong.C[icam].CoolerPower
                    influx.write(ccd_dict,bucket='ccdtemp',measurement=f'ccd_{i}')
                except : continue
        except : pass

        try :
            dict={}
            dict['qhy_thermocouple_1'] = aposong.SW[2].GetSwitchValue(0)
            dict['qhy_thermocouple_2'] = aposong.SW[2].GetSwitchValue(1)
            influx.write(dict,bucket='ccdtemp',measurement='qhy_thermocouple')
        except  : pass

        try :
            # Get iodine cell related data
            pos = aposong.iodine_position()
            temp = aposong.iodine_tget()
            tset = aposong.iodine_tset()
            volt = aposong.iodine_get('voltage')
            curr = aposong.iodine_get('current')
            iodineframe.position.set(pos)
            iodineframe.temp.set(temp+' / '+tset)
            iodineframe.voltage.set(volt)
            iodineframe.current.set(curr)

            # parse temperature channels
            tset1,tset2=tset.split() 
            temp1,temp2=temp.split()
            volt1,volt2=volt.split()
            curr1,curr2=curr.split()
            if float(temp1)>float(tset1)+20 or float(temp2)>float(tset2)+20 :
                # if temp is more than 20 degrees above set temp, disable heaters!
                aposong.iodine_set('enable',0)

            # load into influx database
            iodine_dict={}
            for k,v in zip(['temp1','temp2','volt1','volt2','curr1','curr2'],
                           [temp1,temp2,volt1,volt2,curr1,curr2]) :
                iodine_dict[k] = float(v)
            influx.write(iodine_dict,bucket='iodinetemp',measurement='my_measurement')
        except : pass

        try :
            # get eShel calibration status
            state = ['Off','On']
            eShelframe.mirror.set(state[aposong.SW[1].GetSwitch(3)])
            eShelframe.quartz.set(state[aposong.SW[1].GetSwitch(0)])
            eShelframe.led.set(state[aposong.SW[1].GetSwitch(2)])
            eShelframe.thar.set(state[aposong.SW[1].GetSwitch(1)])
        except : pass

        try :
            t=Time.now()
            t.location=apo
            y,m,d,h,m,s=t.ymdhms
            telframe.ut.set('{:02d}:{:02d}:{:04.1f}'.format(h,m,s))
            h,m,s=t.sidereal_time('mean').hms
            telframe.lst.set('{:02d}:{:02d}:{:04.1f}'.format(int(h),int(m),s))
            telframe.mjd.set('{:.2f}'.format(t.mjd))

            stat = aposong.telescope_status()

            telframe.port.set('{:d}'.format(stat.m3.port))
            ra = stat.mount.ra_j2000_hours
            dec = stat.mount.dec_j2000_degs
            c = SkyCoord(ra=ra*u.h, dec=dec*u.degree)
            radec=c.to_string('hmsdms',sep=':',precision=1) 
            ras=radec.split()[0]
            decs=radec.split()[1]
            telframe.ra.set('{:s}'.format(ras))
            telframe.dec.set('{:s}'.format(decs))
            ha = (t.sidereal_time('mean')-ra*u.hourangle) 
            ha = (ha+12*u.hourangle)%(24*u.hourangle)-(12*u.hourangle)
            h,m,s=ha.hms
            telframe.ha.set('{:02d}:{:02d}:{:04.1f}'.format(int(h),int(m),s))

            telframe.az.set('{:.2f}'.format(stat.Azimuth))
            telframe.alt.set('{:.2f}'.format(stat.Altitude))

            telframe.rot.set('{:.1f}'.format(stat.rotator.mech_position_degs))
            telframe.pa.set('{:.1f}'.format(stat.rotator.field_angle_degs))

            telframe.focus.set('{:d}'.format(aposong.getfoc()))
        except : 
            telframe.ut.set('ERROR')

        try :
            domeframe.coverstate.set('{:s}'.format(coverstate[aposong.Covers.CoverState.value]))

            az, shutterstatus, slewing = aposong.domestatus()
            domeframe.az.set('{:.1f}'.format(az))
            domeframe.shutter.set('{:s}'.format(shutter[shutterstatus]))
            if slewing : domeframe.slewing.set('SLEWING')
            else : domeframe.slewing.set(' ')

            domeframe.stat35m.set(aposong.S.Action('stat35m')+'/'+aposong.S.Action('stat25m'))
            if domeframe.stat35m.get() == 'closed' : domeframe.statcolor.set('red')
            else : domeframe.statcolor.set('green')
        except : pass

        try :
            camframe.filter.set('{:s}'.format(aposong.filtname()))
            #camframe.binning.set('{:d}x{:d}'.format(C.BinX,C.BinY))
            #camframe.state.set('{:s}'.format(camerastate[C.CameraState]))
            camframe.temperature.set('{:.1f}'.format(ccd_dict['camera_0_temp']))
            camframe.cooler.set('{:.1f}'.format(C.CoolerPower))
        except : pass



        root.after(5000,update)

    import signal
    def handler(signum,frame):
        return
    signal.signal(signal.SIGINT,handler)
    update()
    root.mainloop()
