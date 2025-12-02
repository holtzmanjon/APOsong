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
from astropy.table import Table
import astropy.units as u
import pwi4_client

import weather
import spectemps
import influx
import aposong
import database

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
        self.stat35m_label = ttk.Label(self, textvariable=self.stat35m)
        self.stat35m_label.grid(column=4,row=row,sticky=(W,E),padx=10)

     
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
        ttk.Label(self,text="COOLER POWER",width=15).grid(column=3,row=2,sticky=(W))
        self.cooler = StringVar()
        ttk.Label(self, textvariable=self.cooler).grid(column=4,row=2,sticky=(W,E),padx=10)

        ttk.Label(self,text="TEMP",width=8).grid(column=5,row=2,sticky=(W))
        self.spec_temp = StringVar()
        ttk.Label(self, textvariable=self.spec_temp).grid(column=6,row=2,sticky=(W,E),padx=10)
        ttk.Label(self,text="COOLER POWER",width=15).grid(column=7,row=2,sticky=(W))
        self.spec_cooler = StringVar()
        ttk.Label(self, textvariable=self.spec_cooler).grid(column=8,row=2,sticky=(W,E),padx=10)

class IodineWgt(ttk.Frame) :

    def __init__(self,container) :
        super().__init__(container) 

        ttk.Label(self,text="IODINE STAGE",width=16).grid(column=1,row=1,sticky=(W))
        self.iodinestage = StringVar()
        self.iodinestage_label = ttk.Label(self, textvariable=self.iodinestage)
        self.iodinestage_label.grid(column=2,row=1,sticky=(W,E),padx=10)

        ttk.Label(self,text="IODINE TEMP",width=16).grid(column=3,row=1,sticky=(W))
        self.temp = StringVar()
        ttk.Label(self, textvariable=self.temp).grid(column=4,row=1,sticky=(W,E),padx=10)

        ttk.Label(self,text="IODINE VOLTAGE",width=16).grid(column=1,row=2,sticky=(W))
        self.voltage = StringVar()
        ttk.Label(self, textvariable=self.voltage).grid(column=2,row=2,sticky=(W,E),padx=10)

        ttk.Label(self,text="IODINE CURRENT",width=16).grid(column=3,row=2,sticky=(W))
        self.current = StringVar()
        ttk.Label(self, textvariable=self.current).grid(column=4,row=2,sticky=(W,E),padx=10)

class calWgt(ttk.Frame) :

    def __init__(self,container) :
        super().__init__(container) 

        ttk.Label(self,text="CALSTAGE",width=8).grid(column=1,row=1,sticky=(W))
        self.calstage = StringVar()
        self.calstage_label = ttk.Label(self, textvariable=self.calstage)
        self.calstage_label.grid(column=2,row=1,sticky=(W,E),padx=10)

        ttk.Label(self,text="MIRROR",width=8).grid(column=1,row=2,sticky=(W))
        self.mirror = StringVar()
        ttk.Label(self, textvariable=self.mirror).grid(column=2,row=2,sticky=(W,E),padx=10)

        ttk.Label(self,text="QUARTZ",width=8).grid(column=3,row=2,sticky=(W))
        self.quartz = StringVar()
        self.quartz_label = ttk.Label(self, textvariable=self.quartz)
        self.quartz_label.grid(column=4,row=2,sticky=(W,E),padx=10)

        ttk.Label(self,text="LED",width=8).grid(column=5,row=2,sticky=(W))
        self.led = StringVar()
        self.led_label = ttk.Label(self, textvariable=self.led)
        self.led_label.grid(column=6,row=2,sticky=(W,E),padx=10)

        ttk.Label(self,text="ThAr",width=8).grid(column=7,row=2,sticky=(W))
        self.thar = StringVar()
        self.thar_label = ttk.Label(self, textvariable=self.thar)
        self.thar_label.grid(column=8,row=2,sticky=(W,E),padx=10)

def postgres_bool(val) :
    return '1' if val else '0'

def postgres_write(telstatus,domestatus) :

    tab=Table()
    tab['tel_con_state'] = [postgres_bool(telstatus.mount.is_connected)]
    tab['tel_tracking'] = [postgres_bool(telstatus.mount.is_tracking)]
    tab['tel_ra_j2000'] = [telstatus.mount.ra_j2000_hours]
    tab['tel_dec_j2000'] = [telstatus.mount.dec_j2000_degs]
    tab['tel_ra'] = [telstatus.mount.ra_apparent_hours]
    tab['tel_dec'] = [telstatus.mount.dec_apparent_degs]
    tab['tel_alt'] = [telstatus.mount.altitude_degs]
    tab['tel_azm'] = [telstatus.mount.azimuth_degs]
    tab['tel_alt_rms_error'] = [telstatus.mount.axis1.rms_error_arcsec]
    tab['tel_azm_rms_error'] = [telstatus.mount.axis0.rms_error_arcsec]
    tab['m3_pos'] = [telstatus.m3.port]
    tab['dome_shutterstate'] = [domestatus.shutterstate]
    tab['dome_az'] = [domestatus.az]
    tab['dome_slewing'] = [postgres_bool(domestatus.slewing)]
    tab['dome_light_state'] = [postgres_bool(domestatus.lights)]
    tab['temp_m1'] = [None]
    tab['temp_m2'] = [None]
    tab['temp_m3'] = [None]
    tab['temp_back'] = [None]
    tab['temp_amb'] = [None]
    tab['focuser_1_pos'] = [aposong.foc()]
    tab['focuser_1_moving'] = [postgres_bool(aposong.F[0].IsMoving)]
    tab['focuser_2_pos'] = [aposong.specfoc()]
    tab['focuser_2_moving'] = [postgres_bool(aposong.F[4].IsMoving)]
    tab['tel_lst'] = [telstatus.site.lmst_hours]
    d=database.DBSession(host='localhost',database='db_apo')
    d.ingest('public.tel_dome',tab,onconflict='update')
    d.close()
    return tab


niter=0

if __name__ == '__main__' :
    """ Start status window and updater
    """
    aposong.init()

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

    calframe=calWgt(mainframe)
    calframe.grid(column=1,row=9,stick=(W))

    line=ttk.Separator(mainframe,orient='horizontal')
    line.grid(column=1,row=10,stick=(E,W))

    for child in mainframe.winfo_children(): 
        child.grid_configure(padx=5, pady=5)
    mainframe.focus()

    shutter=['Open','Closed','Opening','Closing','Error']
    camerastate=['Idle','Waiting','Exposing','Reading','Download','Error']
    coverstate=['NotPresent','Closed','Moving','Open','Unknown','Error']

    apo=EarthLocation.of_site('APO')

    def update() :
        global niter
        niter=(niter+1)%60
        try :
            # weather status to influxDB
            if niter%10 == 1 :
                wdict=weather.getapo()
                weather.influx_write(wdict)
                weather.postgres_write(wdict)
        except : print('error with weather')

        try :
            # spectrograph temperatures to influxDB
            if niter%10 == 1 :
                sdict=spectemps.get()
        except : print('error with spectemps')

        try :
            if niter%60 == 1 :
                ccd_dict={}
                for i in [0,3] :
                    try :
                        icam = aposong.getcam(i)
                        ccd_dict[f'camera_{i}_temp'] = aposong.C[icam].CCDTemperature
                        ccd_dict[f'camera_{i}_power'] = aposong.C[icam].CoolerPower
                        influx.write(ccd_dict,bucket='ccdtemp',measurement=f'ccd_{i}')
                    except : 
                        print('error with camera: ',i)
                        continue
                #camframe.filter.set('{:s}'.format(aposong.filtname()))
                #camframe.binning.set('{:d}x{:d}'.format(C.BinX,C.BinY))
                #camframe.state.set('{:s}'.format(camerastate[C.CameraState]))
                camframe.temperature.set('{:.1f}'.format(ccd_dict['camera_0_temp']))
                camframe.cooler.set('{:.1f}'.format(ccd_dict['camera_0_power']))
                camframe.spec_temp.set('{:.1f}'.format(ccd_dict['camera_3_temp']))
                camframe.spec_cooler.set('{:.1f}'.format(ccd_dict['camera_3_power']))
        except : print('Error with camera')

        try :
            dict={}
            dict['qhy_thermocouple_1'] = aposong.SW[2].GetSwitchValue(0)
            dict['qhy_thermocouple_2'] = aposong.SW[2].GetSwitchValue(1)
            influx.write(dict,bucket='ccdtemp',measurement='qhy_thermocouple')
        except  : pass

        try :
            # Get iodine cell related data
            pos = aposong.iodine_position()
            temp = aposong.iodine_tget().replace('actual temperature: ','')
            tset = aposong.iodine_tset().replace('set temperature: ','')
            volt = aposong.iodine_get('voltage')
            curr = aposong.iodine_get('current')
            iodineframe.iodinestage.set(pos)
            if abs(pos-aposong.config['iodinestage_in_pos']) < 0.2 : iodineframe.iodinestage_label.config(foreground='green3')
            elif abs(pos-aposong.config['iodinestage_out_pos']) < 0.2 : iodineframe.iodinestage_label.config(foreground='blue')
            else : iodineframe.iodinestage_label.config(foreground='yellow')
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
        except : print('error with iodine')

        try :
            pos =aposong.calstage_position()
            calframe.calstage.set(pos)
            if abs(pos-aposong.config['calstage_in_pos']) < 0.2 : calframe.calstage_label.config(foreground='green3')
            elif abs(pos-aposong.config['calstage_out_pos']) < 0.2 : calframe.calstage_label.config(foreground='blue')
            else : calframe.calstage_label.config(foreground='yellow')
            # get eShel calibration status
            state = ['Off','On']
            calframe.mirror.set(state[aposong.SW[1].GetSwitch(3)])
            calframe.quartz.set(state[aposong.SW[1].GetSwitch(0)])
            if state[aposong.SW[1].GetSwitch(0)] == 'On' : calframe.quartz_label.config(foreground='red')
            else :calframe.quartz_label.config(foreground='black')
            calframe.led.set(state[aposong.SW[1].GetSwitch(2)])
            if state[aposong.SW[1].GetSwitch(2)] == 'On' : calframe.led_label.config(foreground='red')
            else :calframe.led_label.config(foreground='black')
            calframe.thar.set(state[aposong.SW[1].GetSwitch(1)])
            if state[aposong.SW[1].GetSwitch(1)] == 'On' : calframe.thar_label.config(foreground='red')
            else :calframe.thar_label.config(foreground='black')
        except : print('error with cal')

        try :
            t=Time(Time.now(),location=apo)
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
            
            telframe.focus.set('{:d}'.format(aposong.foc()))

        except : 
            telframe.ra.set('ERROR')
            telframe.dec.set('ERROR')
            telframe.pa.set('ERROR')
            telframe.ha.set('ERROR')
            telframe.az.set('ERROR')
            telframe.alt.set('ERROR')
            telframe.rot.set('ERROR')
            telframe.port.set('ERROR')
            telframe.focus.set('{:d}'.format(aposong.foc(port=2)))

        try :
            domestat = aposong.domestatus()
            domeframe.az.set('{:.1f}'.format(domestat.az))
            domeframe.shutter.set('{:s}'.format(shutter[domestat.shutterstate]))
            if domestat.slewing : domeframe.slewing.set('SLEWING')
            else : domeframe.slewing.set(' ')

            stat35m = aposong.S.Action('stat35m')
            stat25m = aposong.S.Action('stat25m')
            domeframe.stat35m.set(stat35m+'/'+stat25m)
            if stat35m == 'open' or stat25m == 'open' : 
                domeframe.stat35m_label.config(foreground='green3')
            else :
                domeframe.stat35m_label.config(foreground='red')

            domeframe.coverstate.set('{:s}'.format(coverstate[aposong.Covers.CoverState.value]))
        except : print('error with dome')

        try : postgres_write(stat,domestat)
        except : pass
        root.after(5000,update)

    import signal
    def handler(signum,frame):
        return
    signal.signal(signal.SIGINT,handler)
    root.after(1000,update)
    #update()
    root.mainloop()

