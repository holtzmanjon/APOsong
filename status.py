#import tkthread
#tkthread.patch()
#tkthread.tkinstall(ensure_root=True)

from tkinter import *
from tkinter import ttk
import tkinter.font
import numpy as np
from numpy.random import randint
try :
    from alpaca.camera import *
    from alpaca.telescope import *
    from alpaca.covercalibrator import *
    from alpaca.dome import *
    from alpaca.focuser import *
    from alpaca.filterwheel import *
    from alpaca.safetymonitor import *
except :
    print('no alpaca in status')

from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
import astropy.units as u


#from coordio import sky, site, time, ICRS

class TelescopeWgt(ttk.Frame) :

    def __init__(self,container) :
        super().__init__(container)
        row=1
        ttk.Label(self,text="RA").grid(column=1,row=row,sticky=(W))
        self.ra = StringVar()
        ttk.Label(self, textvariable=self.ra).grid(column=2,row=row,sticky=(E),padx=10)

        ttk.Label(self,text="AZ").grid(column=3,row=row,sticky=(W))
        self.az = StringVar()
        ttk.Label(self, textvariable=self.az).grid(column=4,row=row,sticky=(E),padx=10) 

        ttk.Label(self,text="UT").grid(column=5,row=row,sticky=(W))
        self.ut = StringVar()
        ttk.Label(self, textvariable=self.ut).grid(column=6,row=row,sticky=(E),padx=10)


        row+=1
        ttk.Label(self,text="DEC").grid(column=1,row=row,sticky=(W))
        self.dec = StringVar()
        ttk.Label(self, textvariable=self.dec).grid(column=2,row=row,sticky=(E),padx=10)

        ttk.Label(self,text="ALT").grid(column=3,row=row,sticky=(W))
        self.alt = StringVar()
        ttk.Label(self, textvariable=self.alt).grid(column=4,row=row,sticky=(E),padx=10)

        ttk.Label(self,text="LST").grid(column=5,row=row,sticky=(W))
        self.lst = StringVar()
        ttk.Label(self, textvariable=self.lst).grid(column=6,row=row,sticky=(E),padx=10)

        row+=1
        ttk.Label(self,text="PA").grid(column=1,row=row,sticky=(W))
        self.pa = StringVar()
        ttk.Label(self, textvariable=self.pa).grid(column=2,row=row,sticky=(E),padx=10)

        ttk.Label(self,text="ROT").grid(column=3,row=row,sticky=(W))
        self.rot = StringVar()
        ttk.Label(self, textvariable=self.rot).grid(column=4,row=row,sticky=(E),padx=10)

        ttk.Label(self,text="MJD").grid(column=5,row=row,sticky=(W))
        self.mjd = StringVar()
        ttk.Label(self, textvariable=self.mjd).grid(column=6,row=row,sticky=(E),padx=10)

        row+=1
        ttk.Label(self,text="FOCUS").grid(column=1,row=row,sticky=(W))
        self.focus = StringVar()
        ttk.Label(self, textvariable=self.focus).grid(column=2,row=row,sticky=(E),padx=10) 

class DomeWgt(ttk.Frame) :

    def __init__(self,container) :
        super().__init__(container) 

        ttk.Label(self,text="SHUTTER ",width=10).grid(column=1,row=1,sticky=(W))
        self.shutter = StringVar()
        ttk.Label(self, textvariable=self.shutter).grid(column=2,row=1,sticky=(W),padx=10)

        ttk.Label(self,text="DOME AZ",width=10).grid(column=3,row=1,sticky=(W))
        self.az = StringVar()
        ttk.Label(self, textvariable=self.az).grid(column=4,row=1,sticky=(W),padx=10)

        self.slewing = StringVar()
        ttk.Label(self, textvariable=self.slewing).grid(column=5,row=1,sticky=(W),padx=10)
     
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


def status(pwi=None, T=None, D=None, Filt=None, F=None, C=None) :
    """ Start status window and updater
    """
    root = Tk()
    default_font=tkinter.font.nametofont('TkDefaultFont')
    default_font.configure(size=12,weight=tkinter.font.BOLD)

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

    for child in mainframe.winfo_children(): 
        child.grid_configure(padx=5, pady=5)
    mainframe.focus()

    shutter=['Open','Closed','Opening','Closing','Error']
    state=['Idle','Waiting','Exposing','Reading','Download','Error']

    apo=EarthLocation.of_site('APO')
    #aposite=site.Site('APO')

    def update() :
        if True :
            t=Time.now()
            t.location=apo
            y,m,d,h,m,s=t.ymdhms
            telframe.ut.set('{:02d}:{:02d}:{:04.1f}'.format(h,m,s))
            h,m,s=t.sidereal_time('mean').hms
            telframe.lst.set('{:02d}:{:02d}:{:04.1f}'.format(int(h),int(m),s))
            telframe.mjd.set('{:.2f}'.format(t.mjd))

            if T is None and pwi is None :
                telframe.ra.set('N/A')
                telframe.dec.set('N/A')
            else :
                if pwi is None :
                    ra = T.RightAscension
                    dec = T.Declination
                # convert from topocentric to ICRS
                #aposite.set_time(time.Time())
                ##print('telescope: ', T.RightAscension, T.Declination)
                ##print('time: ', time.Time())
                ##obs = sky.Observed([[15*T.RightAscension, T.Declination]], site=aposite)
                #obs = sky.Observed([[T.Altitude, T.Azimuth]], site=aposite)
                #icrs = sky.ICRS(obs)
                ##print('ICRS:',icrs)
                #telframe.ra.set('{:f}'.format(icrs[0][0]/15.))
                #telframe.dec.set('{:f}'.format(icrs[0][1]))
                else :
                    stat = pwi.status()
                    ra = stat.mount.ra_j2000_hours
                    dec = stat.mount.dec_j2000_degs
                c = SkyCoord(ra=ra*u.h, dec=dec*u.degree)
                radec=c.to_string('hmsdms',sep=':',precision=1) 
                ras=radec.split()[0]
                decs=radec.split()[1]
                telframe.ra.set('{:s}'.format(ras))
                telframe.dec.set('{:s}'.format(decs))

            if T is not None :
                telframe.az.set('{:.2f}'.format(T.Azimuth))
                telframe.alt.set('{:.2f}'.format(T.Altitude))
            else :
                telframe.az.set('N/A')
                telframe.alt.set('N/A')

            if pwi is not None :
                telframe.rot.set('{:.1f}'.format(stat.rotator.mech_position_degs))
                telframe.pa.set('{:.1f}'.format(stat.rotator.field_angle_degs))
            else :
                telframe.rot.set('missing')
                telframe.pa.set('missing')

            if F is not None :
                telframe.focus.set('{:d}'.format(F.Position))
            else :
                telframe.focus.set('N/A')

        if D is not None :
            domeframe.az.set('{:.1f}'.format(D.Azimuth))
            domeframe.shutter.set('{:s}'.format(shutter[D.ShutterStatus]))
            if D.Slewing : domeframe.slewing.set('SLEWING')
            else : domeframe.slewing.set(' ')
        else :
            domeframe.az.set('N/A')
            domeframe.shutter.set('N/A')
            domeframe.slewing.set('N/A')

        if Filt is not None :
            camframe.filter.set('{:s}'.format(Filt.Names[Filt.Position]))
        else :
            camframe.filter.set('None')

        if C is not None :
            camframe.binning.set('{:d}x{:d}'.format(C.BinX,C.BinY))
            camframe.state.set('{:s}'.format(state[C.CameraState]))
            camframe.temperature.set('{:.1f}/{:.1f}'.format(
                     C.CCDTemperature,C.SetCCDTemperature))
            camframe.cooler.set('{:.1f}'.format( C.CoolerPower))
        else :
            camframe.filter.set('N/A')
            camframe.binning.set('N/A')
            camframe.state.set('N/A')
            camframe.temperature.set('N/A')
            camframe.cooler.set('N/A')
        root.after(1000,update)

    update()
    root.mainloop()
