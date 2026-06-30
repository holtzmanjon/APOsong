from astropy.coordinates import Angle
import astropy.units as u
import spectemps
import weather as weather_apo

def camera(hdu,C,exptime,avg,light) :

    hdu.header['---C----'] = ('------CAMERA------','-------------------------------------')
    hdu.header['DETECTOR'] = (C.Name, 'Detector name')
    hdu.header['CCD_TEMP'] = (C.CCDTemperature, 'CCD Temperature (degrees C)')
    hdu.header['HOR_BIN'] = (C.BinX,'Horizontal binning')
    hdu.header['VER_BIN'] = (C.BinY,'Vertical binning')
    hdu.header['EXPTIME'] = (exptime, 'Exposure time in seconds')
    hdu.header['EXPAVG'] = avg

def mixed(hdu,header) :
    hdu.header['---M----'] = ('-----MIXED--------','-------------------------------------')
    for key,val in zip(header.keys(),header.values()) :
        hdu.header[key] = val

def object (hdu,targ) :
    hdu.header['---O----'] = ('-----OBJECT-------','-------------------------------------')
    if targ is not None : 
        hdu.header['OBJECT'] = (targ, 'Name of object')

def spectrograph(hdu,foc) :
    hdu.header['---SP---'] = ('---SPECTROGRAPH---','-------------------------------------')
    if foc is not None :
        hdu.header['CAMFOCS'] = (foc, 'Spectrograph camera focus')
        hdu.header['SPECFOC'] = (foc, 'Spectrograph camera focus')
    sdict = spectemps.get()
    try :
        hdu.header['TEMPA' ] = (float('{:.2f}'.format(sdict['tempagera'])), 'Temperature A outside of box')
        hdu.header['UPPER' ] = (float('{:.2f}'.format(sdict['uppertable'])), 'Temperature at upper table')
        hdu.header['OAP1' ] = (float('{:.2f}'.format(sdict['OAP1'])), 'Temperature at OAP1')
    except : pass
    try :
        hdu.header['TEMPB' ] = (float('{:.2f}'.format(sdict['tempagerb'])), 'Temperature B outside of box')
        hdu.header['GRATING' ] = (float('{:.2f}'.format(sdict['grating'])), 'Temperature at grating')
        hdu.header['CAMERA' ] = (float('{:.2f}'.format(sdict['camera'])), 'Temperature at camera')
    except : pass

def telescope(hdu,stat,foc,domeaz) :
    hdu.header['---TEL--'] = ('-----TELESCOPE----','-------------------------------------')
    hdu.header['RA'] = (Angle(stat.mount.ra_j2000_hours, unit=u.hour).to_string(sep=':',precision=2), 'Telescope RA (hh:mm:ss)')
    hdu.header['DEC'] = (Angle(stat.mount.dec_j2000_degs, unit=u.degree).to_string(sep=':',precision=1), 'Telescope DEC (dd:mm:ss)')
    hdu.header['AZ'] = (stat.mount.azimuth_degs, 'Telescope azimuth')
    hdu.header['ALT'] = (stat.mount.altitude_degs, 'Telescope altitude')
    hdu.header['ROT'] = (stat.rotator.mech_position_degs, 'Telescope port 1 rotator')
    if foc is not None :
        hdu.header['FOCUS'] = (foc, 'Nasymyth focus')
    hdu.header['DOMEAZ'] = (domeaz, 'Azimuth of dome')
    #hdu.header['TELESCOP'] = 'APO SONG 1m'

def fpu(hdu,pos,i2pos,temp1,temp2,calpos,SW,filt) :
    hdu.header['---FPU--'] = ('---FOCALPLANE-----','-------------------------------------')
    hdu.header['I_POS'] = (float(pos), 'Iodine stage position')
    hdu.header['I2POS'] = (i2pos, 'Code for APO iodine cell')
    hdu.header['I_TEMP1'] = (float(temp1), 'First iodine temperature')
    hdu.header['I_TEMP2'] = (float(temp2), 'Second iodine temperature')
    hdu.header['CAL_POS'] = (calpos, 'Calibration stage position')
    hdu.header['TUNGSTEN'] = (int(SW.GetSwitch(0)),'Calibration source tungstenThAr on/off')
    hdu.header['LED'] = (int(SW.GetSwitch(2)),'Calibration stage LEDThAr on/off')
    hdu.header['THAR'] = (int(SW.GetSwitch(1)),'Calibration stage ThAr on/off')
    hdu.header['FILTER'] = (filt,'Filter')

def weather(hdu) :
    hdu.header['---W----'] = ('-----WEATHER------','-------------------------------------')
    wdict=weather_apo.getapo()
    hdu.header['AIRTEMP'] = (float(wdict['airtemp']),'Air temperature from weather station')
    hdu.header['HUMIDITY'] = (float(wdict['humidity']),'Humidity from weather station')
    hdu.header['WIND'] = (float(wdict['winds']),'Wind speed from weather station')
    hdu.header['WIND_DIR'] = (float(wdict['windd']),'Wind direction from weather station')
    hdu.header['DEWPOINT'] = (float(wdict['dewpoint']),'Dewpoint from weather station')
    hdu.header['PRESSURE'] = (float(wdict['pressure']),'Pressure from weather station')

def sunmoon(hdu) :
    hdu.header['---SM---'] = ('-----SUN/MOON-----','-------------------------------------')

def time(hdu,t) :
    hdu.header['---TI---'] = ('-----TIME---------','-------------------------------------')
    hdu.header['DATE-OBS'] = (t.fits, 'Start time of exposure in UTC')
    hdu.header['JD-DATE'] = (t.jd, 'Start time of exposure in Julian date')
    hdu.header['MJD-DATE'] =( t.mjd, 'Start time of exposure in modified Julian date')
