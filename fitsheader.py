from astropy.coordinates import Angle
import astropy.units as u

def camera(hdu,C,exptime,avg,light) :

    hdu.header['---C----'] = ('------CAMERA------','-------------------------------------')
    hdu.header['DETECTOR'] = (C.Name, 'Detector name')
    hdu.header['CCD_TEMP'] = (C.CCDTemperature, 'CCD Temperature (degrees C)')
    hdu.header['HOR_BIN'] = (C.BinX,'Horizontal binning')
    hdu.header['VER_BIN'] = (C.BinY,'Vertical binning')
    hdu.header['EXPTIME'] = (exptime, 'Exposure time in seconds')
    hdu.header['EXPAVG'] = avg
    if light: hdu.header['IMAGTYP'] = ('LIGHT','Image type')
    else: hdu.header['IMAGTYP'] = ('DARK','Image type')

def mixed(hdu) :
    hdu.header['---M----'] = ('-----MIXED--------','-------------------------------------')

def object (hdu,targ) :
    hdu.header['---O----'] = ('-----OBJECT-------','-------------------------------------')
    if targ is not None : 
        hdu.header['OBJECT'] = (targ, 'Name of object')

def spectrograph(hdu,foc) :
    hdu.header['---SP---'] = ('---SPECTROGRAPH---','-------------------------------------')
    if foc is not None :
        hdu.header['CAMFOCS'] = (foc, 'Spectrograph camera focus')
        hdu.header['SPECFOC'] = (foc, 'Spectrograph camera focus')

def telescope(hdu,stat,foc) :
    hdu.header['---TEL--'] = ('-----TELESCOPE----','-------------------------------------')
    hdu.header['RA'] = (Angle(stat.mount.ra_j2000_hours, unit=u.hour).to_string(sep=':',precision=2), 'Telescope RA (hh:mm:ss)')
    hdu.header['DEC'] = (Angle(stat.mount.dec_j2000_degs, unit=u.degree).to_string(sep=':',precision=1), 'Telescope DEC (dd:mm:ss)')
    hdu.header['AZ'] = (stat.mount.azimuth_degs, 'Telescope azimuth')
    hdu.header['ALT'] = (stat.mount.altitude_degs, 'Telescope altitude')
    hdu.header['ROT'] = (stat.rotator.mech_position_degs, 'Telescope port 1 rotator')
    if foc is not None :
        hdu.header['FOCUS'] = (foc, 'Nasymyth focus')
        hdu.header['TEL-FOC'] = (foc, 'Nasymyth focus')
    #hdu.header['TELESCOP'] = 'APO SONG 1m'

def weather(hdu) :
    hdu.header['---W----'] = ('-----WEATHER------','-------------------------------------')

def sunmoon(hdu) :
    hdu.header['---SM---'] = ('-----SUN/MOON-----','-------------------------------------')

def time(hdu,t) :
    hdu.header['---TI---'] = ('-----TIME---------','-------------------------------------')
    hdu.header['DATE-OBS'] = (t.fits, 'Start time of exposure in UTC')
    hdu.header['JD-DATE'] = (t.jd, 'Start time of exposure in Julian date')
    hdu.header['MJD-DATE'] =( t.mjd, 'Start time of exposure in modified Julian date')
