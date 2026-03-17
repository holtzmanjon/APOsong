from socket import *
import time
import pytz
from datetime import datetime, timezone
from astropy.time import Time
from astropy.table import Table
import influx
import database

def getapo() :
    """ Get APO weather status from network query

        Returns : dictionary with various quantities
    """
    sock = socket(AF_INET, SOCK_DGRAM)
    sock.sendto(b'all', ('10.75.0.152', 6251))
    messin, server = sock.recvfrom(1024)
    sock.close()
    messin=messin.decode()

    # strip off beginning.
    # replace = with =', and replace space with space'
    # for example dewPoint=14.0 becomes dewPoint='14.0'
    start = messin.find('timeStamp')
    stop = messin.find('end')
    stuff = messin[start:stop].replace ("=","='").replace (" ","'; ")
    ldict={}
    exec(stuff,globals(),ldict)
    return ldict

def manual_override() :
    d=database.DBSession(host='song1m.apo.nmsu.edu',database='apo',user='song')
    out=d.query(sql='SELECT pressed_at FROM public.button_presses ORDER BY id DESC LIMIT 1')
    override_time=datetime.fromisoformat(out[0][0])
    now=pytz.timezone('America/Denver').localize(datetime.now())
    if override_time> now :
        return int((override_time-now).total_seconds())
    else :
        return 0

def conditions_ok() :
    d=database.DBSession(host='song1m.apo.nmsu.edu',database='apo',user='song')
    out=d.query('public.button_presses')
    override_time=datetime.fromisoformat(out[0][1])
    now=pytz.timezone('America/Denver').localize(datetime.now())
    if now>override_time :
        print('override not enabled!')
        return False

    wdict=getapo()
    if (float(wdict['winds']) < 25 and float(wdict['winds2'])< 25 and float(wdict['gusts'])<30 and 
        float(wdict['humidity']) < 80 and (float(wdict['outside_airtemp_35m'])-float(wdict['outside_dewpoint_35m']))>3) :
        return True
    else :
        print('weather conditions not satisfied!')
        return False


def postgres_write(wdict) :
    """ Loads subset of weather quantities into table for postgres db loading
    """
    tab=Table()
    tab['ins_at'] = [Time.now().fits]
    tab['dateobs'] = [Time.now().fits]
    tab['temp'] = float(wdict['airtemp'])
    tab['hum'] = float(wdict['humidity'])
    tab['pressure'] = float(wdict['pressure'])
    tab['dewpoint'] = float(wdict['dewpoint'])
    tab['wind_speed'] = float(wdict['winds'])
    tab['wind_dir'] = float(wdict['windd'])
    if wdict['encl35m'] == 'open' : tab['encl35m'] = 1
    else : tab['encl35m'] = 0
    if wdict['encl25m'] == 'open' : tab['encl25m'] = 1
    else : tab['encl25m'] = 0

    d=database.DBSession(host='song1m_db.apo.nmsu.edu',database='db_apo',user='song')
    d.ingest('public.weather',tab,onconflict='update')
    d.close()
    return tab

def influx_write(wdict) :
    """ Loads subset of weather quantities into Influx database
    """
    influx_dict={}
    for key in ['airtemp','humidity','pressure','dewpoint','winds','windd'] :
        influx_dict[key] = float(wdict[key])
    influx.write(influx_dict,bucket='weather',measurement='my_measurement')

