from socket import *
import time
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

