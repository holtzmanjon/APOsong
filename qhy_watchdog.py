import time
from alpaca.camera import *
from alpaca.switch import *
import influx
#import usbrelay

# open thermocouple
S=Switch('172.24.4.202:5555',0)
C=Camera('172.24.4.202:11111',0)
while True :
    try :
        dict={}
        dict['qhy_thermocouple_1'] = S.GetSwitchValue(0)
        dict['qhy_thermocouple_2'] = S.GetSwitchValue(1)
    except  :
        time.sleep(30)
        S=Switch('172.24.4.202:5555',0)
    try :
        dict['qhy_ccdtemp'] = C.CCDTemperature
        dict['qhy_setccdtemp'] = C.SetCCDTemperature
        dict['qhy_coolerpower'] = C.CoolerPower
    except  :
        time.sleep(30)
        C=Camera('172.24.4.202:11111',0)
    influx.write(dict,bucket='ccdtemp',measurement='qhy')
    print(dict)
    time.sleep(10)

# open camera and relay
#C=Camera('172.24.4.202:11111',0)
#C.CoolerOn = True
#C.SetCCDTemperature=-10

#relay=usbrelay.USBRelay()


# check every minute
#for i in range(10) :
#    # if CCDtemp < 30, reset watchdog with relay
#    # if watchdog is not reset, set it to cut power in xxx minutes
#    if True:
#        tccd = C.CCDTemperature
#        t1 =thermo1.get_currentValue()
#        t2 =thermo2.get_currentValue()
#        print('tccd: {:f} {:f} {:f} thermocouple: {:f} {:f}'.format(tccd,C.SetCCDTemperature,C.CoolerPower,t1,t2))
#        if t1 < 30 :
#            # reset watchdog
#            print('resetting watchdog...')
#            relay.on_relay(1)
#            time.sleep(1)
#            relay.off_relay(1)
#        else :
#            print('NOT resetting watchdog...')
#    #except:
#    #    pass
#
#    print('sleeping...')
#    time.sleep(10)
