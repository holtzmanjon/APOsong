import time
from alpaca.camera import *
import usbrelay

from yoctopuce.yocto_api import *
from yoctopuce.yocto_temperature import *

errmsg=YRefParam()
YAPI.RegisterHub("usb",errmsg)

# open camera and relay
C=Camera('172.24.4.202:11111',0)
C.CoolerOn = True
C.SetCCDTemperature=-10
relay=usbrelay.USBRelay()
thermo1=YTemperature.FindTemperature("THRMCPL1-286EEC.temperature1")
thermo1.isOnline()
thermo2=YTemperature.FindTemperature("THRMCPL1-286EEC.temperature2")


# check every minute
for i in range(10) :
    # if CCDtemp < 30, reset watchdog with relay
    # if watchdog is not reset, set it to cut power in xxx minutes
    if True:
        tccd = C.CCDTemperature
        t1 =thermo1.get_currentValue()
        t2 =thermo2.get_currentValue()
        print('tccd: {:f} {:f} {:f} thermocouple: {:f} {:f}'.format(tccd,C.SetCCDTemperature,C.CoolerPower,t1,t2))
        if t1 < 30 :
            # reset watchdog
            print('resetting watchdog...')
            relay.on_relay(1)
            time.sleep(1)
            relay.off_relay(1)
        else :
            print('NOT resetting watchdog...')
    #except:
    #    pass

    print('sleeping...')
    time.sleep(10)
