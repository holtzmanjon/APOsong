import time
from alpaca.camera import *
import usbrelay
from yoctopuce 
from yoctopuce.yocto_api import *
from yoctopuce.yocto_temperature import *

errmsg=YRefParam()
YAPI.RegisterHub("usb",errmsg)

# open camera and relay
C=Camera('localhost:11111',0)
relay=usbrelay.USBRelay()
thermocouple=YTemperature.FindTemperature("THRMCPL1-286EEC.temperature1")
thermocouple.isOnline()

# check every minute
while True :
    # if CCDtemp < 30, reset watchdog with relay
    # if watchdog is not reset, set it to cut power in xxx minutes
    try:
        tccd = C.CCDTemperature
        thermocouple.get_currentValue()
        print('tccd: {:f}  thermocouple: {:f}', tccd,thermocouple)
        if temp < 30 :
            # reset watchdog
            print('resetting watchdog...')
            relay.on_relay(1)
            time.sleep(1)
            relay.on_relay(1)
        else :
            print('NOT resetting watchdog...')
    except:
        pass

    print('sleeping...')
    time.sleep(10)
