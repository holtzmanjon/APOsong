import time
from alpaca.camera import *
import usbrelay

# open camera and relay
C=Camera('localhost:11111',0)
relay=usbrelay.USBRelay()

# check every minute
while True :
    # if CCDtemp < 30, reset watchdog with relay
    # if watchdog is not reset, set it to cut power in xxx minutes
    try:
        temp = C.CCDTemperature
        if temp < 30 :
            # reset watchdog
            print('resetting watchdog...')
            relay.on_relay(1)
            time.sleep(1)
            relay.on_relay(1)
    except:
        pass

    time.sleep(10)
