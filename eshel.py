# routines for use with Shelyak eShel spectograph and calibration channel

import socket
import sys, getopt
try :
    from serial import Serial
except :
    print('failed to import Serial, OK if client')
    pass

from time import sleep

try: import aposong
except :
    print("failed to import aposong, OK if server")
    pass

HOST = "10.75.0.22"  # Standard loopback interface address (localhost)
PORT = 65432  # Port to listen on (non-privileged ports are > 1023)

esock = None

def lamps(close=False,mirror=False,thar=False,quartz=False,led=False) :
    """ Send commands to remote socket for eShel calibration control
    """
    global esock
    if esock is None :
        print('if this hangs, make sure remote program is running')
        esock = socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
        esock.connect(('10.75.0.220', 65432))

    val =  mirror<<7 | led<<6 | thar<<5 | quartz<<4 
  
    try : esock.sendall(str(val).encode())
    except :
        print("communication error: is remote program running?")
    if close :
        esock.shutdown(socket.SHUT_RDWR)
        esock.close()

def cals(display=None,flats=15,thar=15,flat_exptime=8,thar_exptime=10) :
    """ Take series of eShel cals
    """
    lamps(mirror=True,quartz=True,led=True)
    for i in range(flats) :
        aposong.expose(flat_exptime,filt=None,bin=1,display=display,cam=1,name='eShel_flat')
    lamps(mirror=True,thar=True)
    for i in range(thar) :
        aposong.expose(thar_exptime,filt=None,bin=1,display=display,cam=1,name='eShel_ThAr')
    lamps()

def remote() :
  """ Run simple remote socket server to send commands to cal lamp controller
  """
  relaycard=K8056('COM4')
  relaycard.clear(9)
  with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
    s.bind((HOST, PORT))
    while True :
      s.listen()
      conn, addr = s.accept()
      with conn:
        print(f"Connected by {addr}")
        while True:
            data = conn.recv(1024)
            if not data:
                print('done')
                break
            print('received: ', data)
            #conn.sendall(data)
            relaycard.send_byte(int(data))


#    Copyright (c) 2007, 2015 by jerch (brejoe@web.de)
#
#    This program is free software; you can redistribute it and#or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
# -*- coding: utf-8 -*-


class K8056(object):
    """
    K8056 - Class for controlling the velleman K8056 8 channel relay card
    K8056(device, repeat=0, wait=0)
    Give serial port as number or device file.
    For better reliability repeat instructions `repeat` times
    and `wait` seconds between execution.
    """
    def __init__(self, device, repeat=0, wait=0):
        self._serial = Serial(device, 2400)
        self.repeat = repeat
        self.wait = wait
        sleep(0.1)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        '''Close underlying serial device connection.'''
        self._serial.close()
        
    def _process(self, instruction, byte, address):
        cksum = (243 - instruction - byte - address) & 255
        for i in range(self.repeat + 1):
            self._serial.write(bytearray([13, address, instruction, byte, cksum]))
            sleep(self.wait)
        
    def set(self, relay, address=1):
        '''Set `relay` (9 for all) of card at `address` (default 1).'''
        if not 0 < relay < 10:
            raise Exception('invalid relay number')
        self._process(83, relay+48, address&255)
        
    def clear(self, relay, address=1):
        '''Clear `relay` (9 for all) of card at `address` (default 1).'''
        if not 0 < relay < 10:
            raise Exception('invalid relay number')
        self._process(67, relay+48, address&255)
        
    def toggle(self, relay, address=1):
        '''Toggle `relay` (9 for all) of card at `address` (default 1).'''
        if not 0 < relay < 10:
            raise Exception('invalid relay number')
        self._process(84, relay+48, address&255)
        
    def set_address(self, new=1, address=1):
        '''Set address of card at `address` to `new`.'''
        self._process(65, new&255, address&255)
        
    def send_byte(self, num, address=1):
        '''Set relays to `num` (in binary mode).'''
        self._process(66, num&255, address&255)

    def emergency_stop(self):
        '''Clear all relays on all cards. emergency purpose.'''
        self._process(69, 1, 1)
        
    def force_address(self):
        '''Reset all cards to address 1.'''
        self._process(70, 1, 1)
        
    def get_address(self):
        '''Display card address on LEDs.'''
        self._process(68, 1, 1)

def main(argv):
    m = 256
    l = 256
    t = 256
    n = 256
    byte_val = m | l | t | n        # OR the 4 values to one value

    if(len(argv) == 0):
        print("Usage: eShelCalLamps.py --mirror=<0|1> --led=<0|1> --thar=<0|1> --neon=<0|1>")
        print("  or : eShelCalLamps.py -m <0|1> -l <0|1> -t <0|1> -n <0|1>")
        sys.exit(2)
    # print('Commandline Arguments: ', sys.argv[1:])
    try:
        opts, args = getopt.getopt(argv, "hm:l:t:n:", ["help","mirror=","led=","thar=","neon="])
    except getopt.GetoptError:
        print("Usage: eShelCalLamps.py --mirror=<0|1> --led=<0|1> --thar=<0|1> --neon=<0|1>")
        print("  or : eShelCalLamps.py -m <0|1> -l <0|1> -t <0|1> -n <0|1>")
        sys.exit(2)
    # print('opts: ', opts, 'args: ', args)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("Usage: eShelCalLamps.py --mirror=<0|1> --led=<0|1> --thar=<0|1> --neon=<0|1>")
            print("  or : eShelCalLamps.py -m <0|1> -l <0|1> -t <0|1> -n <0|1>")
            sys.exit()
        if opt in ("-m", "--mirror"):
            m = int(arg) << 7         # bitshift to _000 0000
        if opt in ("-l", "--led"):
            l = int(arg) << 6         # bitshift to 0_00 0000
        if opt in ("-t", "--thar"):
            t = int(arg) << 5         # bitshift to 00_0 0000
        if opt in ("-n", "--neon"):
            n = int(arg) << 4         # bitshift to 000_ 0000
        byte_val = m | l | t | n      # OR the 4 values to one value
 
    # print('arg: ', m, l, t, n, 'Byte Value: ', byte_val)

    # Using the Keyspan USB to serial converter
    # device = '/dev/tty.KeySerial1'

    # Using COMM / serial port 2
    device = 'COM4'

    with K8056(device) as relaycard:
        # clear all relays
        relaycard.clear(9)
        sleep(1)
        relaycard.send_byte(int(byte_val))
        sleep(1)

if __name__ == '__main__':
    main(sys.argv[1:])


