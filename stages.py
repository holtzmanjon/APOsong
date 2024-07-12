import lts
import zaber
import socket
from serial import Serial

HOST = "10.75.0.202"  # Standard loopback interface address (localhost)
HOST = "172.24.4.202"  # Standard loopback interface address (localhost)
PORT = 65431  # Port to listen on (non-privileged ports are > 1023)

def remote() :
  """ Run simple remote socket server 
  """
  s1=lts.ThorlabsStage()
  s2=zaber.ZaberStage()
  tc=Serial('COM3',115200,timeout=1)

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
            input = data.decode()
            cmd = input.split()[0]
            try: val = input.split()[1]
            except : val = ''
            if cmd == 'iodine_pos' :
                if len(val) > 0 :
                    s1.move(int(val))
                conn.sendall(str(s1.get_position()).encode())
            elif cmd == 'focus' :
                if len(val) > 0 :
                    s2.move(int(val))
                conn.sendall(str(s2.get_position()).encode())
            elif cmd == 'iodine_tset' :
                if len(val) > 0 :
                    tc.write('TSET1={:s}\r'.format(val).encode())
                    tc.readline()
                    tc.write('TSET2={:s}\r'.format(val).encode())
                    tc.readline()
                tc.write(b'TSET1?\r')
                conn.sendall(tc.readline())
            elif cmd == 'iodine_tact' :
                tc.write(b'TACT1?\r')
                t1=tc.readline()
                tc.write(b'TACT2?\r')
                t2=tc.readline()
                conn.sendall(t1+b' '+t2)
    s1.close()
    s2.close()
    tc.close()

