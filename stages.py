import lts
import zaber
import socket
import serial

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
            input = data.decode()
            cmd = input.split()[0]
            try: val = input.split()[1]
            except : val = ''
            if cmd == 'iodine_pos' :
                if len(val) > 0 :
                    s1.move(int(val))
                return s1.get_position()
            elif cmd == 'focus' :
                if len(val) > 0 :
                    s2.move(int(val))
                return s2.get_position()
            elif cmd == 'iodine_temp' :
                if len(val) > 0 :
                    tc.write(b'TSET1=140\r')
                    tc.readline()
                    tc.write(b'TSET2=140\r')
                    tc.readline()
                tc.write(b'TACT1?\r')
                return tc.readline()

