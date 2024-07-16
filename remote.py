try : 
    import lts
except : 
    print('no lts')
    lts = None
try :
    import zaber
except : 
    print('no zaber')
    zaber = None
try :
    from serial import Serial
except : 
    print('no serial')
    serial = None
import socket

HOST = "10.75.0.202"  # Standard loopback interface address (localhost)
HOST = "172.24.4.202"  # Standard loopback interface address (localhost)
PORT = 65431  # Port to listen on (non-privileged ports are > 1023)

def server() :
  """ Run simple remote socket server 
  """
  if lts is not None :
      lts_stage=lts.ThorlabsStage()
  else :
      lts_stage = None
  if zaber is not None :
      zaber_stage=zaber.ZaberStage()
  else :
      zaber_stage = None
  if serial is not None :
      tc300=Serial('COM3',115200,timeout=1)
  else :
      tc300 = None

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
            device = input.split()[0]
            cmd = input.split()[1]
            try: val = input.split()[2]
            except : val = ''
            if device == 'lts' :
                if lts is None : 
                    conn.sendall(b'Error')
                elif cmd == 'position' :
                    if len(val) > 0 :
                        lts_stage.move(int(val))
                    conn.sendall(str(lts_stage.get_position()).encode())
            elif device == 'zaber' :
                if zaber is None : 
                    conn.sendall(b'Error')
                if cmd == 'position' :
                    if len(val) > 0 :
                        zaber_stage.move(int(val))
                conn.sendall(str(zaber_stage.get_position()).encode())
            elif device == 'tc300' :
                if cmd == 'tset' :
                    if len(val) > 0 :
                        tc300.write('TSET1={:s}\r'.format(val).encode())
                        tc300.readline()
                        tc300.write('TSET2={:s}\r'.format(val).encode())
                        tc300.readline()
                    tc300.write(b'TSET1?\r')
                    conn.sendall(tc300.readline())
                elif cmd == 'tact' :
                    tc300.write(b'TACT1?\r')
                    t1=tc300.readline()
                    tc300.write(b'TACT2?\r')
                    t2=tc300.readline()
                    conn.sendall(t1+b' '+t2)
    if lts_stage is not None : lts_stage.close()
    if zaber_stage is not None : zaber_stage.close()
    if tc300 is not None : tc300.close()

def client(svr,cmd) :
    if svr is None : return 'ERROR'
    remote_socket=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    remote_socket.connect((svr,65431))
    remote_socket.send(cmd.encode())
    response = remote_socket.recv(64).decode().replace('\r\n','').replace('>','')
    remote_socket.close()
    return response

