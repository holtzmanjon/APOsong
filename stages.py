import lts
import zaber
imnport socket

HOST = "10.75.0.202"  # Standard loopback interface address (localhost)
HOST = "172.24.4.202"  # Standard loopback interface address (localhost)
PORT = 76543  # Port to listen on (non-privileged ports are > 1023)

def remote() :
  """ Run simple remote socket server 
  """
  s1=lts.ThorlabsStage()
  s2=zaber.ZaberStage()
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
            cmd = data.decode()
            if cmd == 'iodine 1' :
                s1.move(1000)
            elif cmd == 'iodine 0' :
                s1.move(0)
            elif cmd.split()[0] == 'focus' :
                s2.move(int(cmd.split()[1]))




