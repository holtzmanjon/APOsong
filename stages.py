import lts
import zaber

HOST = "10.75.0.22"  # Standard loopback interface address (localhost)
PORT = 65432  # Port to listen on (non-privileged ports are > 1023)

def remote() :
  """ Run simple remote socket server 
  """
  lts=lts.ThorlabsStage()
  zaber=zaber.ZaberStage()
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



