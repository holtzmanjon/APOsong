import requests
import subprocess
import os

def send(recipients,subject=None,message=None,attachment=None,snapshot=False,html=False) :
    """ Send email to recipients
    """
    f=open('message','w')
    if html : f.write('<HTML><BODY>')
    f.write(message)
    if html : f.write('</BODY></HTML>')
    f.close()
    f=open('message')
    cmd = ['mail']
    if html : cmd.extend(['-M','text/html'])
    if subject is not None : cmd.extend(['-s',subject])
    if snapshot : 
        auth = requests.auth.HTTPDigestAuth('snapshot',os.environ['VIDEOPASS'])
        response=requests.get("http://video1m.apo.nmsu.edu/cgi-bin/snapshot.cgi",stream=True,auth=auth)
        try : os.remove('webcam_snapshot.jpg')
        except : pass
        fimg=open('webcam_snapshot.jpg','wb')
        fimg.write(response.content)
        fimg.close()
        cmd.extend(['-a','webcam_snapshot.jpg'])
    if attachment is not None : 
        if isinstance(attachment,list) :
            for a in attachment :
                cmd.extend(['-a',a])
        else :
            cmd.extend(['-a',attachment])
    cmd.extend(recipients)
    subprocess.run(cmd,stdin=f)
    f.close()
    try : os.remove('webcam_snapshot.jpg')
    except : pass

