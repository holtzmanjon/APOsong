from serial import Serial

dev=Serial('COM3',115200,timeout=1)

dev.write(b'TACT1?\r')
print(dev.readline())
dev.write(b'TSET1?\r')
print(dev.readline())
dev.write(b'TSET1=140\r')
print(dev.readline())
dev.write(b'TSET1?\r')
print(dev.readline())
