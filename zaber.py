from zaber_motion import Units
from zaber_motion.ascii import Connection

class ZaberStage() :

    def __init__(self,port='COM4') :
        self.connection = Connection.open_serial_port(port)
        self.connection.enable_alerts()

        device_list = self.connection.detect_devices()
        print("Found {} devices".format(len(device_list)))

        self.device = device_list[0]
        self.axis = self.device.get_axis(1)
        if not self.axis.is_homed():
            print('homing axis...')
            self.axis.home()

    def home(self) :
        self.axis.home()

    def move(self,pos,units=Units.LENGTH_MILLIMETRES,absolute=True) :
        if absolute : 
            self.axis.move_absolute(pos, Units.LENGTH_MILLIMETRES)
        else :
            axis.move_relative(5, Units.LENGTH_MILLIMETRES)

    def get_position(self) :
        return self.axis.get_position(Units.LENGTH_MILLIMETRES)

    def close(self) :
        self.connection.close()
