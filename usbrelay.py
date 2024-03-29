# class to control USB device, specifically Noyito relay

import pywinusb.hid as hid
from time import sleep

USB_CFG_VENDOR_ID = 0x16c0  # Should suit, if not check ID with a tool like USBDeview
USB_CFG_DEVICE_ID = 0x05DF  # Should suit, if not check ID with a tool like USBDeview


class USBRelay() :
    def __init__(self,vendor=USB_CFG_VENDOR_ID,device=USB_CFG_DEVICE_ID) :
        """ Initialize USB device
        """
        filter = hid.HidDeviceFilter(vendor_id=USB_CFG_VENDOR_ID, product_id=USB_CFG_DEVICE_ID)
        hid_device = filter.get_devices()
        self.device = hid_device[0]
        self.open_device()
        
    def open_device(self):
        """ Open device
        """
        if self.device.is_active():
            if not self.device.is_opened():
                self.device.open()
                self.get_report()
            else:
                print("Device already opened")
        else:
            print("Device is not active")

    def close_device():
        """ Close device
        """
        if self.device.is_active():
            if self.device.is_opened():
                self.device.close()
            else:
                print("Device already closed")
        else:
            print("Device is not active")

    def get_report(self):
        """ Get 'report'
        """
        if not self.device.is_active():
            self.report = None

        for rep in self.device.find_output_reports() + self.device.find_feature_reports():
            self.report = rep

    def read_status_row(self):
        """ Read current status buffer
        """
        if self.report is None:
            print("Cannot read report")
            self.last_row_status = [0, 1, 0, 0, 0, 0, 0, 0, 3]
        else:
            self.last_row_status = self.report.get()
        return self.last_row_status

    def write_row_data(self,buffer):
        """ Write new data buffer
        """
        if self.report is not None:
            self.report.send(raw_data=buffer)
            return True
        else:
            print("Cannot write in the report. check if your device is still plugged")
            return False

    def on_all(self):
        """ Turn all relays on
        """
        if self.write_row_data(buffer=[0, 0xFE, 0, 0, 0, 0, 0, 0, 1]):
            return self.read_relay_status(relay_number=3)
        else:
            print("Cannot put ON relays")

    def off_all(self):
        """ Turn all relays off
        """
        if self.write_row_data(buffer=[0, 0xFC, 0, 0, 0, 0, 0, 0, 1]):
            return self.read_relay_status(relay_number=3)
        else:
            print("Cannot put OFF relays")
            return False

    def on_relay(self,relay_number):
        """ Turn specified relay on
        """
        if self.write_row_data(buffer=[0, 0xFF, relay_number, 0, 0, 0, 0, 0, 1]):
            return self.read_relay_status(relay_number)
        else:
            print("Cannot put ON relay number {}".format(relay_number))

    def off_relay(self,relay_number):
        """ Turn specified relay off
        """
        if self.write_row_data(buffer=[0, 0xFD, relay_number, 0, 0, 0, 0, 0, 1]):
            return self.read_relay_status(relay_number)
        else:
            print("Cannot put OFF relay number {}".format(relay_number))

    def read_relay_status(self,relay_number):
        """ Get relay status for specific relaty
        """
        buffer = self.read_status_row()
        return relay_number & buffer[8]

    def is_relay_on(self,relay_number):
        """ Return if specified relay is on
        """
        return self.read_relay_status(relay_number) > 0

