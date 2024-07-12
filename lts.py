"""
lts_pythonnet
==================

An example of using the LTS integrated stages with python via pythonnet
"""
import os
import time
import sys
import clr

clr.AddReference("C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.DeviceManagerCLI.dll")
clr.AddReference("C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.GenericMotorCLI.dll")
clr.AddReference("C:\\Program Files\\Thorlabs\\Kinesis\\ThorLabs.MotionControl.IntegratedStepperMotorsCLI.dll")
from Thorlabs.MotionControl.DeviceManagerCLI import *
from Thorlabs.MotionControl.GenericMotorCLI import *
from Thorlabs.MotionControl.IntegratedStepperMotorsCLI import *
from System import Decimal  # necessary for real world units

class ThorlabsStage() :

    def __init__(self,serial_no="45441684"):

        DeviceManagerCLI.BuildDeviceList()
        # create new device
        # Connect, begin polling, and enable
        self.device = LongTravelStage.CreateLongTravelStage(serial_no)
        self.device.Connect(serial_no)

        # Ensure that the device settings have been initialized
        if not self.device.IsSettingsInitialized():
            self.device.WaitForSettingsInitialized(10000)  # 10 second timeout
            assert self.device.IsSettingsInitialized() is True

        # Start polling and enable
        self.device.StartPolling(250)  #250ms polling rate
        time.sleep(25)
        self.device.EnableDevice()
        time.sleep(0.25)  # Wait for device to enable

        # Get Device Information and display description
        self.device_info = self.device.GetDeviceInfo()
        print(self.device_info.Description)

        # Load any configuration settings needed by the controller/stage
        motor_config = self.device.LoadMotorConfiguration(serial_no)

        # Get parameters related to homing/zeroing/other
        home_params = self.device.GetHomingParams()
        print(f'Homing velocity: {home_params.Velocity}\n,'
              f'Homing Direction: {home_params.Direction}')
        home_params.Velocity = Decimal(10.0)  # real units, mm/s
        # Set homing params (if changed)
        self.device.SetHomingParams(home_params)

    def home(self) :
        # Home or Zero the device (if a motor/piezo)
        print("Homing Device")
        self.device.Home(60000)  # 60 second timeout
        print("Done")

    def get_velocity(self) :
        return Decimal.ToDouble(self.device.get_Velocity())

    def set_velocity(self,velocity) :
        # Get Velocity Params
        vel_params = self.device.GetVelocityParams()
        vel_params.MaxVelocity = Decimal(50.0)  # This is a bad idea
        self.device.SetVelocityParams(vel_params)

    def move(self,position) :
        # Move the device to a new position
        new_pos = Decimal(position)  # Must be a .NET decimal
        print(f'Moving to {new_pos}')
        self.device.MoveTo(new_pos, 60000)  # 60 second timeout
        print("Done")

    def get_position(self) :
        return Decimal.ToDouble(self.device.get_Position())

    def close(self) :
        # Stop Polling and Disconnect
        self.device.StopPolling()
        self.device.Disconnect()

