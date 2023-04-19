#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import numpy as np
from nomad.metainfo import MSection, Quantity, SubSection


class MeasurementStateMachine(MSection):
    '''Representation of an instrument'''
    use_instruments = Quantity(
        type=bool,
        description='Do you have instruments connected to the computer? True for measurement setups and false for office computers.')
    system_model = Quantity(
        type=str,
        description='What is the Model number of your used Hall Measurement system? 0=[75XX-LVWR(-HS)] 1=[75XX-LVWR-SWT(-HS), 77XX-LVWR] 2=[75XX-HVWR(-HS), 77XX-HVWR] 3=[76XX] 4=[95XX-LVWR(-HS)] 5=[95XX-LVWR-SWT(-HS), 97XX-LVWR] 6=[95XX-HVWR(-HS), 97XX-HVWR] 7=[77XXA] 8=[97XXA]')
    # switch_matrix = Quantity(
    #    type=str,
    #    description='What is the Model number of your used Hall matrix? 0=[Lake Shore 775] 1=[Lake Shore 776]')
    wiring = Quantity(
        type=str,
        description='What is the wiring setup,either 0=[van der Pauw] or 1=[Hall Bar]?')
    reconfirm_ammeter_removal = Quantity(
        type=bool,
        description='Is the safety question on for Amperemeter removal above 2mA?')
    electrometer = Quantity(
        type=str,
        description='If there is a control electrometer installed, which one is it? 0=[Keithley 6512] 1=[Keithley 6514]')
    control_electrometer = Quantity(
        type=bool,
        description='Is a control electrometer installed in the setup?')
    numberofsamples = Quantity(
        type=int,
        description='How many samples are measured at the same time? Only applicable for certain Systems. Options: 1, 2 or 4')
    voltmeter = Quantity(
        type=str,
        description='Which voltmeter is used in your setup? 0=[Keithley 2000] 1=[Keithley 2182] 2=[Keithley 182]')
    currentmeter = Quantity(
        type=str,
        description='Which amperemeter is used in your setup? 0=[Keithley 6485] 1=[Keithley 485/6/7]')
    current_source = Quantity(
        type=str,
        description='Which current source is used in your setup? 0=[Keithley 220] 1=[Keithley 6220]')
    ac_hall = Quantity(
        type=bool,
        description='Does your setup use AC Hall measurements.')
    ac_hall_type = Quantity(
        type=str,
        description='If your setup uses AC Hall measurements, does it use 0=[AC Current only] or 1=[AC Field and Current]?')
    contact_blasting = Quantity(
        type=bool,
        description='Is the contact blasting option turned on? Only for Keithley 2400 at IEEE address 24.')


class Keithley7001(MSection):
    '''Representation of an instrument'''
    resting_state = Quantity(
        type=str,
        description='Are all relays in the switch matrix 0=[closed (defalut)] or 1=[open] in resting state?')


class Keithley6485(MSection):
    '''Representation of an instrument'''
    reading_rate = Quantity(
        type=str,
        description='What is the reading rate of the current meter? 0=[Slow (default)] 1=[Medium] 2=[Fast]')
    median_filter = Quantity(
        type=bool,
        description='Is a median filter turned on? False=[off (default)]')
    median_filter_rank = Quantity(
        type=np.dtype(np.float64),
        description='If the median filter is turned on, what is its rank? Value must be between 1 and 5.')
    digital_filter = Quantity(
        type=bool,
        description='Is a digitial filter turned on? False=[off (default)]')
    digital_filter_type = Quantity(
        type=str,
        description='If the digital filter is turned on, what is its mode? 0=[Repeat] 1=[Moving]')
    digital_filter_count = Quantity(
        type=int,
        description='If the digital filter is turned on, what is its count? Value must be between 1 and 100.')
    advanced_filter = Quantity(
        type=bool,
        description='If the digital filter is turned on, do you use and advanced filter? False=[off (default)]')
    noise_window = Quantity(
        type=np.dtype(np.float64),
        description='What is the cut off noise window in % of the full scale for the advanced filter? Value must be between 0 and 105.')


class Keithley220(MSection):
    '''Representation of an instrument'''
    compliance_voltage = Quantity(
        type=np.dtype(np.float64),
        description='What is the voltage compliance for your current source? in V. (default = 100 V)',
        unit="volt")


class Keithley2000(MSection):
    '''Representation of an instrument'''
    reading_rate = Quantity(
        type=str,
        description='What is the set reading rate for the Keithley2000? 0=[Slow] 1=[Medium(default)] 2=[Fast]')
    digital_filter = Quantity(
        type=bool,
        description='Is a digitial filter turned on? True=[on (default)]')
    digital_filter_type = Quantity(
        type=str,
        description='If the digital filter is turned on, what is its mode? 0=[Repeat] 1=[Moving]')
    digital_filter_count = Quantity(
        type=int,
        description='If the digital filter is turned on, what is its count? Value must be between 1 and 100.')


class Keithley2182(MSection):
    '''Representation of an instrument'''
    line_synchronization = Quantity(
        type=bool,
        description='Is line synchronization turned on? True=[on (default)]')
    reading_rate = Quantity(
        type=str,
        description='What is the set reading rate for the Keithley2182? 0=[Slow(default)] 1=[Medium] 2=[Fast]')
    analog_filter = Quantity(
        type=bool,
        description='Is the analog filter turned on? False=[off (default)]')
    digital_filter = Quantity(
        type=bool,
        description='Is a digitial filter turned on? True=[on (default)]')
    digital_filter_type = Quantity(
        type=str,
        description='If the digital filter is turned on, what is its mode? 0=[Repeat] 1=[Moving]')
    digital_filter_window = Quantity(
        type=str,
        description='If the digital filter is turned on, what is its window? 0=[None] 1=[0.01%] 2=[0.1%] 3=[1%] 4=[10%]')
    digital_filter_count = Quantity(
        type=int,
        description='If the digital filter is turned on, what is its count? Value must be between 1 and 100.')


class Keithley182(MSection):
    '''Representation of an instrument'''
    integration_time = Quantity(
        type=str,
        description='What is the set integration time for the Keithley182? 0=[1 line cycle (default)] 1=[3 ms] 2=[100 ms]')
    analog_filter = Quantity(
        type=bool,
        description='Is the analog filter turned on? False=[off (default)]')
    digital_filter = Quantity(
        type=str,
        description='Is the digital filter turned on, what is it set on? 0=[Off] 1=[Fast] 2=[Medium (default)] 3=[Slow]')


class Keithley2700(MSection):
    '''Representation of an instrument'''
    line_synchronization = Quantity(
        type=bool,
        description='Is line synchronization turned on? True=[on (default)]')
    reading_rate = Quantity(
        type=str,
        description='What is the set reading rate for the Keithley2700? 0=[Slow] 1=[Medium (default)] 2=[Fast]')
    resting_state = Quantity(
        type=str,
        description='Are all relays in the switch matrix 0=[closed (defalut)] or 1=[open] in resting state?')
    digital_filter = Quantity(
        type=bool,
        description='Is a digitial filter turned on? True=[on (default)]')
    digital_filter_type = Quantity(
        type=str,
        description='If the digital filter is turned on, what is its mode? 0=[Repeat] 1=[Moving]')
    digital_filter_window = Quantity(
        type=str,
        description='If the digital filter is turned on, what is its window? 0=[None] 1=[0.01%] 2=[0.1%] 3=[1%] 4=[10%]')
    digital_filter_count = Quantity(
        type=int,
        description='If the digital filter is turned on, what is its count? Value must be between 1 and 100.')


class Keithley2400(MSection):
    '''Representation of an instrument'''
    compliance_voltage = Quantity(
        type=np.dtype(np.float64),
        description='What is the voltage compliance for your current source? in V? Must be a value between 0 V and 210 V. (default = 20 V) Value of 20 V and lower will set current compliance to 1 A. Value of larger than 20 V will set current compliance to 100 mA.')


class FieldController(MSection):
    '''Representation of an instrument'''
    gaussmeter = Quantity(
        type=str,
        description='Which gaussmeter is used in your setup? 0=[LS450] 1=[LS475] 2=[LS736]')
    use_ls475_built_in_pi_controller = Quantity(
        type=bool,
        description='If you use LS475, do you also use the built-in PI Control? True=[yes (default)]')
    setpoint_ramping = Quantity(
        type=bool,
        description='If you use LS475 with the built-in PI Control, do you use setpoint ramping? False=[no (default)]')
    ramp_rate = Quantity(
        type=np.dtype(np.float64),
        description='If you use LS475 with the built-in PI Control and setpoint ramping, what is the ramp rate in gauss/s (default = 500 gauss/s).')
    field_p = Quantity(
        type=np.dtype(np.float64),
        description='What is the loop gain P in the loop control (PID) of the field? (default=1)')
    field_i = Quantity(
        type=np.dtype(np.float64),
        description='What is the reciprocal loop reset time I in the loop control (PID) of the field? (default=25)')
    field_d = Quantity(
        type=np.dtype(np.float64),
        description='What is the loop rate time D in the loop control (PID) of the field? (default=0)')
    power_supply = Quantity(
        type=str,
        description='Which power supply is used in your setup? 0=[642,643, or 662 (Max. 70 Amp)] 1=[665 (Max. 100 Amp)] 2=[648 or 668 (Max. 135 Amp)] 3=[647 (Max. 72 Amp Must be controlled via EXT input terminal)]')
    settle_band = Quantity(
        type=np.dtype(np.float64),
        description='What is the set settle band in gauss [G]? (default=2)')
    settle_time = Quantity(
        type=np.dtype(np.float64),
        description='What is the set settle time in seconds [sec]? (default=1)')
    ls450_display_unit = Quantity(
        type=str,
        description='Which front panel field display unit do you use on your LS450? 0=[Gauss] 1=[Tesla]')
    shutdown_field_at_end_of_measurement = Quantity(
        type=bool,
        description='Do you turn off the field at the end of the measurement? False=[no]')
    log_field_data = Quantity(
        type=bool,
        description='Do you log status/debug data to "fieldcontrollerlog.txt file? (Do not check. This option is for diagnostic purposes only.) False=[no (default)]')
    enable_settle_time_out = Quantity(
        type=bool,
        description='Do you enable a time out time for the case that the field does not settle? False=[no (default)]')
    settle_time_out = Quantity(
        type=np.dtype(np.float64),
        description='If you enable a settle time out, what is the time out period in minutes?')
    time_to_setpoint_for_demo = Quantity(
        type=np.dtype(np.float64),
        description='What is the maximum time to setpoint in minutes?')


class TemperatureController(MSection):
    '''Representation of an instrument'''
    sample_sensor_type = Quantity(
        type=str,
        description='What is the sample sensor type installed on your measuring rod? 0=[Silicon diode] 1=[GaAlAs diode] 2=[Platinum 100/250] 3=[Platinum 100/500] 4=[Platinum 1000] 5=[Rhodium Iron] 6=[Carbon-glass] 7=[Cernox] 8=[RuOx] 9=[Germanium] 10=[Thermocouple]')
    sample_sensor_curve_number = Quantity(
        type=np.dtype(np.float64),
        description='What is the calibration curve number of your used sensor?')
    sample_sensor_filter = Quantity(
        type=bool,
        description='Do you use a filter for the sample temperature sensor? True=[yes/on (default)]')
    sample_sensor_filter_points = Quantity(
        type=np.dtype(np.float64),
        description='If you use a filter for the sample temperature sensor, how many sampling points, do you use? Value must be between 2 and 64. (default=10).')
    sample_sensor_filter_window = Quantity(
        type=np.dtype(np.float64),
        description='If you use a filter for the sample temperature sensor, how large is the sampling window in [%]? Enter a value between 1 and 10.')
    sample_sensor_room_compensation = Quantity(
        type=bool,
        description='Do you use the sample temperature sensor room compensation? True=[yes (default)]')


class TemperatureDomain(MSection):
    '''Representation of an instrument'''
    temperature_low = Quantity(
        type=np.dtype(np.float64),
        unit="kelvin",
        # a_eln=dict(component='NumberEditQuantity', defaultDisplayUnit='celsius'),
        description='The low temperature point of this domain.')
    temperature_high = Quantity(
        type=np.dtype(np.float64),
        unit="kelvin",
        # a_eln=dict(component='NumberEditQuantity', defaultDisplayUnit='celsius'),
        description='The high temperature point of this domain.')
    # direction = Quantity(
    #     type=str,
    #     description='What is the temperature loop direction for this temperature domain? 0=[Both] 1=[Ascending] 2=[Descending]')
    # TODO it doesn't work
    loop_1_p = Quantity(
        type=np.dtype(np.float64),
        # unit="celsius",
        description='What is the loop gain P in the loop 1 control parameters (PID) of the temperature?')
    loop_2_p = Quantity(
        type=np.dtype(np.float64),
        # unit="celsius",
        description='What is the loop gain P in the loop 2 control parameters (PID) of the temperature')
    loop_1_i = Quantity(
        type=np.dtype(np.float64),
        # unit="celsius",
        description='What is the reciprocal loop reset time I in the loop 1 control parameters (PID) of temperature?')
    loop_2_i = Quantity(
        type=np.dtype(np.float64),
        # unit="celsius",
        description='What is the reciprocal loop reset time I in the loop 2 control parameters (PID) of temperature?')
    loop_1_d = Quantity(
        type=np.dtype(np.float64),
        # unit="celsius",
        description='What is the loop rate time D in the loop 1 control parameters (PID) of the temperature?')
    loop_2_d = Quantity(
        type=np.dtype(np.float64),
        # unit="celsius",
        description='What is the loop rate time D in the loop 2 control parameters (PID) of the temperature?')
    heater_range = Quantity(
        type=str,
        description='What is the heater range for this temperature domain? 0=[heater off] 1=[4.9mW] 2=[49mW] 3=[490mW] 4=[4.9W] 5=[49.0W]')
    loop_2_heater = Quantity(
        type=str,
        # unit="celsius",
        description='Is the heater for the loop 2 controller turned on? 0=[off] 1=[on]')
    loop_1_setpoint_offset = Quantity(
        type=np.dtype(np.float64),
        unit="kelvin",
        description='What is the setpoint offset of loop controller 1 in Kelvin?')
    ramp_rate = Quantity(
        type=np.dtype(np.float64),
        unit="kelvin / minute",
        description='What is the setpoint offset of loop controller 1 in Kelvin?')
    settle_wait_time = Quantity(
        type=np.dtype(np.float64),
        unit="minute",
        description='What is the wait time after reaching the settle point in minutes?')
    settle_temperature_drift = Quantity(
        type=np.dtype(np.float64),
        unit="kelvin / minute",
        description='What is the maximum drift of the temperature to settle in Kelvin/minutes?')
    settle_temperature_band = Quantity(
        type=np.dtype(np.float64),
        unit="kelvin",
        description='What is the maximum temperature deviation to settle in Kelvin?')
    settle_time_out = Quantity(
        type=bool,
        description='Is a time out for settling enabled? True=[yes (default)]')
    settle_time_out_period = Quantity(
        type=np.dtype(np.float64),
        unit="minute",
        description='If a time out is enabled, what is the the settle time out period in minutes?')
    needle_valve_during_ramp = Quantity(
        type=np.dtype(np.float64),
        # unit="celsius",
        description='Needle valve position during temperature ramp.')
    needle_valve_after_ramp = Quantity(
        type=np.dtype(np.float64),
        # unit="celsius",
        description='Needle valve position after temperature ramp.')
    sample_space_evacuation_valve_during_ramp = Quantity(
        type=str,
        description='Is the sample space evacuation valve 0=[open] or 1=[closed] during the temperature ramp?')
    sample_space_evacuation_valve_after_ramp = Quantity(
        type=str,
        description='Is the sample space evacuation valve 0=[open] or 1=[closed] after the temperature ramp?')


class Instrument(MSection):
    '''Representation of an instrument'''
    hms_software_version = Quantity(
        type=str,
        description='HMS software version')
    working_directory = Quantity(
        type=str,
        description='Working Directory')
    use_instruments = Quantity(
        type=bool,
        description='Do you have instruments connected to the computer? True for measurement setups and false for office computers.')
    system_model = Quantity(
        type=str,
        description='What is the Model number of your used Hall Measurement system? 0=[75XX-LVWR(-HS)] 1=[75XX-LVWR-SWT(-HS), 77XX-LVWR] 2=[75XX-HVWR(-HS), 77XX-HVWR] 3=[76XX] 4=[95XX-LVWR(-HS)] 5=[95XX-LVWR-SWT(-HS), 97XX-LVWR] 6=[95XX-HVWR(-HS), 97XX-HVWR] 7=[77XXA] 8=[97XXA]')
    wiring = Quantity(
        type=str,
        description='What is the wiring setup,either 0=[van der Pauw] or 1=[Hall Bar]?')
    reconfirm_ammeter_removal = Quantity(
        type=bool,
        description='Is the safety question on for Amperemeter removal above 2mA?')
    electrometer = Quantity(
        type=str,
        description='If there is a control electrometer installed, which one is it? 0=[Keithley 6512] 1=[Keithley 6514]')
    control_electrometer = Quantity(
        type=bool,
        description='Is a control electrometer installed in the setup?')
    numberofsamples = Quantity(
        type=int,
        description='How many samples are measured at the same time? Only applicable for certain Systems. Options: 1, 2 or 4')
    voltmeter = Quantity(
        type=str,
        description='Which voltmeter is used in your setup? 0=[Keithley 2000] 1=[Keithley 2182] 2=[Keithley 182]')
    currentmeter = Quantity(
        type=str,
        description='Which amperemeter is used in your setup? 0=[Keithley 6485] 1=[Keithley 485/6/7]')
    current_source = Quantity(
        type=str,
        description='Which current source is used in your setup? 0=[Keithley 220] 1=[Keithley 6220]')
    ac_hall = Quantity(
        type=bool,
        description='Does your setup use AC Hall measurements.')
    ac_hall_type = Quantity(
        type=str,
        description='If your setup uses AC Hall measurements, does it use 0=[AC Current only] or 1=[AC Field and Current]?')
    contact_blasting = Quantity(
        type=bool,
        description='Is the contact blasting option turned on? Only for Keithley 2400 at IEEE address 24.')

    keithley_182 = SubSection(section_def=Keithley182)
    keithley_220 = SubSection(section_def=Keithley220)
    keithley_2000 = SubSection(section_def=Keithley2000)
    keithley_2182 = SubSection(section_def=Keithley2182)
    keithley_2400 = SubSection(section_def=Keithley2400)
    keithley_2700 = SubSection(section_def=Keithley2700)
    keithley_6485 = SubSection(section_def=Keithley6485)
    keithley_7001 = SubSection(section_def=Keithley7001)
    field_controller = SubSection(section_def=FieldController)
    temperature_controller = SubSection(section_def=TemperatureController)
    temperature_domain = SubSection(section_def=TemperatureDomain, repeats=True)
