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
from nomad.metainfo import MSection, Quantity, SubSection, Datetime, Section


class IVData(MSection):
    '''Container for IV-Curve measurement data'''
    m_def = Section(
        a_eln=dict(lane_width='600px'),
        a_plot=[
            {
                'label': 'IV-Curve',
                'x': 'current',
                'y': [
                    'voltage',
                    'best_fit_values'
                ],
                'layout': {'yaxis': {'type': 'lin'}},
                'lines': [{'mode': 'markers'}],
            }
        ]
    )

    contact_set = Quantity(
        type=str,
        description='The contact set being used'
    )
    best_fit_resistance = Quantity(
        type=np.float64,
        unit='ohm',
        description='The resistance value from the IV-curve best fit.'
    )
    best_fit_offset = Quantity(
        type=np.float64,
        unit='volt',
        description='The offset of the IV-curve best fit.'
    )
    correlation = Quantity(
        type=np.float64,
        description='The best fit correlation.'
    )
    best_fit_values = Quantity(
        type=np.float64,
        unit='volt',
        shape=['*'],
        description='The discrete values of the best fit.'
    )
    current = Quantity(
        type=np.dtype(np.float64),
        unit='ampere',
        shape=['*'],
        description='The applied current steps.',
    )
    voltage = Quantity(
        type=np.dtype(np.float64),
        unit='volt',
        shape=['*'],
        description='The measured voltages.')
    field = Quantity(
        type=np.dtype(np.float64),
        unit='gauss',
        shape=['*'],
        description='The applied magnetic field steps of the experiment.')
    temperature = Quantity(
        type=np.dtype(np.float64),
        unit='kelvin',
        shape=['*'],
        description='The temperatures values of the experiment.')


class VariableFieldData(MSection):
    '''Container for variable magnetic field data'''
    m_def = Section(
        a_eln=dict(lane_width='600px'),
        a_plot=[
            {
                'label': 'Carrier density',
                'x': 'field',
                'y': 'carrier_density',
                'layout': {'yaxis': {'type': 'lin'}},
                'lines': [{'mode': 'markers'}],
            },
            {
                'label': 'Hall coefficient',
                'x': 'field',
                'y': 'hall_coefficient',
                'layout': {'yaxis': {'type': 'lin'}},
                'lines': [{'mode': 'markers'}],
            },
            {
                'label': 'Hall mobility',
                'x': 'field',
                'y': 'hall_mobility',
                'layout': {'yaxis': {'type': 'lin'}},
                'lines': [{'mode': 'markers'}],
            },
            {
                'label': 'Resistivity',
                'x': 'field',
                'y': 'resistivity',
                'layout': {'yaxis': {'type': 'lin'}},
                'lines': [{'mode': 'markers'}],
            }
        ]
    )

    field = Quantity(
        type=np.dtype(np.float64),
        unit='gauss',
        shape=['*'],
        description='The applied magnetic field steps of the experiment.'
    )
    resistivity = Quantity(
        type=np.dtype(np.float64),
        unit='ohm * cm',
        shape=['*'],
        description='The measured resistivities.'
    )
    hall_coefficient = Quantity(
        type=np.dtype(np.float64),
        unit='cm**3 / C',
        shape=['*'],
        description='The measured hall coefficients.'
    )
    carrier_density = Quantity(
        type=np.dtype(np.float64),
        unit='1 / cm**3',
        shape=['*'],
        description='The measured carrier densities.'
    )
    hall_mobility = Quantity(
        type=np.dtype(np.float64),
        unit='cm**2 / volt / second',
        shape=['*'],
        description='The measured hall mobilities.'
    )
    temperature = Quantity(
        type=np.dtype(np.float64),
        unit='kelvin',
        shape=['*'],
        description='The temperature steps of the experiment.'
    )


class VariableTemperatureData(MSection):
    '''Container for variable hall temperature data'''
    m_def = Section(
        a_eln=dict(lane_width='600px'),
        a_plot=[
            {
                'label': 'Carrier density',
                'x': 'temperature',
                'y': 'carrier_density',
                'layout': {'yaxis': {'type': 'lin'}},
                'lines': [{'mode': 'markers'}],
            },
            {
                'label': 'Hall coefficient',
                'x': 'temperature',
                'y': 'hall_coefficient',
                'layout': {'yaxis': {'type': 'lin'}},
                'lines': [{'mode': 'markers'}],
            },
            {
                'label': 'Hall mobility',
                'x': 'temperature',
                'y': 'hall_mobility',
                'layout': {'yaxis': {'type': 'lin'}},
                'lines': [{'mode': 'markers'}],
            },
            {
                'label': 'Resistivity',
                'x': 'temperature',
                'y': 'resistivity',
                'layout': {'yaxis': {'type': 'lin'}},
                'lines': [{'mode': 'markers'}],
            }
        ]
    )

    carrier_density = Quantity(
        type=np.dtype(np.float64),
        unit='1 / centimeter ** 3',
        shape=['*'],
        description='The measured carrier density.',
    )
    field = Quantity(
        type=np.dtype(np.float64),
        unit='gauss',
        shape=['*'],
        description='The magnetic field steps of the experiment.')
    hall_coefficient = Quantity(
        type=np.dtype(np.float64),
        unit='centimeter ** 3 / coulomb',
        shape=['*'],
        description='The measured hall coefficients.')
    hall_mobility = Quantity(
        type=np.dtype(np.float64),
        unit='centimeter ** 2 / volt / second',
        shape=['*'],
        description='The measured hall mobilities.')
    resistivity = Quantity(
        type=np.dtype(np.float64),
        unit='centimeter * ohm',
        shape=['*'],
        description='The measured resistivity values.')
    temperature = Quantity(
        type=np.dtype(np.float64),
        unit='kelvin',
        shape=['*'],
        description='The temperature steps of the experiment.')
    zero_field_resistivity = Quantity(
        type=np.dtype(np.float64),
        unit='centimeter * ohm',
        shape=['*'],
        description='The measured zero field resistivities.')


class Measurement(MSection):
    '''A general hall measaurement representation'''
    data: SubSection


class IVCurveMeasurement(Measurement):
    '''Representation of an IV curve measurement'''
    current_step = Quantity(
        type=np.dtype(np.float64),
        unit='nA',
        description='The current steps used in this measurement.')
    dwell_time = Quantity(
        type=np.dtype(np.float64),
        unit='seconds',
        description='The dwell time of the measurement.')
    elapsed_time = Quantity(
        type=np.dtype(np.float64),
        unit='seconds',
        description='The duration of the experiment.')
    starting_current = Quantity(
        type=np.dtype(np.float64),
        unit='uA',
        description='The start current of the experiment.')
    ending_current = Quantity(
        type=np.dtype(np.float64),
        unit='uA',
        description='The ending current of the measurement.')
    resistance_range = Quantity(
        type=str,
        description='The resistance range set in the experiment.')
    start_time = Quantity(
        type=Datetime,
        description='The start time of the experiment.')
    time_completed = Quantity(
        type=Datetime,
        description='The end time of the experiment.')
    contact_sets = SubSection(section_def=IVData, repeats=True)


class VariableFieldMeasurement(Measurement):
    '''Representation of an variable magnetic field hall measurement'''
    start_time = Quantity(
        type=Datetime,
        description='The start time of the measurement.'
    )
    time_completed = Quantity(
        type=Datetime,
        description='The end time of the measurement.'
    )
    elapsed_time = Quantity(
        type=np.dtype(np.float64),
        unit='seconds',
        description='The duration of the experiment.'
    )
    field_profile = Quantity(
        type=str,
        description='The magnetic field profile.'
    )
    maximum_field = Quantity(
        type=np.dtype(np.float64),
        unit='kG',
        description='The maximum magnetic field during the measurement.'
    )
    minimum_field = Quantity(
        type=np.dtype(np.float64),
        unit='kG',
        description='The minimum magnetic field during the measurement.'
    )
    field_step = Quantity(
        type=np.dtype(np.float64),
        unit='kG',
        description='The magnetic field steps of the measurement.'
    )
    direction = Quantity(
        type=str,
        description='The direction of the magnetic field scan.'
    )
    measurement_type = Quantity(
        type=str,
        description='The measurement type.'
    )
    excitation_current = Quantity(
        type=np.dtype(np.float64),
        unit='uA',
        description='The excitation current of the measurement.'
    )
    resistance_range = Quantity(
        type=str,
        description='The resistiance range during the measurement.'
    )
    dwell_time = Quantity(
        type=np.dtype(np.float64),
        unit='seconds',
        description='The dwell time during the measurement.'
    )
    current_reversal = Quantity(
        type=bool,
        description='Current reversal setting of the measurement.'
    )
    geometry_selection = Quantity(
        type=str,
        description='The selected measurement geometry.'
    )
    use_zero_field_resistivity = Quantity(
        type=bool,
        description='Has the zero field resistivity been used for the measurement?'
    )
    zero_field_resistivity = Quantity(
        type=np.dtype(np.float64),
        unit='ohm * centimeter',
        description='The measured zero-field resistivity.'
    )
    field_at_zero_resistivity = Quantity(
        type=np.dtype(np.float64),
        unit='gauss',
        description='The magnetic field while measuring zero-field resistivity.'
    )
    temperature_at_zero_resistivity = Quantity(
        type=np.dtype(np.float64),
        unit='kelvin',
        description='The temperature while measuring zero-field resistivity.'
    )

    data = SubSection(section_def=VariableFieldData, repeats=True)


class VariableTemperatureMeasurement(Measurement):
    '''Representation of a variable temperature hall measurement'''
    current_reversal = Quantity(
        type=bool,
        description='Was the current reversed during the measurement?')
    dwell_time = Quantity(
        type=np.dtype(np.float64),
        unit='seconds',
        description='The dwell time of the measurement.')
    elapsed_time = Quantity(  # This is already present in base_classes.Activity
        type=np.dtype(np.float64),
        unit='seconds',
        description='The duration of the measurement.')
    ending_temperature = Quantity(
        type=np.dtype(np.float64),
        unit='kelvin',
        description='The ending temperature of the measurement.')
    excitation_current = Quantity(
        type=np.dtype(np.float64),
        unit='nA',
        description='The excitation current during the measurement.')
    field_at = Quantity(
        type=np.dtype(np.float64),
        unit='kG',
        description='The magnetic field during the measurement.')
    geometry_selection = Quantity(
        type=str,
        description='The selected measurement geometry.')
    measurement_type = Quantity(  # This is already present in base_classes.Activity
        type=str,
        description='The selected measurement type.')
    resistance_range = Quantity(
        type=str,
        description='The selected resistance range.')
    spacing = Quantity(
        type=str,
        description='The selected spacing.')
    start_time = Quantity(  # This is already present in base_classes.Activity
        type=Datetime,
        description='The start time of the measurement.')
    starting_temperature = Quantity(
        type=np.dtype(np.float64),
        unit='kelvin',
        description='The starting temperature of the measurement.')
    temperature_step = Quantity(
        type=np.dtype(np.float64),
        unit='kelvin',
        description='The temperature step of the measurement.')
    time_completed = Quantity(  # This is already present in base_classes.Activity
        type=Datetime,
        description='The end time of the measurement.')

    data = SubSection(section_def=VariableTemperatureData, repeats=True)
