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
from typing import Generator, Tuple, Optional
import re

from .measurement import (
    Measurement,
    VariableTemperatureMeasurement,
    VariableTemperatureData,
    VariableFieldMeasurement,
    VariableFieldData,
    IVCurveMeasurement,
    IVData
)

from .hall_instrument import (
    Instrument,
    Keithley7001,
    Keithley6485,
    Keithley182,
    Keithley2000,
    Keithley2182,
    Keithley220,
    Keithley2400,
    Keithley2700,
    FieldController,
    TemperatureController,
    TemperatureDomain
)

import nexusutils.dataconverter.readers.hall.helpers as hall_helpers  # pylint: disable=import-error
from nomad.units import ureg


def get_measurement_object(measurement_type: str) -> Measurement:
    '''
    Gets a measurement MSection object from the given measurement type.

    Args:
        measurement_type (str): The measurement type.

    Returns:
        Measurement: A MSection representing a Hall measurement.
    '''
    if measurement_type == 'Variable Temperature Measurement':
        return VariableTemperatureMeasurement()
    if measurement_type == 'Variable Field Measurement':
        return VariableFieldMeasurement()
    if measurement_type == 'IV Curve Measurement':
        return IVCurveMeasurement()
    return Measurement()


def get_data_object(measurement_type: str):
    '''
    Gets a measurement data MSection object from the given measurement type.

    Args:
        measurement_type (str): The measurement type.

    Returns:
        A MSection representing a Hall measurement data object.
    '''
    if measurement_type == 'Variable Temperature Measurement':
        return VariableTemperatureData()
    if measurement_type == 'Variable Field Measurement':
        return VariableFieldData()
    if measurement_type == 'IV Curve Measurement':
        return IVData()
    return None


def clean_key_str(key: str) -> str:
    '''
    Cleans a key, i.e.e replaces spaces and `-` with `_` and converts to lower case

    Args:
        key (str): The `dirty` key

    Returns:
        str: The clean key
    '''
    return re.sub(r'[\s_\-]+', '_', key).lower()


def split_value_unit(expr: str) -> Tuple[str, Optional[str]]:
    '''
    Searches for a value unit pair and returns the values for a combination
    `value [unit]`.

    Args:
        expr (str): The expression to search for a value unit pair.

    Returns:
        Tuple[str, Optional[str]]:
            A tuple of value and unit.
            Returns the expr, where spaces are replaced with `_` and None when no
            value unit expression is found.
    '''
    is_value_unit = re.search(r'([^\[]+)\s\[(.*)\]', expr)
    if is_value_unit:
        value = clean_key_str(is_value_unit.group(1))
        unit = is_value_unit.group(2)
        return value.lower(), hall_helpers.clean(unit)
    return clean_key_str(expr), None


def rename_key(key: str) -> str:
    '''
    Renames the key from the file to the eln

    Args:
        key (str): They key as read from the file.

    Returns:
        str: The key replaced with its eln counterpart
    '''
    key_map = {
        'use_zero_field_resistivity_to_calculate_hall_mobility':
            'use_zero_field_resistivity',
        'at_field': 'field_at_zero_resistivity',
        'at_temperature': 'temperature_at_zero_resistivity'
    }
    return key_map.get(key, key)


def calc_best_fit_values(iv_measurement: IVData) -> IVData:
    '''
    Calculates the best fit voltage values from the provided
    fitting data.

    Args:
        iv_measurement (IVData): The IVData without discrete best fit values.

    Returns:
        IVData: The IVdata with discret best fit values
    '''
    iv_measurement.best_fit_values = (
        iv_measurement.current * iv_measurement.best_fit_resistance
        + iv_measurement.best_fit_offset
    )

    return iv_measurement


def get_measurements(data_template: dict) -> Generator[Measurement, None, None]:
    '''
    Returns a hall measurement MSection representation form its corresponding
    nexus data_template.

    Args:
        data_template (dict): The nomad-parser-nexus data template.

    Yields:
        Generator[Measurement, None, None]:
            A generator yielding the single hall measurements.
    '''
    highest_index = 1

    for key in data_template:
        if bool(re.search(f'^/entry/measurement/{highest_index}_.+/', key)):
            highest_index += 1

    for measurement_index in range(1, highest_index):
        first = True
        data_entries: dict = {}
        contact_sets: dict = {}

        for key in data_template:
            if not key.startswith(f'/entry/measurement/{measurement_index}_'):
                continue

            if first:
                measurement_type = re.search(
                    f'measurement/{measurement_index}_([^/]+)/', key
                ).group(1)
                first = False
                eln_measurement = get_measurement_object(measurement_type)

            clean_key = clean_key_str(key.split(f'{measurement_type}/')[1])

            regexp = re.compile('/data(\\d+)/')
            if bool(regexp.search(key)):
                data_index = regexp.search(key).group(1)
                if data_index not in data_entries:
                    data_entries[data_index] = get_data_object(measurement_type)
                clean_dkey = clean_key.split(f'data{data_index}/')[1]
                if hasattr(data_entries[data_index], clean_dkey):
                    if f'{key}/@units' in data_template:
                        setattr(
                            data_entries[data_index],
                            clean_dkey,
                            data_template[key] * ureg(data_template[f'{key}/@units'])
                        )
                    else:
                        setattr(data_entries[data_index], clean_dkey, data_template[key])
                continue

            if '/Contact Sets/' in key:
                contact_set = re.search('/Contact Sets/([^/]+)/', key).group(1)
                if contact_set not in contact_sets:
                    contact_sets[contact_set] = get_data_object(measurement_type)
                clean_dkey = clean_key.split(f'{contact_set.lower()}/')[1]
                contact_sets[contact_set].contact_set = contact_set

                if 'data0' in key:
                    data = data_template[key]

                    for column in data.columns:
                        col, unit = split_value_unit(column)
                        clean_col = col.lower().replace(' ', '_')
                        if hasattr(contact_sets[contact_set], clean_col):
                            if unit is not None:
                                setattr(
                                    contact_sets[contact_set],
                                    clean_col,
                                    data[column] * ureg(unit)
                                )
                            else:
                                setattr(
                                    contact_sets[contact_set], clean_col, data[column]
                                )
                    continue

                clean_dkey, unit = split_value_unit(key.split(f'{contact_set}/')[1])
                if hasattr(contact_sets[contact_set], clean_dkey):
                    if unit is not None:
                        setattr(
                            contact_sets[contact_set],
                            clean_dkey,
                            data_template[key] * ureg(unit)
                        )
                    elif f'{key}/@units' in data_template:
                        setattr(
                            contact_sets[contact_set],
                            clean_dkey,
                            data_template[key] * ureg(data_template[f'{key}/@units'])
                        )
                    else:
                        setattr(
                            contact_sets[contact_set], clean_dkey, data_template[key]
                        )
                continue

            clean_key, unit = split_value_unit(key.split(f'{measurement_type}/')[1])
            clean_key = rename_key(clean_key)
            if hasattr(eln_measurement, clean_key):
                if f'{key}/@units' in data_template:
                    setattr(
                        eln_measurement,
                        clean_key,
                        data_template[key] * ureg(data_template[f'{key}/@units'])
                    )
                elif unit is not None:
                    setattr(
                        eln_measurement,
                        clean_key,
                        data_template[key] * ureg(unit)
                    )
                else:
                    setattr(eln_measurement, clean_key, data_template[key])

        eln_measurement.data = []
        for data_entry in data_entries.values():
            eln_measurement.data.append(data_entry)

        for data_entry in contact_sets.values():
            eln_measurement.contact_sets.append(calc_best_fit_values(data_entry))

        yield eln_measurement


def get_instrument(data_template: dict, logger):
    '''
    Returns a hall instrument MSection representation form its corresponding
    nexus data_template.

    Args:
        data_template (dict): The nomad-parser-nexus data template.

    Yields:
        an Instrument object according to the schema in hall_instrument.py
    '''

    instrument_sections = ['/systemparameters/',
                           '/measurement_state_machine/',
                           '/keithley_182/',
                           '/keithley_220/',
                           '/keithley_2000/',
                           '/keithley_2182/',
                           '/keithley_2400/',
                           '/keithley_2700/',
                           '/keithley_6485/',
                           '/keithley_7001/',
                           '/temperature_controller/',
                           '/field_controller/']

    hall_instrument = Instrument()
    hall_instrument.keithley_7001 = Keithley7001()
    hall_instrument.keithley_6485 = Keithley6485()
    hall_instrument.keithley_220 = Keithley220()
    hall_instrument.keithley_2000 = Keithley2000()
    hall_instrument.keithley_2182 = Keithley2182()
    hall_instrument.keithley_182 = Keithley182()
    hall_instrument.keithley_2700 = Keithley2700()
    hall_instrument.keithley_2400 = Keithley2400()
    hall_instrument.temperature_controller = TemperatureController()
    hall_instrument.field_controller = FieldController()

    temperature_domains: dict = {}

    for key in data_template:

        clean_key = clean_key_str(key)

        regexp = re.compile('temperature_domain_(\\d+)/')
        if bool(regexp.search(clean_key)):
            data_index = regexp.search(clean_key).group(1)
            if data_index not in temperature_domains:
                temperature_domains[data_index] = TemperatureDomain()
            clean_dkey = clean_key.split(f'temperature_domain_{data_index}/')[1]
            if hasattr(temperature_domains[data_index], clean_dkey):
                setattr(temperature_domains[data_index], clean_dkey, data_template[key])
            continue

        for instrument_section in instrument_sections:
            if instrument_section in clean_key:
                field_key = clean_key.split(str(instrument_section))[1]
                if hasattr(hall_instrument, field_key):
                    setattr(hall_instrument, field_key, data_template[key])
                elif (hasattr(hall_instrument, str(instrument_section.replace('/', '')))
                      and hasattr(getattr(hall_instrument, str(instrument_section.replace('/', ''))), field_key)):
                    setattr(getattr(hall_instrument, str(instrument_section.replace('/', ''))), field_key, data_template[key])

    for t_domain in temperature_domains.values():
        hall_instrument.temperature_domain.append(t_domain)

    return hall_instrument
