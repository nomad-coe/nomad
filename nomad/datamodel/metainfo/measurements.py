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

from nomad.metainfo import (
    MSection, Package, Quantity, SubSection, Datetime)
from nomad.metainfo.metainfo import Reference, SectionProxy


m_package = Package(name='measurements')


class Sample(MSection):
    sample_id = Quantity(
        type=str, description='Identification number or signatures of the sample used.')

    name = Quantity(
        type=str, description='A human readable free text name for the sample.')

    description = Quantity(
        type=str, description='A description of the sample.')

    sample_state = Quantity(
        type=str, description='The physical state of the sample.')

    sample_temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin',
        description='The temperature of the sample during the measurement.')

    sample_microstructure = Quantity(
        type=str, description='The sample microstructure.')

    sample_constituents = Quantity(
        type=str, description='The constituents.')

    elements = Quantity(type=str, shape=["*"])
    chemical_formula = Quantity(type=str)


class Experiment(MSection):
    name = Quantity(
        type=str, description='A human readable free text name for the experiment.')

    description = Quantity(
        type=str, description='A description of the experiment.')

    steps = Quantity(
        type=str, shape=['*'], description='Human readable experiment steps.')

    sample = SubSection(section_def=Sample, description='The used sample.')
    sample_ref = Quantity(type=Reference(Sample.m_def), description='Reference to the used sample.')

    measurement = SubSection(
        section_def=SectionProxy('Measurement'),
        repeats=True,
        description='Measurements performed in this experiment.')


class Instrument(MSection):
    instrument_id = Quantity(
        type=str, description='Identification number or signatures of the instrument used.')

    name = Quantity(
        type=str, description='A human readable free text name for the instrument.')

    description = Quantity(
        type=str, description='A description of the instrument.')


class Measurement(MSection):
    measurement_id = Quantity(type=str)
    name = Quantity(type=str)
    description = Quantity(type=str)

    method_name = Quantity(type=str)
    method_abbreviation = Quantity(type=str)

    start_time = Quantity(
        type=Datetime, description='The datetime of the beginning of the measurement.')
    end_time = Quantity(
        type=Datetime, description='The datetime of the measurement end.')

    facility = Quantity(
        type=str, description='''
        Description of the facility (e.g. in full or an acronym) where
        the measurement was conducted.''')

    sample = SubSection(section_def=Sample, repeats=True)
    instrument = SubSection(section_def=Instrument, repeats=True)


class Spectrum(MSection):
    '''
    Generic spectrum data with energies and counts. May include additional
    channels.
    '''

    class SpectrumChannel(MSection):
        '''
        Provides the metadata for a generic additional spectrum channel. Do not
        use it for energy or count; they have their predefined channels.
        '''
        channel_id = Quantity(type=str)
        label = Quantity(type=str)
        unit = Quantity(type=str)

    def derive_n_values(self):
        if self.count is not None:
            return len(self.count)
        if self.energy is not None:
            return len(self.energy)
        else:
            return 0

    n_values = Quantity(type=int, derived=derive_n_values)
    count = Quantity(
        type=np.dtype(np.float64), shape=['n_values'],
        description='The count at each energy value, dimensionless')
    energy = Quantity(
        type=np.dtype(np.float64), shape=['n_values'], unit='J',
        description='The energy range of the spectrum')

    additional_channel_data = Quantity(
        type=np.dtype(np.float64), shape=['n_channels', 'n_values'],
        description='''
            Data from additional channels. The channels are described in `additional channels`.
        ''')
    additional_channels = SubSection(
        sub_section=SpectrumChannel, repeats=True,
        description='''
            Metadata for additional channels. The order is the same as the channel data
            appears in `additional_channel_data`.
        ''')
    n_additional_channels = Quantity(type=int, derived=lambda spectrum: len(spectrum.additional_channels))


m_package.__init_metainfo__()
