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
from numpy.core.fromnumeric import repeat

from nomad.metainfo import (
    MSection, MCategory, Package, Quantity, SubSection, Datetime)


m_package = Package(name='experimental_common')


class UserProvided(MCategory):
    pass


class SectionWithNotes(MSection):
    '''
    A common base-class for sections that should contain optional user provided
    notes.
    '''
    notes = Quantity(type=str, categories=[UserProvided])


class DeviceSettings(MSection):
    device_name = Quantity(type=str)
    analysis_method = Quantity(type=str)
    analyzer_lens = Quantity(type=str)
    analyzer_slit = Quantity(type=str)
    scan_mode = Quantity(type=str)
    detector_voltage = Quantity(type=str)
    workfunction = Quantity(type=str)
    channel_id = Quantity(type=str)
    max_energy = Quantity(type=str)
    min_energy = Quantity(type=str)
    guntype = Quantity(type=str)
    beam_energy = Quantity(type=str)
    resolution = Quantity(type=str)
    step_size = Quantity(type=str)
    acquisition_mode = Quantity(type=str)
    beam_current = Quantity(type=str)
    detector_type = Quantity(type=str)
    dark_current = Quantity(type=str)


class SampleMaterial(MSection):
    ''' This section describes a sample's material. '''
    elements = Quantity(
        type=str, shape=['*'],
        description='A list of element symbols for chemical elements in the material.')

    formula = Quantity(
        type=str, description='The chemical formula that describes the material.')

    name = Quantity(
        type=str, description='The name that describes the material.')

    space_group_number = Quantity(
        type=int, description='Space group of the material (if crystalline).')


class Sample(SectionWithNotes):
    sample_id = Quantity(
        type=str, description='Identification number or signatures of the sample used.')

    sample_name = Quantity(
        type=str, description='A human readable free text name for the sample.')

    sample_description = Quantity(
        type=str, description='A description of the sample.')

    sample_state = Quantity(
        type=str, description='The physical state of the sample.')

    sample_temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin',
        description='The temperature of the sample during the experiment.')

    sample_microstructure = Quantity(
        type=str, description='The sample microstructure.')

    sample_constituents = Quantity(
        type=str, description='The constituents.')

    section_material = SubSection(sub_section=SampleMaterial, repeats=True)


class ExperimentLocation(MSection):
    address = Quantity(
        type=str, description='''
        The address where the experiment took place, format 'Country, City, Street'
        ''')

    institution = Quantity(
        type=str, description='''
        Name of the institution hosting the experimental facility (e.g. in full or an
        acronym).
        ''')

    facility = Quantity(
        type=str, description='''
        Name of the experimental facility (e.g. in full or an acronym).
        ''')


class Experiment(SectionWithNotes):
    method_name = Quantity(type=str)
    method_abbreviation = Quantity(type=str)
    experiment_id = Quantity(type=str)
    experiment_publish_time = Quantity(
        type=Datetime, description='The datetime when this experiment was published.')
    experiment_start_time = Quantity(
        type=Datetime, description='The datetime of the beginning of the experiment.')
    experiment_end_time = Quantity(
        type=Datetime, description='The datetime of the experiment end.')
    edges = Quantity(type=str, shape=['*'])
    description = Quantity(type=str)

    experiment_location = SubSection(sub_section=ExperimentLocation)


class Instrument(SectionWithNotes):
    n_scans = Quantity(type=str)
    dwell_time = Quantity(type=str)
    excitation_energy = Quantity(type=str)
    source_label = Quantity(type=str)

    section_device_settings = SubSection(sub_section=DeviceSettings, repeats=True)


class Author(MSection):
    author_name = Quantity(type=str)
    author_profile_url = Quantity(type=str)
    author_profile_api_url = Quantity(type=str)
    group_name = Quantity(type=str)


class Origin(MSection):
    '''
    A section that describes a potential foreign origin of an entry. An example would
    be data found on zenodo.org or a specialized database like eelsdb.org.
    '''
    # TODO This is not experiment specific and should be added to the general NOMAD
    # entry metadata.

    authors = SubSection(sub_section=Author, repeats=True)
    permalink = Quantity(type=str)
    api_permalink = Quantity(type=str)
    repository_name = Quantity(
        type=str, description='The name of the repository, where the data is stored.')

    repository_url = Quantity(
        type=str, description='An URL to the repository, where the data is stored.')

    preview_url = Quantity(
        type=str, description='An URL to an image file that contains a preview.')

    entry_repository_url = Quantity(
        type=str, description='An URL to the entry on the repository, where the data is stored.')


class DataHeader(MSection):
    channel_id = Quantity(type=str)
    label = Quantity(type=str)
    unit = Quantity(type=str)


class Metadata(MSection):
    section_sample = SubSection(sub_section=Sample)
    section_experiment = SubSection(sub_section=Experiment)
    section_instrument = SubSection(sub_section=Instrument)
    section_origin = SubSection(sub_section=Origin)


class Spectrum(MSection):
    n_values = Quantity(type=int)
    count = Quantity(type=np.dtype(np.float64), shape=['n_values'], description='The count at each energy value, dimensionless')
    energy = Quantity(type=np.dtype(np.float64), shape=['n_values'], unit='J', description='The energy range of the spectrum')

    spectrum_region = Quantity(type=str, shape=[])

    more_channel_data = Quantity(
        type=np.dtype(np.float64), shape=['n_channels', 'n_values'],
        description='Additional channel data according to `more_channel_data_header`.')
    n_more_channels = Quantity(type=int, derived=lambda spectrum: len(spectrum.more_channel_data))
    more_channel_data_header = SubSection(sub_section=DataHeader, repeats=True)


class Data(MSection):
    section_spectrum = SubSection(sub_section=Spectrum)


class Measurement(SectionWithNotes):
    section_metadata = SubSection(sub_section=Metadata)
    section_data = SubSection(sub_section=Data)


m_package.__init_metainfo__()
