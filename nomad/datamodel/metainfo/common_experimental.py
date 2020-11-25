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
    MSection, Package, Quantity, SubSection, SectionProxy, Section,
    Datetime, JSON)


m_package = Package(name='experimental_common')


class Experiment(MSection):
    '''
    The root section for all (meta)data that belongs to a single experiment.
    '''
    m_def = Section(validate=False)

    experiment_summary = Quantity(
        type=str, description='A descriptive summary of the content of the experiment.')

    experiment_location = SubSection(sub_section=SectionProxy('Location'))

    experiment_publish_time = Quantity(
        type=Datetime, description='The datetime when this experiment was published.')

    experiment_time = Quantity(
        type=Datetime, description='The datetime of the beginning of the experiment.')

    experiment_end_time = Quantity(
        type=Datetime, description='The datetime of the experiment end.')

    raw_metadata = Quantity(
        type=JSON, description='The whole or partial metadata in its original source JSON format.')

    section_data = SubSection(sub_section=SectionProxy('Data'))

    section_method = SubSection(sub_section=SectionProxy('Method'))

    section_sample = SubSection(sub_section=SectionProxy('Sample'))


class Location(MSection):
    m_def = Section(validate=False)

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


class Data(MSection):
    '''
    This section contains information about the stored data.
    '''
    m_def = Section(validate=False)

    repository_name = Quantity(
        type=str, description='The name of the repository, where the data is stored.')

    repository_url = Quantity(
        type=str, description='An URL to the repository, where the data is stored.')

    preview_url = Quantity(
        type=str, description='An URL to an image file that contains a preview.')

    entry_repository_url = Quantity(
        type=str, description='An URL to the entry on the repository, where the data is stored.')


class Method(MSection):
    '''
    This section contains information about the applied experimental method.
    '''
    m_def = Section(validate=False)

    data_type = Quantity(
        type=str, description='Name of the type of data that the experiment produced.')

    method_name = Quantity(
        type=str, description='Full name of the experimental method in use')

    method_abbreviation = Quantity(
        type=str,
        description='Abbreviated name (i.e. acronym) of the experimental method')

    probing_method = Quantity(
        type=str, description='The probing method used')

    instrument_description = Quantity(
        type=str, description='A description of the instrumentation used for the experiment.')


class Sample(MSection):
    '''
    The section for all sample related (meta)data that was used in the experiment.
    '''
    m_def = Section(validate=False)

    sample_description = Quantity(
        type=str, description='Description of the sample used in the experiment.')

    sample_id = Quantity(
        type=str, description='Identification number or signatures of the sample used.')

    sample_state = Quantity(
        type=str, description='The physical state of the sample.')

    sample_temperature = Quantity(
        type=np.dtype(np.float64), unit='kelvin',
        description='The temperature of the sample during the experiment.')

    sample_microstructure = Quantity(
        type=str, description='The sample microstructure.')

    sample_constituents = Quantity(
        type=str, description='The constituents.')

    section_material = SubSection(sub_section=SectionProxy('Material'))


class Material(MSection):
    ''' This section describes a sample's material. '''
    m_def = Section(validate=False)

    chemical_formula = Quantity(
        type=str, description='The chemical formula that describes the sample.')

    chemical_name = Quantity(
        type=str, description='The chemical name that describes the sample.')

    atom_labels = Quantity(
        type=str, shape=['number_of_elements'],
        description='Atom labels for distinct elements in the sample.')

    number_of_elements = Quantity(
        type=int, derived=lambda m: len(m.atom_labels) if m.atom_labels else 0,
        description='Number of distinct chemical elements in the sample.')

    space_group = Quantity(
        type=int, description='Space group of the sample compound (if crystalline).')


m_package.__init_metainfo__()
