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
import pytest
import yaml

from nomad.utils import strip
from nomad.metainfo import (
    Package, MSection, Quantity, Reference, SubSection, Section, MProxy, MetainfoError,
    Context)

m_package = Package()


class Sample(MSection):
    sample_id = Quantity(
        type=str,
        a_eln=dict(component="StringEditQuantity"),
        description='''
        This is a description with *markup* using [markdown](https://markdown.org).
        It can have multiple lines, because yaml allows to easily do this.
        ''')


class Process(MSection):
    samples = Quantity(type=Reference(Sample.m_def), shape=["*"])
    layers = SubSection(sub_section=Sample.m_def, repeats=True)


class SpecialProcess(Process):
    values = Quantity(type=np.dtype(np.float64), shape=[3, 3])


m_package.__init_metainfo__()


def yaml_to_package(yaml_str):
    class MyContext(Context):
        def create_reference(
                self, section: MSection, quantity_def: Quantity, value: MSection,
                global_reference: bool = False) -> str:
            if section.m_root() == value.m_root():  # type: ignore
                return super().create_reference(section, quantity_def, value)
            return None

    yaml_obj = yaml.safe_load(yaml_str)
    package = Package.m_from_dict(yaml_obj, m_context=MyContext())
    return package


yaml_schema_example = strip('''
    m_def: 'nomad.metainfo.metainfo.Package'
    sections:
        Sample:
            base_section: 'nomad.datamodel.metainfo.measurements.Sample'
            quantities:
                sample_id:
                    type: str
                    description: |
                        This is a description with *markup* using [markdown](https://markdown.org).
                        It can have multiple lines, because yaml allows to easily do this.
                    m_annotations:
                        eln:
                            component: StringEditQuantity
        Process:
            quantities:
                samples:
                    type: '#/Sample'
                    shape: ['*']
            sub_sections:
                layers:
                    section_def: '#/Sample'
                    repeats: true
        SpecialProcess:
            base_section: '#/Process'
            quantities:
                values:
                    type: np.float64
                    shape: [3, 3]
''')


def test_yaml_deserialization():
    des_m_package = yaml_to_package(yaml_schema_example)

    sample = m_package['section_definitions'][0]
    process = m_package['section_definitions'][1]
    special_process = m_package['section_definitions'][2]

    des_sample = des_m_package['section_definitions'][0]
    des_process = des_m_package['section_definitions'][1]
    des_special_process = des_m_package['section_definitions'][2]

    def assert_referenced_section(section):
        assert section is not None
        section = section.m_resolved()
        assert isinstance(section, Section)
        assert not isinstance(section, MProxy)

    assert sample.name == des_sample.name == "Sample"

    sample_id = sample['quantities'][0]
    des_sample_id = des_sample['quantities'][0]
    assert sample_id.name == des_sample_id.name == 'sample_id'
    assert sample_id.type == des_sample_id.type == str
    assert sample_id.shape == des_sample_id.shape
    assert sample_id.m_annotations["eln"].component == des_sample_id.m_annotations["eln"].component
    assert sample_id.description == des_sample_id.description.rstrip('\n')

    assert process.name == des_process.name == "Process"

    samples = process['quantities'][0]
    des_samples = des_process['quantities'][0]
    assert samples.name == des_samples.name == 'samples'
    assert samples.shape == des_samples.shape == ["*"]
    assert_referenced_section(des_samples.type.target_section_def)
    assert samples.type.target_section_def.name == des_samples.type.target_section_def.name == 'Sample'

    layers = process['sub_sections'][0]
    des_layers = des_process['sub_sections'][0]
    assert layers.name == des_layers.name == 'layers'
    assert layers.repeats is des_layers.repeats is True
    assert_referenced_section(des_layers.section_def)
    assert layers.section_def.name == des_layers.section_def.name == 'Sample'

    assert special_process.name == des_special_process.name == "SpecialProcess"

    values = special_process['quantities'][0]
    des_values = des_special_process['quantities'][0]
    assert values.name == des_values.name == 'values'
    assert values.shape == des_values.shape == [3, 3]

    base_section = special_process['base_sections'][0]
    des_base_section = des_special_process['base_sections'][0]
    assert_referenced_section(des_base_section)
    assert base_section.name == des_base_section.name == 'Process'

    des_m_package.m_to_dict()


yaml_schema_example_extended_types = strip('''
    m_def: 'nomad.metainfo.metainfo.Package'
    sections:
        Sample:
            base_section: 'nomad.datamodel.metainfo.measurements.Sample'
            quantities:
                method:
                    type: string
                    m_annotations:
                        eln:
                            component: StringEditQuantity
                spin:
                    type: boolean
                    m_annotations:
                        eln:
                            component: BoolEditQuantity
''')


def test_yaml_extended_types_deserialization():
    des_m_package = yaml_to_package(yaml_schema_example_extended_types)

    des_sample = des_m_package['section_definitions'][0]

    assert des_sample.name == "Sample"

    method = des_sample['quantities'][0]
    assert method.name == 'method'
    assert method.type == str

    spin = des_sample['quantities'][1]
    assert spin.name == 'spin'
    assert spin.type == bool

    des_m_package.m_to_dict()


@pytest.mark.parametrize('yaml_schema, expected_error', [
    pytest.param(strip('''
        m_def: 'nomad.metainfo.metainfo.Package'
        sections:
            Process:
                quantities:
                    samples:
                        type: np.float6
    '''), 'float6 is not a valid numpy type.', id='wrong np type'),
    pytest.param(strip('''
        m_def: 'nomad.metainfo.metainfo.Package'
        sections:
            Process:
                quantities:
                    samples:
                        type: numpy.int3
    '''), 'int3 is not a valid numpy type.', id='wrong numpy type'),
    pytest.param(strip('''
        m_def: 'nomad.metainfo.metainfo.Package'
        sections:
            Process:
                quantities:
                    samples:
                        type: float
                        m_annotations: eln
    '''), 'The provided m_annotations is of a wrong type. str was provided.', id='wrong m_annotations')
])
def test_errors(yaml_schema, expected_error):
    with pytest.raises(Exception) as exception:
        yaml_to_package(yaml_schema)

    assert isinstance(exception.value, MetainfoError)
    assert exception.value.args[0] == expected_error


def test_sub_section_tree():
    yaml = yaml_to_package('''
      sections:
          Parent:
              sub_sections:
                  the_child:
                      sub_section:
                          quantities:
                              quantity:
                                  type: str
    ''')
    reference = yaml_to_package('''
        sections:
            Parent:
                sub_sections:
                    the_child:
                        section: TheChild
                sections:
                    TheChild:
                        quantities:
                            quantity:
                                type: str
    ''')

    assert yaml.m_to_dict() == reference.m_to_dict()


@pytest.mark.parametrize('source_type, target_type', [
    pytest.param(
        'nomad.datamodel.data.ArchiveSection',
        'nomad.datamodel.data.ArchiveSection',
        id='python'),
    pytest.param('MySection', '/section_definitions/0', id='yaml'),
])
def test_references(source_type, target_type):
    yaml = yaml_to_package(f'''
        sections:
            MySection:
                quantities:
                    reference:
                       type: {source_type}
    ''')

    assert yaml.m_to_dict()['section_definitions'][0]['quantities'][0]['type'] == {
        'type_kind': 'reference',
        'type_data': target_type
    }
