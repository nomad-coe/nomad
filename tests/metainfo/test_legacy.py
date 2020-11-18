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

import pytest
import numpy as np

from nomadcore.local_meta_info import InfoKindEl, InfoKindEnv

from nomad.metainfo.metainfo import Section, Quantity, SubSection
from nomad.metainfo.legacy import LegacyMetainfoEnvironment, convert, python_package_mapping


@pytest.fixture(scope='session')
def legacy_example():
    return {
        'type': 'nomad_meta_info_1_0',
        'metaInfos': [
            dict(name='section_run', kindStr='type_section'),
            dict(
                name='program_name', dtypeStr='C', superNames=['program_info']),
            dict(name='bool_test', dtypeStr='b', superNames=['section_run']),
            dict(name='section_system', kindStr='type_section', superNames=['section_run']),
            dict(
                name='atom_labels', dtypeStr='C', shape=['number_of_atoms'],
                superNames=['section_system']),
            dict(
                name='atom_positions', dtypeStr='f', shape=['number_of_atoms', 3],
                superNames=['section_system']),
            dict(
                name='number_of_atoms', dtypeStr='i', kindStr='type_dimension',
                superNames=['section_system', 'system_info']),
            dict(name='section_method', kindStr='type_section', superNames=['section_run']),
            dict(
                name='method_system_ref', dtypeStr='i',
                referencedSections=['section_system'],
                superNames=['section_method']),
            dict(
                name='program_info', kindStr='type_abstract_document_content',
                superNames=['section_run']),
            dict(name='system_info', kindStr='type_abstract_document_content')
        ]
    }


@pytest.fixture(scope='session')
def legacy_env(legacy_example):
    env = InfoKindEnv()
    env.name = 'test.nomadmetainfo.json'
    for definition in legacy_example.get('metaInfos'):
        env.addInfoKindEl(InfoKindEl(
            description='test_description', package='test_package', **definition))
    return env


@pytest.fixture(scope='session')
def env(legacy_env):
    return convert(legacy_env)


@pytest.mark.parametrize('package,path,name', [
    (
        'vasp',
        'dependencies/parsers/vasp/vaspparser/metainfo/vasp.py',
        'vaspparser.metainfo.vasp'),
    (
        'common',
        'nomad/datamodel/metainfo/common.py',
        'nomad.datamodel.metainfo.common'),
    (
        'fhi_aims',
        'dependencies/parsers/fhi-aims/fhiaimsparser/metainfo/fhi_aims.py',
        'fhiaimsparser.metainfo.fhi_aims')
])
def test_package_mapping(package, path, name):
    assert python_package_mapping(package) == (name, path)


def test_environment(env: LegacyMetainfoEnvironment, no_warn):
    assert env.packages[0].name == 'test_package'
    assert 'section_system' in env.packages[0].all_definitions
    section_run_def = env.resolve_definition('section_run', Section)
    assert section_run_def is not None
    assert section_run_def.section_cls is not None
    assert env.resolve_definition('section_run', Section) is not None
    assert env.resolve_definition('section_system', Section) is not None
    assert env.resolve_definition('section_system', SubSection) is not None

    assert env.resolve_definition('method_system_ref', Quantity).type.target_section_def == \
        env.resolve_definition('section_system', Section)
    assert env.resolve_definition('atom_labels', Quantity).type == str
    assert env.resolve_definition('atom_positions', Quantity).type == np.dtype('float64')
    assert env.resolve_definition('bool_test', Quantity).type == bool
    assert env.resolve_definition('program_name', Quantity).m_parent.name == 'section_run'
