# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest
import numpy as np

from nomadcore.local_meta_info import InfoKindEl, InfoKindEnv

from nomad.metainfo.metainfo import (
    MSection, MCategory, Section, Quantity, SubSection, Definition, Package, DeriveError,
    MetainfoError, Environment, MResource, Datetime, units, Annotation, SectionAnnotation,
    DefinitionAnnotation, Reference, MProxy, derived)
from nomad.metainfo.legacy import LegacyMetainfoEnvironment, convert, python_package_mapping
from nomad.parsing.metainfo import MetainfoBackend


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
        'vasp_incars',
        'dependencies/parsers/vasp/vaspparser/metainfo/vasp_incars.py',
        'vaspparser.metainfo.vasp_incars')
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


def test_backend(env, no_warn):
    backend = MetainfoBackend(env)
    run = backend.openSection('section_run')
    backend.addValue('program_name', 'vasp')

    system_0 = backend.openSection('section_system')
    assert system_0 == 0

    backend.addValue('number_of_atoms', 3)
    backend.addArrayValues('atom_labels', np.array(['H', 'H', 'O']))
    backend.addArrayValues('atom_positions', [[0, 0, 0], [0, 0, 0], [0, 0, 0]])
    backend.closeSection('section_system', system_0)

    system_1 = backend.openSection('section_system')
    assert system_1 == 1
    backend.closeSection('section_system', system_1)

    method = backend.openSection('section_method')
    backend.addValue('method_system_ref', system_0)
    backend.closeSection('section_method', method)

    backend.closeSection('section_run', run)

    assert backend.get_sections('section_system') == [0, 1]
    assert len(backend.get_value('atom_labels', 0)) == 3
    assert backend.get_value('method_system_ref', 0) == 0
    assert backend.get_value('program_name', 0) == 'vasp'

    backend.openContext('section_run/0/section_system/1')
    backend.addValue('number_of_atoms', 10)
    backend.closeContext('section_run/0/section_system/1')
    assert backend.get_value('number_of_atoms', 1) == 10
