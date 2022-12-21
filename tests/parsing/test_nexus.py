"""This is a code that performs several tests on nexus tool

"""
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

from typing import Any

import pytest

from nomad.datamodel import EntryArchive
from nomad.metainfo import Section
from nomad.metainfo.nexus import nexus_metainfo_package
from nomad.parsing.nexus import NexusParser
from nomad.units import ureg
from nomad.utils import get_logger


@pytest.mark.parametrize('path,value', [
    pytest.param('name', 'nexus'),
    pytest.param('NXobject.name', 'NXobject'),
    pytest.param('NXentry.nx_kind', 'group'),
    pytest.param('NXentry.NXdata', '*'),
    pytest.param('NXdetector.real_time__field', '*'),
    pytest.param('NXentry.NXdata.nx_optional', True),
    pytest.param('NXentry.DATA.section_def.nx_kind', 'group'),
    pytest.param('NXentry.DATA.section_def.nx_optional', True),
    pytest.param('NXentry.DATA.section_def.name', 'NXdata'),
    pytest.param('NXdetector.real_time__field.name', 'real_time__field'),
    pytest.param('NXdetector.real_time__field.nx_type', 'NX_NUMBER'),
    pytest.param('NXdetector.real_time__field.nx_units', 'NX_TIME'),
    pytest.param('NXarpes.NXentry.NXdata.nx_optional', False),
    pytest.param('NXentry.nx_category', 'base'),
    pytest.param('NXapm.nx_category', 'application')
])
def test_assert_nexus_metainfo(path: str, value: Any):
    '''
    Test the existence of nexus metainfo
    '''
    current = nexus_metainfo_package
    for name in path.split('.'):
        for content in current.m_contents():
            if getattr(content, 'name', None) == name:
                current = content  # type: ignore
                break

        else:
            current = getattr(current, name, None)

        if current is None:
            assert False, f'{path} does not exist'

    if value == '*':
        assert current is not None, f'{path} does not exist'
    elif value is None:
        assert current is None, f'{path} does exist'
    else:
        assert current == value, f'{path} has wrong value'

    if isinstance(current, Section):
        assert current.nx_kind is not None
        for base_section in current.all_base_sections:
            assert base_section.nx_kind == current.nx_kind


def test_nexus_example():
    archive = EntryArchive()

    example_data = 'tests/data/parsers/nexus/201805_WSe2_arpes.nxs'
    NexusParser().parse(example_data, archive, get_logger(__name__))
    assert archive.nexus.NXarpes.ENTRY[0].SAMPLE[0].pressure__field == ureg.Quantity('3.27e-10*millibar')

    instrument = archive.nexus.NXarpes.ENTRY[0].INSTRUMENT[0]

    assert instrument.monochromator.energy__field == ureg.Quantity('36.49699020385742*electron_volt')
    assert instrument.analyser.entrance_slit_size__field == ureg.Quantity('750 micrometer')
    # good ENUM - x-ray
    assert instrument.SOURCE[0].probe__field == 'x-ray'
    # wrong inherited ENUM - Burst
    assert instrument.SOURCE[0].mode__field is None
    # wrong inherited ENUM for extended field - 'Free Electron Laser'
    assert instrument.SOURCE[0].type__field is None

    data = archive.nexus.NXarpes.ENTRY[0].DATA[0]
    assert data.angles__field is not None
    # assert data.delays is not None
    assert data.energies__field is not None
    assert data.angles__field.check("1/Å")
    # assert data.delays.check("fs")
    assert data.energies__field.check("eV")


def test_same_name_field_and_group():
    archive = EntryArchive()
    example_data = 'tests/data/parsers/nexus/SiO2onSi.ellips.nxs'
    NexusParser().parse(example_data, archive, get_logger(__name__))
    archive.m_to_dict(with_out_meta=True)
