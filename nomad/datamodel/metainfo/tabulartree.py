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

import typing

from nomad.metainfo import (
    MSection, Package, Quantity, SubSection)


m_package = Package(name='tabulartree')


class TabularTreeNodeInfo(MSection):
    value = Quantity(type=typing.Any)
    description = Quantity(type=str)
    unit = Quantity(type=str)


class TabularTreeLevel3(MSection):
    name = Quantity(type=str, default='<node name?>')
    info = SubSection(sub_section=TabularTreeNodeInfo)


class TabularTreeLevel2(MSection):
    name = Quantity(type=str, default='<node name?>')
    info = SubSection(sub_section=TabularTreeNodeInfo)
    nodes = SubSection(sub_section=TabularTreeLevel3, repeats=True)


class TabularTreeLevel1(MSection):
    name = Quantity(type=str, default='<node name?>')
    info = SubSection(sub_section=TabularTreeNodeInfo)
    nodes = SubSection(sub_section=TabularTreeLevel2, repeats=True)


class TabularTree(MSection):
    name = Quantity(type=str, default='<node name?>')
    info = SubSection(sub_section=TabularTreeNodeInfo)
    nodes = SubSection(sub_section=TabularTreeLevel1, repeats=True)


m_package.__init_metainfo__()
