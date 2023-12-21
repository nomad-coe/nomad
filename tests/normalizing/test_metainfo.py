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

from nomad.datamodel import EntryData, EntryArchive
from nomad.metainfo import Quantity, SubSection
from nomad.client import normalize_all


def test_normalizer_level():
    class SubSection0(EntryData):
        number = Quantity(type=np.int32)

        def normalize(self, archive, logger):
            archive.data.numbers = np.append(archive.data.numbers, self.number)

    class SubSection1(SubSection0):
        sub_section = SubSection(sub_section=SubSection0)
        normalizer_level = 2

    class SubSection2(SubSection0):
        sub_section = SubSection(sub_section=SubSection0)
        normalizer_level = 1

    class Root(EntryData):
        numbers = Quantity(type=np.int32, shape=['*'])
        sub_section_0 = SubSection(sub_section=SubSection0)
        sub_section_1 = SubSection(sub_section=SubSection1)
        sub_section_2 = SubSection(sub_section=SubSection2, repeats=True)

    archive = EntryArchive(data=Root(numbers=[]))
    archive.data.sub_section_0 = SubSection0(number=5)
    archive.data.sub_section_1 = SubSection1(
        sub_section=SubSection0(number=1), number=2
    )
    archive.data.sub_section_2.append(
        SubSection2(sub_section=SubSection0(number=3), number=4)
    )
    archive.data.sub_section_2.append(
        SubSection2(sub_section=SubSection0(number=6), number=7)
    )

    normalize_all(archive)
    assert (archive.data.numbers == [5, 3, 4, 6, 7, 1, 2]).all()
