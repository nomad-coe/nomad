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

from nomad.datamodel.data import EntryData
from nomad.datamodel.metainfo.annotations import ELNAnnotation, ELNComponentEnum
from nomad.metainfo import (
    Datetime,
    MEnum,
    MSection,
    Package,
    Quantity,
    Section,
    SectionProxy,
    SubSection,
)

m_package = Package()


class MySection(MSection):
    m_def = Section()
    name = Quantity(
        type=str,
        a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity),
        description='For testing subsection string quantity.',
    )
    count = Quantity(
        type=int,
        a_eln=ELNAnnotation(component=ELNComponentEnum.NumberEditQuantity),
        description='For testing subsection integer quantity.',
    )
    frequency = Quantity(
        type=float,
        unit='1/s',
        a_eln=ELNAnnotation(component=ELNComponentEnum.NumberEditQuantity),
        description='For testing subsection floating point quantity.',
    )


class MySectionRecursiveA(MSection):
    child = SubSection(section_def=SectionProxy('MySectionRecursiveB'))
    name_a = Quantity(type=str)


class MySectionRecursiveB(MSection):
    child = SubSection(section_def=MySectionRecursiveA)
    name_b = Quantity(type=str)


class MyBaseSchemaA(EntryData):
    inherited_a = Quantity(
        type=str,
        a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity),
        description='For testing inherited quantities.',
    )


class MyBaseSchemaB(MSection):
    inherited_b = Quantity(
        type=str,
        a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity),
        description='For testing multiple inheritance.',
    )


class MySchema(MyBaseSchemaA, MyBaseSchemaB):
    name = Quantity(
        type=str,
        a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity),
        description='For testing string field.',
    )
    message = Quantity(
        type=MEnum(['A', 'B']),
        a_eln=ELNAnnotation(component=ELNComponentEnum.EnumEditQuantity),
        description='For testing enum field.',
    )
    empty = Quantity(
        type=str,
        a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity),
        description='For testing empty field.',
    )
    valid = Quantity(
        type=bool,
        a_eln=ELNAnnotation(component=ELNComponentEnum.BoolEditQuantity),
        description='For testing boolean field.',
    )
    count = Quantity(
        type=int,
        a_eln=ELNAnnotation(component=ELNComponentEnum.NumberEditQuantity),
        description='For testing integer field.',
    )
    frequency = Quantity(
        type=float,
        unit='1/s',
        a_eln=ELNAnnotation(component=ELNComponentEnum.NumberEditQuantity),
        description='For testing floating point field.',
    )
    timestamp = Quantity(
        type=Datetime,
        a_eln=ELNAnnotation(component=ELNComponentEnum.DateTimeEditQuantity),
        description='For testing datetime field.',
    )
    reference_section = Quantity(
        type=MySection,
        a_eln=ELNAnnotation(component=ELNComponentEnum.ReferenceEditQuantity),
        description='For testing section reference.',
    )
    non_scalar = Quantity(
        type=np.float64, shape=[3, 3], description='For testing non-scalar field.'
    )

    child = SubSection(section_def=MySection, repeats=False)
    child_repeating = SubSection(section_def=MySection, repeats=True)
    child_recursive = SubSection(section_def=MySectionRecursiveA)

    def normalize(self, archive, logger):
        super().normalize(archive, logger)


m_package.__init_metainfo__()
