from nomad.metainfo import Quantity, Package, Section, MEnum, Datetime, MSection, SubSection, SectionProxy
from nomad.datamodel.data import EntryData
from nomad.datamodel.metainfo.annotations import ELNAnnotation, ELNComponentEnum
import numpy as np

m_package = Package()


class MySection(MSection):
    m_def = Section()
    name = Quantity(
        type=str,
        a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity),
        description='For testing subsection quantity.'
    )


class MySectionRecursiveA(MSection):
    child = SubSection(section_def=SectionProxy("MySectionRecursiveB"))
    name_a = Quantity(type=str)


class MySectionRecursiveB(MSection):
    child = SubSection(section_def=MySectionRecursiveA)
    name_b = Quantity(type=str)


class MySchema(EntryData):
    name = Quantity(
        type=str,
        a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity),
        description='For testing string field.'
    )
    message = Quantity(
        type=MEnum(['A', 'B']),
        a_eln=ELNAnnotation(component=ELNComponentEnum.EnumEditQuantity),
        description='For testing enum field.'
    )
    empty = Quantity(
        type=str,
        a_eln=ELNAnnotation(component=ELNComponentEnum.StringEditQuantity),
        description='For testing empty field.'
    )
    valid = Quantity(
        type=bool,
        a_eln=ELNAnnotation(component=ELNComponentEnum.BoolEditQuantity),
        description='For testing boolean field.'
    )
    count = Quantity(
        type=int,
        a_eln=ELNAnnotation(component=ELNComponentEnum.NumberEditQuantity),
        description='For testing integer field.'
    )
    frequency = Quantity(
        type=float,
        unit='1/s',
        a_eln=ELNAnnotation(component=ELNComponentEnum.NumberEditQuantity),
        description='For testing floating point field.'
    )
    timestamp = Quantity(
        type=Datetime,
        a_eln=ELNAnnotation(component=ELNComponentEnum.DateTimeEditQuantity),
        description='For testing datetime field.'
    )
    reference_section = Quantity(
        type=MySection,
        a_eln=ELNAnnotation(component=ELNComponentEnum.ReferenceEditQuantity),
        description='For testing section reference.'
    )
    non_scalar = Quantity(
        type=np.float64,
        shape=[3, 3],
        description='For testing non-scalar field.'
    )

    child = SubSection(section_def=MySection, repeats=False)
    child_repeating = SubSection(section_def=MySection, repeats=True)
    child_recursive = SubSection(section_def=MySectionRecursiveA)

    def normalize(self, archive, logger):
        super(MySchema, self).normalize(archive, logger)


m_package.__init_metainfo__()
