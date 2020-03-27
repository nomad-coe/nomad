import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import general_experimental

m_package = Package(
    name='general_experimental_method_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='general.experimental.method.nomadmetainfo.json'))


class section_method(MSection):
    '''
    This section contains information about the applied experimental method.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_method'))

    experiment_method_name = Quantity(
        type=str,
        shape=[],
        description='''
        Full name of the experimental method in use
        ''',
        a_legacy=LegacyDefinition(name='experiment_method_name'))

    experiment_method_abbreviation = Quantity(
        type=str,
        shape=[],
        description='''
        Abbreviated name (i.e. acronym) of the experimental method
        ''',
        a_legacy=LegacyDefinition(name='experiment_method_abbreviation'))

    equipment_description = Quantity(
        type=str,
        shape=[],
        description='''
        Name or model of the equipment (e.g. in full or an acronym).
        ''',
        a_legacy=LegacyDefinition(name='equipment_description'))

    probing_method = Quantity(
        type=str,
        shape=[],
        description='''
        The probing method used
        ''',
        a_legacy=LegacyDefinition(name='probing_method'))


class section_experiment(general_experimental.section_experiment):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_experiment'))

    section_method = SubSection(
        sub_section=SectionProxy('section_method'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_method'))


m_package.__init_metainfo__()
