import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import general_experimental

m_package = Package(
    name='general_experimental_sample_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='general.experimental.sample.nomadmetainfo.json'))


class section_sample(MSection):
    '''
    The section for all sample related (meta)data that was used in the experiment.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_sample'))

    sample_description = Quantity(
        type=str,
        shape=[],
        unit='dimensionless',
        description='''
        Description of the sample used in the experiment.
        ''',
        a_legacy=LegacyDefinition(name='sample_description'))

    sample_id = Quantity(
        type=str,
        shape=[],
        unit='dimensionless',
        description='''
        Identification number or signatures of the sample used.
        ''',
        a_legacy=LegacyDefinition(name='sample_id'))

    sample_state = Quantity(
        type=str,
        shape=[],
        description='''
        The physical state of the sample.
        ''',
        a_legacy=LegacyDefinition(name='sample_state'))

    sample_chemical_formula = Quantity(
        type=str,
        shape=[],
        description='''
        The chemical formula that describes the sample
        ''',
        a_legacy=LegacyDefinition(name='sample_chemical_formula'))

    sample_chemical_name = Quantity(
        type=str,
        shape=[],
        description='''
        The chemical name that describes the sample
        ''',
        a_legacy=LegacyDefinition(name='sample_chemical_name'))

    sample_atom_labels = Quantity(
        type=str,
        shape=['n'],
        description='''
        The chemical name that describes the sample
        ''',
        a_legacy=LegacyDefinition(name='sample_atom_labels'))

    number_of_elements = Quantity(
        type=int,
        shape=[],
        description='''
        Number of distinct chemical elements in the sample.
        ''',
        a_legacy=LegacyDefinition(name='number_of_elements'))

    sample_space_group = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        unit='dimensionless',
        description='''
        Space group of the sample compound (if crystalline).
        ''',
        a_legacy=LegacyDefinition(name='sample_space_group'))

    sample_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='kelvin',
        description='''
        The temperature of the sample during the experiment in K.
        ''',
        a_legacy=LegacyDefinition(name='sample_temperature'))

    sample_microstructure = Quantity(
        type=str,
        shape=[],
        description='''
        The sample microstructure
        ''',
        a_legacy=LegacyDefinition(name='sample_microstructure'))

    sample_constituents = Quantity(
        type=str,
        shape=[],
        description='''
        The constituents
        ''',
        a_legacy=LegacyDefinition(name='sample_constituents'))


class section_experiment(general_experimental.section_experiment):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_experiment'))

    section_sample = SubSection(
        sub_section=SectionProxy('section_sample'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_sample'))


m_package.__init_metainfo__()
