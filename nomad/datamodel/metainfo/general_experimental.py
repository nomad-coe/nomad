import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference, Datetime, JSON
)
from nomad.metainfo.legacy import LegacyDefinition


m_package = Package(
    name='general_experimental_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='general.experimental.nomadmetainfo.json'))


class section_experiment(MSection):
    '''
    The root section for all (meta)data that belongs to a single experiment.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_experiment'))

    experiment_summary = Quantity(
        type=str,
        shape=[],
        description='''
        A descriptive summary of the content of the experiment.
        ''',
        a_legacy=LegacyDefinition(name='experiment_summary'))

    experiment_location = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the city and country the experiment took place, format 'Country, City'
        ''',
        a_legacy=LegacyDefinition(name='experiment_location'))

    experiment_facility_institution = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the institution hosting the experimental facility (e.g. in full or an
        acronym).
        ''',
        a_legacy=LegacyDefinition(name='experiment_facility_institution'))

    experiment_facility_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the experimental facility (e.g. in full or an acronym).
        ''',
        a_legacy=LegacyDefinition(name='experiment_facility_name'))

    experiment_publish_time = Quantity(
        type=Datetime,
        description='''
        The datetime when this experiment was published.
        ''')

    experiment_time = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        unit='second',
        description='''
        The datetime of the beginning of the experiment in seconds since epoch.
        ''',
        a_legacy=LegacyDefinition(name='experiment_time'))

    experiment_end_time = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        unit='second',
        description='''
        The datetime of the experiment end in seconds since epoch.
        ''',
        a_legacy=LegacyDefinition(name='experiment_end_time'))

    raw_metadata = Quantity(
        type=JSON,
        description='''
        The whole or partial metadata in its original source JSON format.
        ''')

    section_data = SubSection(
        sub_section=SectionProxy('section_data'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_data'))

    section_method = SubSection(
        sub_section=SectionProxy('section_method'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_method'))

    section_sample = SubSection(
        sub_section=SectionProxy('section_sample'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_sample'))


class section_data(MSection):
    '''
    This section contains information about the stored data.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_data'))

    data_repository_name = Quantity(
        type=str,
        shape=[],
        description='''
        The name of the repository, where the data is stored.
        ''',
        a_legacy=LegacyDefinition(name='data_repository_name'))

    data_repository_url = Quantity(
        type=str,
        shape=[],
        description='''
        An URL to the repository, where the data is stored.
        ''',
        a_legacy=LegacyDefinition(name='data_repository_url'))

    data_preview_url = Quantity(
        type=str,
        shape=[],
        description='''
        An URL to an image file that contains a preview.
        ''',
        a_legacy=LegacyDefinition(name='data_preview_url'))

    entry_repository_url = Quantity(
        type=str,
        shape=[],
        description='''
        An URL to the entry on the repository, where the data is stored.
        ''',
        a_legacy=LegacyDefinition(name='entry_repository_url'))


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


m_package.__init_metainfo__()
