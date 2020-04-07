import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
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

    section_data = SubSection(
        sub_section=SectionProxy('section_data'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_data'))


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


m_package.__init_metainfo__()
