import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition


m_package = Package(
    name='general_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='general.nomadmetainfo.json'))


class section_entry_info(MSection):
    '''
    General information about this entry that is independent from its domain, field, or
    used parser
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_entry_info'))

    entry_upload_time = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        Upload datetime, given as total number of seconds is the elapsed since the unix
        epoch (1 January 1970)
        ''',
        a_legacy=LegacyDefinition(name='entry_upload_time'))

    entry_uploader_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the uploader, given as lastname, firstname.
        ''',
        a_legacy=LegacyDefinition(name='entry_uploader_name'))

    entry_uploader_id = Quantity(
        type=str,
        shape=[],
        description='''
        The id of the uploader.
        ''',
        a_legacy=LegacyDefinition(name='entry_uploader_id'))

    upload_id = Quantity(
        type=str,
        shape=[],
        description='''
        Nomad upload id
        ''',
        a_legacy=LegacyDefinition(name='upload_id'))

    calc_id = Quantity(
        type=str,
        shape=[],
        description='''
        Nomad calc id.
        ''',
        a_legacy=LegacyDefinition(name='calc_id'))

    calc_hash = Quantity(
        type=str,
        shape=[],
        description='''
        Calculation hash based on raw file contents.
        ''',
        a_legacy=LegacyDefinition(name='calc_hash'))

    mainfile = Quantity(
        type=str,
        shape=[],
        description='''
        Path to the main file within the upload.
        ''',
        a_legacy=LegacyDefinition(name='mainfile'))

    parser_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the parser used to extract this information.
        ''',
        a_legacy=LegacyDefinition(name='parser_name'))

    filepaths = Quantity(
        type=str,
        shape=['number_of_files'],
        description='''
        Filepaths of files that belong to this entry, i.e. files in the same directory.
        Filepaths are relative to the upload.
        ''',
        a_legacy=LegacyDefinition(name='filepaths'))

    number_of_files = Quantity(
        type=int,
        shape=[],
        description='''
        Number of files that belong to this entry.
        ''',
        a_legacy=LegacyDefinition(name='number_of_files'))

    section_archive_processing_info = SubSection(
        sub_section=SectionProxy('section_archive_processing_info'),
        repeats=True,
        a_legacy=LegacyDefinition(name='section_archive_processing_info'))


class section_archive_processing_info(MSection):
    '''
    Information about the used archive processing steps and their execution.
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='section_archive_processing_info'))

    archive_processor_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the applied archive processing program.
        ''',
        a_legacy=LegacyDefinition(name='archive_processor_name'))

    archive_processor_error = Quantity(
        type=str,
        shape=[],
        description='''
        The main error during execution of the archive processing program that failed the
        program.
        ''',
        a_legacy=LegacyDefinition(name='archive_processor_error'))

    number_of_archive_processor_warnings = Quantity(
        type=int,
        shape=[],
        description='''
        Number of warnings during execution of the archive processing program.
        ''',
        a_legacy=LegacyDefinition(name='number_of_archive_processor_warnings'))

    archive_processor_warnings = Quantity(
        type=str,
        shape=['number_of_archive_processor_warnings'],
        description='''
        Warnings during execution of the archive processing program.
        ''',
        a_legacy=LegacyDefinition(name='archive_processor_warnings'))

    archive_processor_status = Quantity(
        type=str,
        shape=[],
        description='''
        Status returned by archive processing program.
        ''',
        a_legacy=LegacyDefinition(name='archive_processor_status'))


m_package.__init_metainfo__()
