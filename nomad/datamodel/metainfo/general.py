import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)


m_package = Package(name='general', description='None')


class section_entry_info(MSection):
    '''
    General information about this entry that is independent from its domain, field, or
    used parser
    '''

    m_def = Section(validate=False)

    entry_upload_time = Quantity(
        type=np.dtype(np.int64),
        shape=[],
        description='''
        Upload datetime, given as total number of seconds is the elapsed since the unix
        epoch (1 January 1970)
        ''')

    entry_uploader_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the uploader, given as lastname, firstname.
        ''')

    entry_uploader_id = Quantity(
        type=str,
        shape=[],
        description='''
        The id of the uploader.
        ''')

    upload_id = Quantity(
        type=str,
        shape=[],
        description='''
        Nomad upload id
        ''')

    calc_id = Quantity(
        type=str,
        shape=[],
        description='''
        Nomad calc id.
        ''')

    calc_hash = Quantity(
        type=str,
        shape=[],
        description='''
        Calculation hash based on raw file contents.
        ''')

    mainfile = Quantity(
        type=str,
        shape=[],
        description='''
        Path to the main file within the upload.
        ''')

    parser_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the parser used to extract this information.
        ''')

    filepaths = Quantity(
        type=str,
        shape=['number_of_files'],
        description='''
        Filepaths of files that belong to this entry, i.e. files in the same directory.
        Filepaths are relative to the upload.
        ''')

    number_of_files = Quantity(
        type=int,
        shape=[],
        description='''
        Number of files that belong to this entry.
        ''')

    section_archive_processing_info = SubSection(
        sub_section=SectionProxy('section_archive_processing_info'),
        repeats=True)


class section_archive_processing_info(MSection):
    '''
    Information about the used archive processing steps and their execution.
    '''

    m_def = Section(validate=False)

    archive_processor_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name of the applied archive processing program.
        ''')

    archive_processor_error = Quantity(
        type=str,
        shape=[],
        description='''
        The main error during execution of the archive processing program that failed the
        program.
        ''')

    number_of_archive_processor_warnings = Quantity(
        type=int,
        shape=[],
        description='''
        Number of warnings during execution of the archive processing program.
        ''')

    archive_processor_warnings = Quantity(
        type=str,
        shape=['number_of_archive_processor_warnings'],
        description='''
        Warnings during execution of the archive processing program.
        ''')

    archive_processor_status = Quantity(
        type=str,
        shape=[],
        description='''
        Status returned by archive processing program.
        ''')


m_package.__init_metainfo__()
