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

import pytest
import os
import os.path
import re
import datetime
import yaml

from nomad import config
from nomad.datamodel.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.context import ClientContext
from nomad.utils import generate_entry_id, strip
from nomad.parsing.parser import ArchiveParser
from tests.normalizing.conftest import run_normalize
from nomad.processing import Upload, Entry
from nomad.processing import ProcessStatus
from nomad import files


def quantity_generator(quantity_name, header_name, shape='shape: [\'*\']', to_dict=False):
    if to_dict:
        return {quantity_name: {
            'type': 'str',
            'shape': shape,
            'm_annotations': {
                'tabular': {
                    'name': header_name
                }
            }
        }}
    else:
        base_case = f'''{quantity_name}:
                type: str
                {shape}
                m_annotations:
                  tabular:
                    name: {header_name}'''
        return re.sub(r'\n\s*\n', '\n', base_case)


@pytest.mark.parametrize('schema', [
    pytest.param(
        strip('''
            definitions:
                name: 'testing_modes'
                sections:
                    My_schema:
                        base_section: nomad.datamodel.data.EntryData
                        sub_sections:
                            process:
                                section:
                                    base_section: nomad.parsing.tabular.TableData
                                    quantities:
                                        data_file:
                                            type: str
                                            m_annotations:
                                                tabular_parser:
                                                    parsing_options:
                                                        comment: '#'
                                        quantity_1:
                                            type: str
                                            m_annotations:
                                                tabular:
                                                    name: column_1
            data:
                m_def: My_schema
                process:
                    data_file: Test.xlsx
        '''), id='w/o_sheetName_colMode'),
    pytest.param(
        strip('''
            definitions:
                name: 'testing_modes'
                sections:
                    My_schema:
                        base_section: nomad.datamodel.data.EntryData
                        sub_sections:
                            process:
                                section:
                                    base_section: nomad.parsing.tabular.TableData
                                    quantities:
                                        data_file:
                                            type: str
                                            m_annotations:
                                                tabular_parser:
                                                    parsing_options:
                                                        comment: '#'
                                        quantity_1:
                                            type: str
                                            m_annotations:
                                                tabular:
                                                    name: sheet_1/column_1
            data:
                m_def: My_schema
                process:
                    data_file: Test.xlsx
        '''), id='w_sheetName_colMode'),
    pytest.param(
        strip('''
            definitions:
                name: 'testing_modes'
                sections:
                    My_schema:
                        base_section: nomad.datamodel.data.EntryData
                        sub_sections:
                            process:
                                section:
                                    base_section: nomad.parsing.tabular.TableData
                                    quantities:
                                        data_file:
                                            type: str
                                            description: |
                                                A reference to an uploaded .xlsx
                                            m_annotations:
                                                tabular_parser:
                                                    parsing_options:
                                                        comment: '#'
                                        quantity_1:
                                            type: str
                                            m_annotations:
                                                tabular:
                                                    name: sheet_1/column_1
                                        quantity_2:
                                            type: np.float64
                                            shape: ['*']
                                            unit: K
                                            m_annotations:
                                                tabular:
                                                    name: sheet_2/column_2
            data:
                m_def: My_schema
                process:
                    data_file: Test.xlsx
        '''), id='w_sheetName_colMode_multy_quantities'),
    pytest.param(
        strip('''
            definitions:
                name: 'testing_similar_cols_in_datafile'
                sections:
                    My_schema:
                        base_section: nomad.parsing.tabular.TableData
                        quantities:
                            data_file:
                                type: str
                                m_annotations:
                                    tabular_parser:
                                        parsing_options:
                                            comment: '#'
                                        mapping_options:
                                            -   mapping_mode: row
                                                file_mode: current_entry
                                                sections:
                                                    -  test_subsection
                        sub_sections:
                            test_subsection:
                                repeats: true
                                section:
                                    sub_sections:
                                        my_subsection:
                                            repeats: true
                                            section:
                                                quantities:
                                                    my_quantity_1:
                                                        type: str
                                                        m_annotations:
                                                            tabular:
                                                                name: sheet_3/quantity_1
                                                    my_quantity_2:
                                                        type: str
                                                        m_annotations:
                                                            tabular:
                                                                name: sheet_3/quantity_2
            data:
                m_def: My_schema
                data_file: Test.xlsx
        '''), id='row mode with similar multiple columns in the excel sheet'),
    pytest.param(
        strip('''
        definitions:
            name: 'testing_mixed_mode'
            sections:
                MyColumnMode:
                    quantities:
                        my_quantity_1:
                            type: np.float64
                            shape: ['*']
                            m_annotations:
                                tabular:
                                    name: sheet_4/Header 1
                    sub_sections:
                        MyColumnModeSubsection:
                            section:
                                quantities:
                                    my_quantity_2:
                                        type: str
                                        shape: ['*']
                                        m_annotations:
                                            tabular:
                                                name: sheet_4/Header 2
                                    my_quantity_3:
                                        type: str
                                        shape: ['*']
                                        m_annotations:
                                            tabular:
                                                name: sheet_4/Header 3
                MyRowMode:
                    sub_sections:
                        MyRowModeSubsection:
                            repeats: true
                            section:
                                quantities:
                                    my_quantity_4:
                                        type: str
                                        m_annotations:
                                            tabular:
                                                name: sheet_4/Header 4
                MyMixedMode:
                    base_sections:
                        - nomad.parsing.tabular.TableData
                    quantities:
                        data_file:
                            type: str
                            default: datafile.xlsx
                            m_annotations:
                                tabular_parser:
                                    parsing_options:
                                        comment: '#'
                                    mapping_options:
                                        -   mapping_mode: row
                                            file_mode: current_entry
                                            sections:
                                                -  MySubsection2
                                        -   mapping_mode: column
                                            file_mode: current_entry
                                            sections:
                                                -  MySubsection1
                    sub_sections:
                        MySubsection1:
                            section: '#/MyColumnMode'
                        MySubsection2:
                            repeats: true
                            section:
                                base_section: '#/MyRowMode'
        data:
            m_def: MyMixedMode
            data_file: Test.xlsx
    '''), id='mixed mode')
])
def test_tabular_complex_schema(raw_files, monkeypatch, schema):
    '''
    Testing TabularParser parser. This feature creates an entry out of each row from the given excel/csv file
    '''
    _, schema_file = get_files(schema)
    excel_file = os.path.join(os.path.dirname(__file__), '../../tests/data/parsers/tabular/Test.xlsx')

    class MyContext(ClientContext):
        def raw_file(self, path, *args, **kwargs):
            return open(excel_file, *args, **kwargs)
    context = MyContext(local_dir='')

    main_archive, _ = get_archives(context, schema_file, None)
    ArchiveParser().parse(schema_file, main_archive)
    run_normalize(main_archive)

    assert main_archive.data is not None
    if 'testing_modes' in schema:
        assert 'quantity_1' in main_archive.data.process
        assert main_archive.data.process.quantity_1 == 'value_1'
        if 'quantity_2' in main_archive.data.process:
            assert len(main_archive.data.process['quantity_2']) == 6
    elif 'testing_similar_cols_in_datafile' in schema:
        assert len(main_archive.data.test_subsection) == 2  # 2 rows in sheet_3 of the Excel file
        for row_index, test_subsection in enumerate(main_archive.data.test_subsection):
            for data_index, my_subsection in enumerate(test_subsection.my_subsection):
                assert my_subsection['my_quantity_1'] == f'q1_d{data_index}_r{row_index}'
                assert my_subsection['my_quantity_2'] == f'q2_d{data_index}_r{row_index}'
    elif 'testing_mixed_mode' in schema:
        assert len(main_archive.data.MySubsection2) == 3
        assert main_archive.data.MySubsection2[0].MyRowModeSubsection[0].my_quantity_4 == 'my_str1'
        assert len(main_archive.data.MySubsection1.my_quantity_1) == 3
        assert len(main_archive.data.MySubsection1.MyColumnModeSubsection.my_quantity_2) == 3
        assert len(main_archive.data.MySubsection1.MyColumnModeSubsection.my_quantity_3) == 3


def test_tabular_entry_mode(mongo, test_user, raw_files, monkeypatch, proc_infra):

    upload = Upload(
        upload_id='test_upload_id',
        main_author=test_user.user_id)
    upload.save()
    files.StagingUploadFiles(upload_id=upload.upload_id, create=True)
    upload.staging_upload_files.add_rawfiles('tests/data/parsers/tabular/')
    upload.process_upload()
    upload.block_until_complete()
    assert upload is not None
    assert upload.processed_entries_count == 2

    for entry in Entry.objects(upload_id='test_upload_id'):
        assert entry.process_status == ProcessStatus.SUCCESS


@pytest.mark.parametrize('test_case, sub_sections_placeholder, mapping_options_placeholder, csv_content', [
    pytest.param('test_1', {}, [], 'header_0,header_1\n0_0,0_1\n1_0,1_1', id='simple'),
    pytest.param(
        'test_2', {'my_substance': {'section': '#/Mysection'}},
        [{'mapping_mode': 'column', 'file_mode': 'current_entry', 'sections': ['my_substance']}],
        'header_0,header_1\n0_0,0_1\n1_0,1_1', id='nested'),
])
def test_tabular_column_mode(raw_files, monkeypatch, test_case, sub_sections_placeholder, mapping_options_placeholder,
                             csv_content):
    '''
    Testing the TableData normalizer using default mode (column mode). This feature creates a list of values
    out of the given column in the excel/csv file for the given quantity.
    '''
    base_schema = {
        'definitions': {
            'name': 'Eln',
            'sections': {
                'Mysection': {'quantities': quantity_generator('quantity_1', 'header_1', shape=['*'], to_dict=True)},
                'My_schema': {
                    'base_section': 'nomad.parsing.tabular.TableData',
                    'quantities': {
                        'data_file': {
                            'type': 'str', 'description': test_case,
                            'm_annotations': {
                                'tabular_parser': {
                                    'parsing_options': {
                                        'comment': '#'
                                    },
                                    'mapping_options': mapping_options_placeholder
                                }
                            }
                        },
                        'quantity_0': {'type': 'str', 'shape': ['*'], 'm_annotations': {'tabular': {'name': 'header_0'}}}
                    },
                    'sub_sections': sub_sections_placeholder
                }
            }
        },
        'data': {'m_def': 'My_schema', 'data_file': 'test.my_schema.archive.csv'}
    }
    schema = yaml.dump(base_schema)
    csv_file, schema_file = get_files(schema, csv_content)

    class MyContext(ClientContext):
        def raw_file(self, path, *args, **kwargs):
            return open(csv_file, *args, **kwargs)
    context = MyContext(local_dir='')

    main_archive, _ = get_archives(context, schema_file, None)
    ArchiveParser().parse(schema_file, main_archive)
    run_normalize(main_archive)

    assert main_archive.data is not None
    if 'test_1' in schema:
        assert main_archive.data.quantity_0 == ['0_0', '1_0']
    elif 'test_2' in schema:
        assert main_archive.data.my_substance.quantity_1 == ['0_1', '1_1']


@pytest.mark.parametrize('test_case, section_placeholder, target_sub_section_placeholder, sub_sections_placeholder, csv_content', [
    pytest.param('test_1', None, ['my_substance1'],
                 {'my_substance1': {'repeats': 'true', 'section': {'base_section': 'Substance1'}}},
                 'header_0,header_1\n0_0,0_1\n1_0,1_1',
                 id='simple_1_section'),
    pytest.param('test_2',
                 {'Substance2': {'quantities': quantity_generator('quantity_2', 'header_2', shape='', to_dict=True)}},
                 ['my_substance1', 'my_substance2'],
                 {
                     'my_substance1': {'repeats': 'true', 'section': {'base_section': 'Substance1'}},
                     'my_substance2': {'repeats': 'true', 'section': {'base_section': 'Substance2'}}
                 },
                 'header_0,header_1,header_2\n0_0,0_1,0_2\n1_0,1_1,1_2',
                 id='simple_2_sections'),
    pytest.param('test_3', None,
                 ['subsection_1/my_substance1'],
                 {'subsection_1': {'section': {'sub_sections': {
                     'my_substance1': {'repeats': 'true', 'section': {'base_section': 'Substance1'}}
                 }}}},
                 'header_0,header_1\n0_0,0_1\n1_0,1_1',
                 id='nested')])
def test_tabular_row_mode(raw_files, monkeypatch, test_case, section_placeholder, target_sub_section_placeholder,
                          sub_sections_placeholder, csv_content):
    '''
    Testing the TableData normalizer with mode set to row. This feature is used to create a section out of each row in a
    given sheet_name of an excel file or a csv file, and append it to the repeating (sub)section(s).
    '''
    base_schema = {
        'definitions': {
            'name': 'Eln',
            'sections': {
                'Substance1': {'quantities': quantity_generator('quantity_4', 'header_0', shape='', to_dict=True)},
                'Substance2': {'quantities': quantity_generator('quantity_2', 'header_2', shape='', to_dict=True)},
                'My_schema': {
                    'base_section': 'nomad.parsing.tabular.TableData',
                    'quantities': {'data_file': {
                        'type': 'str', 'description': test_case,
                        'm_annotations': {
                            'tabular_parser': {
                                'parsing_options': {
                                    'comment': '#'
                                },
                                'mapping_options': [
                                    {'mapping_mode': 'row', 'file_mode': 'current_entry', 'sections': target_sub_section_placeholder}
                                ]
                            }
                        }
                    }},
                    'sub_sections': sub_sections_placeholder
                }
            }
        },
        'data': {'m_def': 'My_schema', 'data_file': 'test.my_schema.archive.csv'}
    }
    schema = yaml.dump(base_schema)
    csv_file, schema_file = get_files(schema, csv_content)

    class MyContext(ClientContext):
        def raw_file(self, path, *args, **kwargs):
            return open(csv_file, *args, **kwargs)
    context = MyContext(local_dir='')

    main_archive, _ = get_archives(context, schema_file, None)
    ArchiveParser().parse(schema_file, main_archive)
    run_normalize(main_archive)

    assert main_archive.data is not None
    if 'test_1' in schema:
        assert len(main_archive.data.my_substance1) == 2
        ii = 0
        for item in main_archive.data.my_substance1:
            assert item.quantity_4 == f'{ii}_0'
            ii += 1
    elif 'test_2' in schema:
        assert len(main_archive.data.my_substance2) == 2
        ii = 0
        for item in main_archive.data.my_substance2:
            assert item.quantity_2 == f'{ii}_2'
            ii += 1
    elif 'test_3' in schema:
        assert len(main_archive.data.subsection_1.my_substance1) == 2
        ii = 0
        for item in main_archive.data.subsection_1.my_substance1:
            assert item.quantity_4 == f'{ii}_0'
            ii += 1


@pytest.mark.parametrize('schema,content', [
    pytest.param(
        strip('''
            definitions:
                name: 'space in header'
                sections:
                    MyTable:
                        base_section: nomad.parsing.tabular.TableData
                        quantities:
                            data_file:
                              type: str
                              m_annotations:
                                tabular_parser:
                                    parsing_options:
                                      comment: '#'
                            header_0:
                                type: str
                            header_1:
                                type: str
            data:
                m_def: MyTable
                data_file: test.my_schema.archive.csv
        '''),
        strip('''
            header_0, header_1
            a,b
        '''), id='space in header'
    ),
    pytest.param(
        strip('''
            definitions:
                name: 'datetime in row mode'
                sections:
                    MyTableWithDatetime:
                        base_section: nomad.parsing.tabular.TableData
                        m_annotations:
                            eln:
                        quantities:
                            data_file:
                                type: str
                                m_annotations:
                                    tabular_parser:
                                        parsing_options:
                                            comment: '#'
                                        mapping_options:
                                        -   mapping_mode: row
                                            file_mode: current_entry
                                            sections:
                                            -   MySubsection
                        sub_sections:
                            MySubsection:
                                repeats: true
                                section:
                                    quantities:
                                        header_0:
                                            type: Datetime
                                            m_annotations:
                                                tabular:
                                                    name: "Header_0"
            data:
                m_def: MyTableWithDatetime
                data_file: test.my_schema.archive.csv
        '''),
        strip('''
            Header_0, Header_1
            2018-03-19 13:01:47+01, 1
            2018-03-19 13:01:48+01, 2
            2018-03-19 13:01:49+01, 3
        '''), id='datetime in row mode'
    ),
    pytest.param(
        strip('''
        definitions:
          name: 'data file as reference'
          sections:
            Base_Class:
              quantities:
                name:
                  type: str
            Main_Class:
              base_sections:
                - nomad.parsing.tabular.TableData
              quantities:
                csv_file:
                  type: str
                  default: 'placeholder'
                  m_annotations:
                    tabular_parser:
                        parsing_options:
                            comment: '#'
                        mapping_options:
                        -   mapping_mode: row
                            file_mode: current_entry
                            sections:
                            -   subsection_1
              sub_sections:
                subsection_1:
                  repeats: True
                  section: '#/Class_1'
            Class_1:
              base_sections:
              - '#/Base_Class'
              sub_sections:
                subsection_1:
                  repeats: true
                  section:
                    quantities:
                      header_0:
                        type: str
                        m_annotations:
                          tabular:
                            name: header_0
        data:
          m_def: Main_Class
          csv_file: test.my_schema.archive.csv
          subsection_1:
          - m_def: Class_1
        '''),
        strip('''
            header_0,header_1
            0_0,0_1
            1_0,1_1
        '''), id='data file as reference'
    )
])
def test_tabular_csv(raw_files, monkeypatch, schema, content):
    '''Tests that missing data is handled correctly. Pandas by default
    interprets missing numeric values as NaN, which are incompatible with
    metainfo.
    '''
    csv_file, schema_file = get_files(schema, content)

    class MyContext(ClientContext):
        def raw_file(self, path, *args, **kwargs):
            return open(csv_file, *args, **kwargs)
    context = MyContext(local_dir='')

    main_archive, _ = get_archives(context, schema_file, None)
    ArchiveParser().parse(schema_file, main_archive)
    run_normalize(main_archive)
    if re.search('space in header', schema):
        assert main_archive.data.header_1 is not None
    elif re.search('datetime in row mode', schema):
        for instance in [0, 1, 2]:
            assert isinstance(main_archive.data.MySubsection[instance].header_0, datetime.datetime)
    elif re.search('data file as reference', schema):
        assert len(main_archive.data.subsection_1) == 2
        assert len(main_archive.data.subsection_1[0].subsection_1) == 1
        assert main_archive.data.subsection_1[0].subsection_1[0].header_0 == '0_0'
        assert main_archive.data.subsection_1[1].subsection_1[0].header_0 == '1_0'


@pytest.mark.parametrize('schema,content', [
    pytest.param(
        strip('''
        definitions:
            name: 'TableData checkbox'
            sections:
                Main_Class:
                    base_section: nomad.parsing.tabular.TableData
                    quantities:
                        data_file:
                            type: str
                            m_annotations:
                                tabular_parser:
                                    parsing_options:
                                        comment: '#'
                                    mapping_options:
                                    -   mapping_mode: row
                                        file_mode: current_entry
                                        sections:
                                        -   MySubsection
                    sub_sections:
                        MySubsection:
                            repeats: true
                            section:
                                quantities:
                                    header_0:
                                        type: str
                                        m_annotations:
                                            tabular:
                                                name: header_0
        data:
          m_def: Main_Class
          data_file: test.my_schema.archive.csv
        '''),
        strip('''
            header_0,header_1
            a,b
        '''), id='checkoing checkbox'
    ),
    pytest.param(
        strip('''
        definitions:
            name: 'hiding checkbox'
            sections:
                Main_Class:
                    base_section: nomad.parsing.tabular.TableData
                    m_annotations:
                        eln:
                            hide: ['fill_archive_from_datafile']
                    quantities:
                        data_file:
                            type: str
                            m_annotations:
                                tabular_parser:
                                    parsing_options:
                                        comment: '#'
                                    mapping_options:
                                    -   mapping_mode: row
                                        file_mode: current_entry
                                        sections:
                                        -   MySubsection
                    sub_sections:
                        MySubsection:
                            repeats: true
                            section:
                                quantities:
                                    header_0:
                                        type: str
                                        m_annotations:
                                            tabular:
                                                name: header_0
        data:
          m_def: Main_Class
          data_file: test.my_schema.archive.csv
        '''),
        strip('''
            header_0,header_1
            a,b
        '''), id='hide checkbox')
])
def test_tabular_checkbox(raw_files, monkeypatch, schema, content):
    '''Tests that missing data is handled correctly. Pandas by default
        interprets missing numeric values as NaN, which are incompatible with
        metainfo.
        '''
    csv_file, schema_file = get_files(schema, content)

    class MyContext(ClientContext):
        def raw_file(self, path, *args, **kwargs):
            return open(csv_file, *args, **kwargs)

    context = MyContext(local_dir='')

    main_archive, _ = get_archives(context, schema_file, None)
    ArchiveParser().parse(schema_file, main_archive)
    run_normalize(main_archive)

    if 'TableData checkbox' in schema:
        assert main_archive.data.MySubsection[0].header_0 == 'a'
        assert main_archive.data.fill_archive_from_datafile is False
        main_archive.data.MySubsection[0].header_0 = 'c'
        run_normalize(main_archive)
        assert main_archive.data.MySubsection[0].header_0 == 'c'
    elif 'fill_archive_from_datafile' in schema:
        assert main_archive.data.MySubsection[0].header_0 == 'a'
        assert main_archive.data.fill_archive_from_datafile is True
        main_archive.data.MySubsection[0].header_0 = 'c'
        run_normalize(main_archive)
        assert main_archive.data.MySubsection[0].header_0 == 'a'


def get_files(schema=None, content=None):
    '''Prepares files containing schema and data in the temporary file
    directory.
    '''
    if schema:
        schema_file = os.path.join(config.fs.tmp, 'my_schema.archive.yaml')
        with open(schema_file, 'wt') as f:
            f.write(schema)
    else:
        schema_file = None
    if content:
        mainfile = os.path.join(config.fs.tmp, 'test.my_schema.archive.csv')
        with open(mainfile, 'wt') as f:
            f.write(content)
    else:
        mainfile = None
    return mainfile, schema_file


def get_context(upload_id, schema_file):
    '''Prepares a custom context for testing.
    '''
    class MyContext(ClientContext):
        def load_raw_file(self, path: str, upload_id: str, installation_url: str, url: str = None):
            archive = super().load_raw_file(path, upload_id, installation_url, url)
            archive.metadata = EntryMetadata(
                upload_id=upload_id,
                entry_id=generate_entry_id(upload_id, schema_file))
            return archive

    return MyContext(local_dir='')


def get_archives(context, mainfile, upload_id, keys=None):
    '''Prepares a main archive and child archives for parsing.
    '''
    main_archive = EntryArchive(m_context=context, metadata=EntryMetadata(
        upload_id=None,
        mainfile=mainfile,
        entry_id=generate_entry_id(upload_id, mainfile)))
    if keys:
        child_archives = {
            key: EntryArchive(m_context=context, metadata=EntryMetadata(
                upload_id=upload_id,
                mainfile=mainfile,
                mainfile_key=key,
                entry_id=generate_entry_id(upload_id, mainfile, key)))
            for key in keys}
    else:
        child_archives = None

    return main_archive, child_archives
