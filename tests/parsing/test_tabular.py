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
import pandas as pd

from nomad import config
from nomad.datamodel.datamodel import EntryArchive, EntryMetadata
from nomad.datamodel.context import ClientContext
from nomad.utils import generate_entry_id, strip
from nomad.parsing.tabular import TabularDataParser
from nomad.parsing.parser import ArchiveParser
from tests.normalizing.conftest import run_normalize


@pytest.mark.parametrize('schema,content', [
    pytest.param(
        strip('''
            definitions:
                sections:
                    MyTable:
                        base_section: nomad.parsing.tabular.TableRow
                        quantities:
                            header_0:
                                type: str
                            header_1:
                                type: str
        '''),
        strip('''
            header_0,header_1
            0_0,0_1
            1_0,1_1
        '''), id='simple'),
    pytest.param(
        strip('''
            definitions:
                sections:
                    MyTable:
                        base_section: nomad.parsing.tabular.TableRow
                        quantities:
                            header_0:
                                type: str
                            quantity:
                                type: str
                                m_annotations:
                                    tabular:
                                        name: header_1
                        sub_sections:
                            my_sub_section:
                                sub_section: MySubSection
                    MySubSection:
                        quantities:
                            quantity:
                                type: str
        '''),
        strip('''
            header_0,header_1,my_sub_section.quantity
            0_0,0_1,0_2
            1_0,1_1,1_2
        '''), id='nested'),
    pytest.param(
        strip('''
            definitions:
                sections:
                    MyTable:
                        base_section: nomad.parsing.tabular.TableRow
                        quantities:
                            header_0:
                                type: np.float64
                            header_1:
                                type: np.float64
                                unit: s
                                m_annotations:
                                    tabular:
                                        unit: ms
        '''),
        strip('''
            header_0,header_1
            0.0,0.1
            1.0,1.1
        '''), id='units')
])
def test_tabular(raw_files, monkeypatch, schema, content):
    mainfile = os.path.join(config.fs.tmp, 'test.my_schema.archive.csv')
    schema_file = os.path.join(config.fs.tmp, 'my_schema.archive.yaml')
    with open(mainfile, 'wt') as f:
        f.write(content)
    with open(schema_file, 'wt') as f:
        f.write(schema)

    data = pd.read_csv(mainfile)

    parser = TabularDataParser()
    keys = parser.is_mainfile(mainfile, 'text/application', bytes(), '')

    assert isinstance(keys, list)
    assert len(keys) == data.shape[0]

    class MyContext(ClientContext):
        def load_raw_file(self, path, upload_id, installation_url):
            archive = super().load_raw_file(path, upload_id, installation_url)
            archive.metadata = EntryMetadata(
                upload_id='test_upload',
                entry_id=generate_entry_id('test_upload', schema_file))
            return archive

    context = MyContext(local_dir='')
    main_archive = EntryArchive(m_context=context, metadata=EntryMetadata(
        upload_id=None,
        mainfile=mainfile,
        entry_id=generate_entry_id('test_upload', mainfile)))
    child_archives = {
        key: EntryArchive(m_context=context, metadata=EntryMetadata(
            upload_id='test_upload',
            mainfile=mainfile,
            mainfile_key=key,
            entry_id=generate_entry_id('test_upload', mainfile, key)))
        for key in keys}

    parser.parse(mainfile, main_archive, None, child_archives)
    main_archive.metadata.upload_id = 'test_upload'

    assert main_archive.data is not None
    for child_archive in child_archives.values():
        child_archive.data is not None

    # print('# main: ', json.dumps(main_archive.m_to_dict(), indent=2))
    # for key in keys:
    #     print(f'# {key}: ', json.dumps(child_archives[key].m_to_dict(), indent=2))


@pytest.mark.parametrize('schema', [
    pytest.param(
        strip('''
            definitions:
              name: 'A test schema for excel file parsing'
              sections:  
                MovpeSto_schema:  
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
                                comment: '#'
                              browser:
                                adaptor: RawFileAdaptor 
                              eln:
                                component: FileEditQuantity  
                          experiment_identifier:
                            type: str
                            m_annotations:
                              tabular:
                                name: Experiment Identifier
                              eln:
                                component: StringEditQuantity
            data:  
              m_def: MovpeSto_schema  
              process:
                data_file: Test.xlsx
        '''), id='w/o_sheetName_rowMode'),
    pytest.param(
        strip('''
            definitions:
              name: 'A test schema for excel file parsing'
              sections:  
                MovpeSto_schema:  
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
                                comment: '#'
                              browser:
                                adaptor: RawFileAdaptor 
                              eln:
                                component: FileEditQuantity  
                          experiment_identifier:
                            type: str
                            m_annotations:
                              tabular:
                                name: Overview/Experiment Identifier
                              eln:
                                component: StringEditQuantity
            data:  
              m_def: MovpeSto_schema  
              process:
                data_file: Test.xlsx
        '''), id='w_sheetName_rowMode'),
    pytest.param(
        strip('''
            definitions:
              name: 'A test schema for excel file parsing'
              sections:  
                MovpeSto_schema:  
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
                                comment: '#'
                              browser:
                                adaptor: RawFileAdaptor 
                              eln:
                                component: FileEditQuantity  
                          experiment_identifier:
                            type: str
                            m_annotations:
                              tabular:
                                name: Overview/Experiment Identifier
                              eln:
                                component: StringEditQuantity
                          pyrotemperature:
                            type: np.float64
                            shape: ['*']
                            unit: K
                            description: My test description here
                            m_annotations:
                              tabular:
                                name: Deposition Control/Pyrotemperature
            data:  
              m_def: MovpeSto_schema  
              process:
                data_file: Test.xlsx
        '''), id='w_sheetName_colMode')
])
def test_xlsx_tabular(raw_files, monkeypatch, schema):
    schema_file = os.path.join(config.fs.tmp, 'excel_parser.archive.yaml')
    with open(schema_file, 'wt') as f:
        f.write(schema)

    excel_file = os.path.join(os.path.dirname(__file__), '../../tests/data/parsers/tabular/Test.xlsx')

    class MyContext(ClientContext):
        def raw_file(self, path, *args, **kwargs):
            return open(excel_file, *args, **kwargs)
    context = MyContext(local_dir='')
    main_archive = EntryArchive(m_context=context, metadata=EntryMetadata(
        upload_id=None,
        mainfile=schema_file,
        entry_id=generate_entry_id('test_upload', schema_file)))
    ArchiveParser().parse(schema_file, main_archive)
    run_normalize(main_archive)

    assert main_archive.data is not None
    assert 'experiment_identifier' in main_archive.data.process
    assert main_archive.data.process.experiment_identifier == '22-01-21-MA-255'
    if 'pyrotemperature' in main_archive.data.process:
        assert len(main_archive.data.process['pyrotemperature']) == 6
