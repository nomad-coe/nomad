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
import json
import os
from zipfile import ZipFile

import pytest

from nomad.config import config
from nomad.app.v1.routers.metainfo import store_package_definition
from nomad.datamodel import EntryArchive, ClientContext
from nomad.metainfo import MSection, MetainfoReferenceError
from nomad.utils import generate_entry_id, create_uuid
from tests.processing.test_data import run_processing


@pytest.mark.parametrize(
    'metainfo_data',
    [
        pytest.param(
            {
                'm_def': 'nomad.metainfo.metainfo.Package',
                'name': 'test.Package',
                'section_definitions': [{'name': 'MySection'}],
            },
            id='python',
        )
    ],
)
def test_metainfo_section_id_endpoint(metainfo_data, mongo_module, client):
    assert (
        MSection.from_dict(metainfo_data).m_to_dict(
            with_root_def=True, with_out_meta=True
        )
        == metainfo_data
    )

    package = MSection.from_dict(metainfo_data)

    store_package_definition(package, with_root_def=True, with_out_meta=True)

    section_id = package.section_definitions[0].definition_id

    response = client.get(f'metainfo/{section_id}')
    assert response.status_code == 200
    pkg_definition = response.json()['data']
    del pkg_definition['entry_id_based_name']
    assert pkg_definition == metainfo_data

    response = client.get(f'metainfo/{section_id[::-1]}')
    assert response.status_code == 404


def test_upload_and_download(
    client, test_user, proc_infra, mongo_module, no_warn, monkeypatch, tmp
):
    monkeypatch.setattr('nomad.config.process.store_package_definition_in_mongo', True)
    monkeypatch.setattr('nomad.config.process.add_definition_id_to_reference', True)
    monkeypatch.setattr('nomad.config.process.write_definition_id_to_archive', True)

    monkeypatch.setattr('requests.get', getattr(client, 'get'))
    monkeypatch.setattr('requests.post', getattr(client, 'post'))

    m_def = '../upload/raw/schema.archive.json#/definitions/section_definitions/1'

    def client_context():
        return ClientContext(installation_url='')

    def join_tmp_dir(fn: str) -> str:
        return os.path.join(tmp, fn)

    schema_file_name = 'schema.archive.json'
    data_file_name = 'sample.archive.json'
    archive_name = 'example_versioned_metainfo.zip'

    def pack_and_publish(
        property_name: str, property_value: str, with_schema: bool, def_id: str = None
    ):
        """
        Generates and publishes an example upload with customizable schema and
        data using the schema with different m_def flavors.
        """
        schema_path = join_tmp_dir(schema_file_name)
        data_path = join_tmp_dir(data_file_name)
        archive_path = join_tmp_dir(archive_name)

        if with_schema:
            with open(schema_path, 'w') as f:
                json.dump(
                    {
                        'name': 'test schema package',
                        'definitions': {
                            'section_definitions': [
                                {
                                    'base_sections': ['nomad.datamodel.data.EntryData'],
                                    'name': 'Chemical',
                                },
                                {
                                    'base_sections': ['nomad.datamodel.data.EntryData'],
                                    'name': 'Sample',
                                    'quantities': [
                                        {
                                            'name': property_name,
                                            'type': {
                                                'type_kind': 'python',
                                                'type_data': 'str',
                                            },
                                        }
                                    ],
                                },
                            ]
                        },
                    },
                    f,
                )

        with open(data_path, 'w') as f:
            use_m_def = m_def
            if def_id:
                use_m_def = f'{use_m_def}@{def_id}'
            json.dump({'data': {'m_def': use_m_def, property_name: property_value}}, f)

        with ZipFile(archive_path, 'w') as zipObj:
            if with_schema:
                zipObj.write(schema_path, arcname=schema_file_name)
            zipObj.write(data_path, arcname=data_file_name)

        processed = run_processing(
            (create_uuid(), archive_path), test_user, publish_directly=True
        )
        return processed.upload_id

    # 1. create a first upload
    upload_1_id = pack_and_publish(
        property_name='test_quantity',
        property_value='test_value',
        with_schema=True,
        def_id=None,
    )

    response = client.get(
        f'entries/{generate_entry_id(upload_1_id, data_file_name)}/archive'
    )
    assert response.status_code == 200, response.ext
    entry_data = response.json()['data']['archive']['data']

    # check if 'test_quantity' quantity is in the entry
    assert 'test_quantity' in entry_data
    assert entry_data['test_quantity'] == 'test_value'

    # 2. check if package is stored in mongo
    original_def_id = entry_data['m_def_id']
    response = client.get(f'metainfo/{original_def_id}')
    assert response.status_code == 200

    # 3. prepare a new entry refers to the previously uploaded package
    upload_2_id = pack_and_publish(
        property_name='test_quantity',
        property_value='new_value',
        with_schema=False,
        def_id=entry_data['m_def_id'],
    )

    response = client.get(
        f'entries/{generate_entry_id(upload_2_id, data_file_name)}/archive'
    )
    assert response.status_code == 200, response.text
    entry_data = response.json()['data']['archive']['data']

    # 4. check if 'test_quantity' quantity is in the entry and has the correct value
    assert 'test_quantity' in entry_data
    assert entry_data['test_quantity'] == 'new_value'

    # 5. test if client side can read the package using versioned package
    entry = EntryArchive.m_from_dict(entry_data, m_context=client_context())
    assert entry.test_quantity == 'new_value'

    # 6. test if client side can detect wrong package version
    definition_reference, definition_id = entry_data['m_def'].split('@')
    entry_data['m_def'] = f'{definition_reference}@{definition_id[::-1]}'
    with pytest.raises(MetainfoReferenceError):
        EntryArchive.m_from_dict(entry_data, m_context=client_context())

    # 7. now test if client side can read the package using non-versioned package
    entry_data[
        'm_def'
    ] = f'/upload/{upload_1_id}/raw/{schema_file_name}#/definitions/section_definitions/1'
    entry_data = EntryArchive.m_from_dict(entry_data, m_context=client_context())
    assert entry_data.test_quantity == 'new_value'

    # 8. generate version two with 'updated_quantity' quantity
    # TODO this test is kinda pointless, because it does not test different version, but a
    # fully new schema
    upload_3_id = pack_and_publish(
        property_name='updated_quantity',
        property_value='new_value',
        with_schema=True,
        def_id=None,
    )
    response = client.get(
        f'entries/{generate_entry_id(upload_3_id, data_file_name)}/archive'
    )
    entry_data = response.json()['data']['archive']['data']

    # check if 'updated_quantity' quantity is in the entry
    assert 'updated_quantity' in entry_data

    # check two sections shall have different id
    assert entry_data['m_def_id'] != original_def_id


@pytest.fixture(scope='function')
def example_upload_two_schemas():
    return {
        'schema_1': {
            'name': 'test schema package',
            'definitions': {
                'section_definitions': [
                    {
                        'base_sections': ['nomad.datamodel.data.EntryData'],
                        'name': 'Chemical',
                    },
                    {
                        'base_sections': ['nomad.datamodel.data.EntryData'],
                        'name': 'Sample',
                    },
                ]
            },
        },
        'schema_2': {
            'definitions': {
                'section_definitions': [
                    {
                        'base_sections': [
                            '../upload/raw/schema_1.archive.json#/definitions/section_definitions/1'
                        ],
                        'name': 'Sample',
                        'quantities': [
                            {
                                'name': 'chemical',
                                'type': {
                                    'type_kind': 'reference',
                                    'type_data': '../upload/raw/schema_1.archive.json#/definitions/section_definitions/0',
                                },
                            }
                        ],
                    }
                ]
            }
        },
        'chemical': {
            'data': {
                'm_def': '../upload/raw/schema_1.archive.json#/definitions/section_definitions/0',
                'name': 'NaCl',
            }
        },
        'sample': {
            'data': {
                'm_def': '../upload/raw/schema_2.archive.json#/definitions/section_definitions/0',
                'name': 'MySample',
                'chemical': '../upload/raw/chemical.archive.json#/data',
            }
        },
    }


def test_two_schemas(
    example_upload_two_schemas, client, test_user, proc_infra, no_warn, monkeypatch
):
    monkeypatch.setattr('nomad.config.process.store_package_definition_in_mongo', True)
    monkeypatch.setattr('nomad.config.process.add_definition_id_to_reference', True)
    monkeypatch.setattr('nomad.config.process.write_definition_id_to_archive', True)

    def tmp(fn: str) -> str:
        return os.path.join(config.fs.tmp, fn)

    def public(fn: str) -> str:
        return os.path.join(
            config.fs.public, f"ex/{fn.replace('.zip', '')}/raw-public.plain.zip"
        )

    archive_name = 'example_upload_two_schemas.zip'

    # 1. pack and process initial archive containing two schemas
    with ZipFile(tmp(archive_name), 'w') as zipObj:
        for k, v in example_upload_two_schemas.items():
            zipObj.writestr(f'{k}.archive.json', json.dumps(v))

    processed = run_processing(
        (archive_name.replace('.zip', ''), tmp(archive_name)),
        test_user,
        publish_directly=True,
    )

    # 2. manually remove schema files
    with ZipFile(public(archive_name), 'w') as zipObj:
        for k, v in example_upload_two_schemas.items():
            if 'schema' in k:
                continue
            zipObj.writestr(f'{k}.archive.json', json.dumps(v))

    # retrieve the archive
    response = client.get(
        f'entries/{generate_entry_id(processed.upload_id, "sample.archive.json")}/archive'
    )

    entry_data = response.json()['data']['archive']['data']
    entry_data[
        'chemical'
    ] = f'/entries/{generate_entry_id(processed.upload_id, "chemical.archive.json")}/archive#/data'
    # get data in absence of original schema files
    entry = EntryArchive.m_from_dict(entry_data, m_context=ClientContext())

    assert entry.m_def.name == 'Sample'

    chemical = entry.chemical.m_proxy_resolve()

    assert chemical.m_def.name == 'Chemical'
