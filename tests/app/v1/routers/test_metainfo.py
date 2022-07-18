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

from nomad import config
from nomad.app.v1.routers.metainfo import store_package_definition
from nomad.datamodel import EntryArchive, ClientContext
from nomad.metainfo import MSection, MetainfoReferenceError
from nomad.utils import generate_entry_id
from tests.processing.test_data import run_processing


@pytest.mark.parametrize('metainfo_data', [
    pytest.param({
        'm_def': 'nomad.metainfo.metainfo.Package',
        'name': 'test.Package',
        'section_definitions': [
            {
                'name': 'MySection'
            }
        ]
    }, id='python')
])
def test_metainfo_section_id_endpoint(metainfo_data, mongo_infra, client):
    assert MSection.from_dict(metainfo_data).m_to_dict(with_root_def=True, with_out_meta=True) == metainfo_data

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


def simple_schema(name: str):
    return {
        "name": "test schema package",
        "definitions": {
            "section_definitions": [
                {
                    "base_sections": [
                        "nomad.datamodel.data.EntryData"
                    ],
                    "name": "Chemical"
                },
                {
                    "base_sections": [
                        "nomad.datamodel.data.EntryData"
                    ],
                    "name": "Sample",
                    "quantities": [
                        {
                            "name": name,
                            "type": {
                                "type_kind": "python",
                                "type_data": "str"
                            }
                        }
                    ]
                }
            ]
        }
    }


def simple_data(name: str):
    return {
        "data": {
            "m_def": "../upload/raw/schema.json#/definitions/section_definitions/1",
            name: "this is my name"
        }
    }


def test_upload_and_download(client, test_user, proc_infra, mongo_infra, no_warn, monkeypatch, tmp):
    monkeypatch.setattr('nomad.config.process.store_package_definition_in_mongo', True)
    monkeypatch.setattr('nomad.config.process.add_definition_id_to_reference', True)
    monkeypatch.setattr('nomad.config.process.write_definition_id_to_archive', True)

    def j(fn: str) -> str:
        return os.path.join(tmp, fn)

    schema_file_name = 'schema.json'
    data_file_name = 'sample.archive.json'
    archive_name = 'example_versioned_metainfo.zip'

    # 1. generate version one with 'chemicals' quantity
    # 2. upload and record version one
    def pack_and_publish(name: str):
        jschema = j(schema_file_name)
        jdata = j(data_file_name)
        jarchive = j(archive_name)

        with open(jschema, 'w') as f:
            json.dump(simple_schema(name), f)

        with open(jdata, 'w') as f:
            json.dump(simple_data(name), f)

        with ZipFile(jarchive, 'w') as zipObj:
            zipObj.write(jschema, arcname=schema_file_name)
            zipObj.write(jdata, arcname=data_file_name)

        return run_processing((name, jarchive), test_user, publish_directly=True)

    processed = pack_and_publish('chemicals')

    upload_id = processed.upload_id

    response = client.get(f'entries/{generate_entry_id(upload_id, data_file_name)}/archive')

    entry_data = response.json()['data']['archive']['data']

    # check if 'chemicals' quantity is in the entry
    assert 'chemicals' in entry_data

    # check if package is stored in mongo
    response = client.get(f'metainfo/{entry_data["m_def_id"]}')

    assert response.status_code == 200

    # 3. prepare a new entry refers to the previously uploaded package
    data_file_name = 'new_' + data_file_name
    with open(j(data_file_name), 'w') as f:
        data = simple_data('chemicals')
        data['data']['m_def'] += f'@{entry_data["m_def_id"]}'
        data['data']['chemicals'] = 'this is my new name'
        json.dump(data, f)

    processed = run_processing(
        (data_file_name.replace('.json', ''), j(data_file_name)), test_user,
        publish_directly=True)

    response = client.get(f'entries/{generate_entry_id(processed.upload_id, data_file_name)}/archive')

    new_entry_data = response.json()['data']['archive']['data']

    # 4. check if 'chemicals' quantity is in the entry and has the correct value
    assert 'chemicals' in new_entry_data
    assert new_entry_data['chemicals'] == 'this is my new name'

    # 5. test if client side can read the package using versioned package
    new_entry_data = EntryArchive.m_from_dict(data['data'], m_context=ClientContext())
    assert new_entry_data.chemicals == 'this is my new name'

    # 6. test if client side can detect wrong package version
    definition_reference, definition_id = data['data']['m_def'].split('@')
    data['data']['m_def'] = f'{definition_reference}@{definition_id[::-1]}'
    with pytest.raises(MetainfoReferenceError):
        EntryArchive.m_from_dict(data['data'], m_context=ClientContext())

    # 7. now test if client side can read the package using non-versioned package
    data['data']['m_def'] = f'/upload/{upload_id}/raw/schema.json#/definitions/section_definitions/1'
    new_entry_data = EntryArchive.m_from_dict(data['data'], m_context=ClientContext())
    assert new_entry_data.chemicals == 'this is my new name'

    # 8. generate version two with 'toxicchemicals' quantity
    processed = pack_and_publish('toxicchemicals')
    response = client.get(f'entries/{generate_entry_id(processed.upload_id, data_file_name)}/archive')

    new_entry_data = response.json()['data']['archive']['data']

    # check if 'chemicals' quantity is in the entry
    assert 'toxicchemicals' in new_entry_data

    # check two sections shall have different id
    assert entry_data["m_def_id"] != new_entry_data["m_def_id"]


@pytest.fixture(scope='function')
def example_upload_two_schemas():
    return {
        "schema_1": {
            "name": "test schema package",
            "definitions": {
                "section_definitions": [
                    {
                        "base_sections": [
                            "nomad.datamodel.data.EntryData"
                        ],
                        "name": "Chemical"
                    },
                    {
                        "base_sections": [
                            "nomad.datamodel.data.EntryData"
                        ],
                        "name": "Sample"
                    }
                ]
            }
        },
        "schema_2": {
            "definitions": {
                "section_definitions": [
                    {
                        "base_sections": [
                            "../upload/raw/schema_1.archive.json#/definitions/section_definitions/1"
                        ],
                        "name": "Sample",
                        "quantities": [
                            {
                                "name": "chemical",
                                "type": {
                                    "type_kind": "reference",
                                    "type_data": "../upload/raw/schema_1.archive.json#/definitions/section_definitions/0"
                                }
                            }
                        ]
                    }
                ]
            }
        },
        "chemical": {
            "data": {
                "m_def": "../upload/raw/schema_1.archive.json#/definitions/section_definitions/0",
                "name": "NaCl",
            }
        },
        "sample": {
            "data": {
                "m_def": "../upload/raw/schema_2.archive.json#/definitions/section_definitions/0",
                "name": "MySample",
                "chemical": "../upload/raw/chemical.archive.json#/data"
            }
        }
    }


def test_two_schemas(
        example_upload_two_schemas, client, test_user, proc_infra, mongo_infra, no_warn, monkeypatch, raw_files_infra):
    monkeypatch.setattr('nomad.config.process.store_package_definition_in_mongo', True)
    monkeypatch.setattr('nomad.config.process.add_definition_id_to_reference', True)
    monkeypatch.setattr('nomad.config.process.write_definition_id_to_archive', True)

    def tmp(fn: str) -> str:
        return os.path.join(config.fs.tmp, fn)

    def public(fn: str) -> str:
        return os.path.join(config.fs.public, f"ex/{fn.replace('.zip', '')}/raw-public.plain.zip")

    archive_name = 'example_upload_two_schemas.zip'

    # 1. pack and process initial archive containing two schemas
    with ZipFile(tmp(archive_name), 'w') as zipObj:
        for k, v in example_upload_two_schemas.items():
            zipObj.writestr(f'{k}.archive.json', json.dumps(v))

    processed = run_processing(
        (archive_name.replace('.zip', ''), tmp(archive_name)), test_user,
        publish_directly=True)

    # 2. manually remove schema files
    with ZipFile(public(archive_name), 'w') as zipObj:
        for k, v in example_upload_two_schemas.items():
            if 'schema' in k:
                continue
            zipObj.writestr(f'{k}.archive.json', json.dumps(v))

    # retrieve the archive
    response = client.get(f'entries/{generate_entry_id(processed.upload_id, "sample.archive.json")}/archive')

    entry_data = response.json()['data']['archive']['data']
    entry_data['chemical'] = f'/entries/{generate_entry_id(processed.upload_id, "chemical.archive.json")}/archive#/data'
    # get data in absence of original schema files
    entry = EntryArchive.m_from_dict(entry_data, m_context=ClientContext())

    assert entry.m_def.name == 'Sample'

    chemical = entry.chemical.m_proxy_resolve()

    assert chemical.m_def.name == 'Chemical'
