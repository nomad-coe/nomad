# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest
from datetime import datetime

from nomad import infrastructure, config
from nomad.datamodel import EntryMetadata
from nomad.app.dcat.mapping import Mapping

from tests.conftest import clear_elastic
from tests.app.test_app import BlueprintClient


@pytest.fixture(scope='session')
def api(session_client):
    return BlueprintClient(session_client, '/dcat')


@pytest.fixture(scope='module')
def example_entry(test_user, other_test_user):

    entry = EntryMetadata(
        calc_id='test-id',
        upload_id='upload-id',
        upload_time=datetime.now(),
        last_processing=datetime.now(),
        uploader=test_user,
        coauthors=[other_test_user],
        comment='this is a calculation comment',
        formula='H20',
        published=True)

    yield entry


def test_mapping(example_entry):
    mapping = Mapping()
    mapping.map_entry(example_entry)
    assert mapping.g is not None
    # print(mapping.g.serialize(format='ttl').decode('utf-8'))


def test_get_dataset(elastic_infra, api, example_entry):
    clear_elastic(elastic_infra)
    example_entry.a_elastic.index()
    calc_id = 'test-id'
    rv = api.get('/datasets/%s' % calc_id)
    assert rv.status_code == 200

    clear_elastic(elastic_infra)


def test_get_catalog(elastic_infra, api, example_entry):
    clear_elastic(elastic_infra)

    for i in range(1, 11):
        example_entry.calc_id = 'test-id-%d' % i
        example_entry.upload_time = datetime(2000, 1, 1)
        example_entry.last_processing = datetime(2020, 1, i)
        example_entry.a_elastic.index()

    infrastructure.elastic_client.indices.refresh(index=config.elastic.index_name)

    rv = api.get('/catalog/?after=test-id-3&modified_since=2020-01-07&format=nt')
    assert rv.status_code == 200

    clear_elastic(elastic_infra)
