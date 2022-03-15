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
from datetime import datetime
from fastapi.testclient import TestClient

from nomad.app.dcat.main import app
from nomad.app.dcat.mapping import Mapping
from nomad.datamodel.results import Material, Results
from nomad.datamodel import Dataset
from nomad.utils.exampledata import ExampleData


@pytest.fixture(scope='session')
def api():
    return TestClient(app, base_url='http://testserver/')


def create_dataset(**kwargs):
    dataset = Dataset(
        dataset_create_time=datetime.now(), dataset_modified_time=datetime.now(),
        **kwargs)
    dataset.m_get_annotations('mongo').save()
    return dataset


@pytest.fixture(scope='module')
def data(test_user, other_test_user, elastic_infra, mongo_infra):
    example_attrs = dict(
        entry_id='test-id',
        upload_id='upload-id',
        last_processing_time=datetime.now(),
        entry_coauthors=[other_test_user],
        datasets=[
            create_dataset(
                dataset_id='dataset_1',
                user_id=test_user.user_id,
                dataset_name='test dataset 1',
                dataset_type='owned',
                doi='TEST/DOI'),
            create_dataset(
                dataset_id='dataset_2',
                user_id=test_user.user_id,
                dataset_name='test dataset 2',
                dataset_type='owned')
        ],
        comment='this is an entry comment')

    data = ExampleData(main_author=test_user)
    data.create_upload(
        upload_id='upload-id', upload_create_time=datetime(2000, 1, 1), published=True, embargo_length=0)
    archive = data.create_entry(**example_attrs)
    archive.m_create(Results).m_create(Material).chemical_formula_descriptive = 'H2O'

    for i in range(1, 11):
        example_attrs.update(
            entry_id='test-id-%d' % i,
            last_processing_time=datetime(2020, 1, i))
        data.create_entry(**example_attrs)

    data.save(with_files=False)

    return data


@pytest.fixture(scope='module')
def example_entry(data):
    return data.entries['test-id']


def test_mapping(example_entry):
    mapping = Mapping()
    mapping.map_entry(example_entry)
    assert mapping.g is not None


def test_get_dataset(api, example_entry):
    entry_id = 'test-id'
    rv = api.get('/datasets/%s' % entry_id)
    assert rv.status_code == 200


@pytest.mark.parametrize('after,modified_since', [
    (None, None),
    (None, '2020-01-07'),
    ('test-id-3', '2020-01-07')])
def test_get_catalog(api, data, after, modified_since):
    url = '/catalog/?format=turtle'
    if after:
        url += '&after=' + after
    if modified_since:
        url += '&modified_since=' + modified_since
    rv = api.get(url)
    assert rv.status_code == 200
