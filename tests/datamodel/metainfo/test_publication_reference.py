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

import datetime
import json

import pytest
import requests

from tests.normalizing.conftest import run_processing


class MockResponse:
    def __init__(self, json_data: dict):
        self._json = json_data
        self.status_code = 200

    def json(self):
        return self._json

    @property
    def ok(self):
        return True


@pytest.fixture(scope='function')
def crossref_api_fixture(monkeypatch):
    with open('tests/data/datamodel/metainfo/external_api_mock_response.json') as fp:
        mock_responses = json.load(fp)

    def mock_requests_get(url, **kwargs):
        return MockResponse(mock_responses['crossref_api'])

    monkeypatch.setattr(requests, 'get', mock_requests_get)


def test_publication_reference(crossref_api_fixture):
    directory = 'tests/data/datamodel/metainfo'
    mainfile = 'test_publication_reference.archive.yaml'
    test_archive = run_processing(directory, mainfile)

    entry_data = test_archive.data

    # DOI should be prefixed with https://doi.org/
    assert entry_data.DOI_number == 'https://doi.org/10.1002/pol.1959.1203512832'

    # Authors should be populated
    assert entry_data.publication_authors is not None
    assert len(entry_data.publication_authors) == 2
    assert 'S. Banerjee' in entry_data.publication_authors
    assert 'M. S. Muthana' in entry_data.publication_authors

    # Journal should be populated
    assert entry_data.journal == 'Journal of Polymer Science'

    # Title should be populated
    assert (
        entry_data.publication_title == 'Studies on copolymerization of vinyl benzoate'
    )

    # Publication date 1959-02 should be stored as 1959-02-01 (day defaults to 1)
    assert isinstance(entry_data.publication_date, datetime.datetime)
    assert entry_data.publication_date.replace(tzinfo=None) == datetime.datetime(
        1959, 2, 1
    )

    # DOI should appear in archive metadata references
    assert test_archive.metadata.references is not None
    assert (
        'https://doi.org/10.1002/pol.1959.1203512832'
        in test_archive.metadata.references
    )
