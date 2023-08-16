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

from typing import Iterable
import pytest
import json

from tests.normalizing.conftest import run_processing
from nomad.datamodel.metainfo import basesections


class MockResponse:

    def __init__(self, json: dict):
        self._json = json
        self.status_code = 200

    def json(self):
        return self._json

    @property
    def ok(self):
        return True


@pytest.fixture(scope='function')
def external_api_fixture(monkeypatch):

    with open('tests/data/datamodel/metainfo/external_api_mock_response.json', 'r') as fp:
        mock_responses = json.load(fp)

    def pub_chem_api_get_properties(cid: int, properties: Iterable[str]):
        return MockResponse(mock_responses['pub_chem_api_get_properties'][str(cid)])

    def pub_chem_api_get_synonyms(cid: int):
        return MockResponse(mock_responses['pub_chem_api_get_synonyms'][str(cid)])

    def pub_chem_api_search(path: str, search: str):
        return MockResponse(mock_responses['pub_chem_api_search'])

    def cas_api_search(search: str):
        return MockResponse(mock_responses['cas_api_search'])

    def cas_api_details(cas_rn: str):
        return MockResponse(mock_responses['cas_api_details'])

    monkeypatch.setattr(basesections, 'pub_chem_api_get_properties', pub_chem_api_get_properties)
    monkeypatch.setattr(basesections, 'pub_chem_api_get_synonyms', pub_chem_api_get_synonyms)
    monkeypatch.setattr(basesections, 'pub_chem_api_search', pub_chem_api_search)
    monkeypatch.setattr(basesections, 'cas_api_search', cas_api_search)
    monkeypatch.setattr(basesections, 'cas_api_details', cas_api_details)


@pytest.mark.parametrize('mainfile', [
    pytest.param('test_substance.archive.yaml', id='Substance'),
    pytest.param('test_cas_substance.archive.yaml', id='CASSubstance'),
    pytest.param('test_pub_chem_substance.archive.yaml', id='PubChemSubstance'),
    pytest.param('test_not_overwrite_substance.archive.yaml', id='not-overwrite'),
])
def test_substance(mainfile, external_api_fixture):
    directory = 'tests/data/datamodel/metainfo'
    test_archive = run_processing(directory, mainfile)

    # Check that the material results section was populated by the normalizer
    assert "I" in test_archive.results.material.elements
    assert "Pb" in test_archive.results.material.elements
    assert len(test_archive.results.material.elemental_composition) == 2
    for composition in test_archive.results.material.elemental_composition:
        if composition.element == 'I':
            assert composition.atomic_fraction == 2 / 3
        elif composition.element == 'Pb':
            assert composition.atomic_fraction == 1 / 3
        else:
            raise ValueError(
                f'Unknown element "{composition.element}" in'
                ' results.material.elemental_composition'
            )
