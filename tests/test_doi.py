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
from nomad import config
from nomad.doi import DOI, DOIException
import pytest
from unittest.mock import MagicMock


def test_create(mongo, test_user, no_warn):
    doi = DOI.create('the_title', test_user)

    assert DOI.objects(doi=doi.doi).first() is not None
    assert doi.metadata_xml is not None


def test_create_doi_counter(mongo, test_user, no_warn):
    DOI.create('the_title', test_user)
    doi = DOI.create('the_title', test_user)
    assert doi.doi.endswith('-2')


def test_create_draft_doi(mongo, test_user, no_warn):
    if config.datacite.enabled:
        doi = DOI.create('the_title', test_user)
        doi.create_draft()
        doi.delete()

        assert DOI.objects(doi=doi.doi).first() is None


@pytest.mark.parametrize('status_code,response_ok,is_findable,text', [
    pytest.param(200, True, False, 'Success', id='pass-with-200-only-drafted'),
    pytest.param(200, True, True, 'Success', id='pass-with-200-with-findable'),
    pytest.param(400, False, False, 'Bad Request', id='fail-with-400'),
])
def test_datacite_requests(mongo, monkeypatch, test_user, status_code, response_ok, is_findable, text):
    if config.datacite.enabled:
        def mock_datacite_request(*args, **kwargs):
            mock_response = MagicMock()
            mock_response.status_code = status_code
            mock_response.ok = response_ok
            mock_response.text = text
            return mock_response

        doi = DOI.create('the_title', test_user)

        monkeypatch.setattr(doi.create_draft.__globals__['requests'], 'post', mock_datacite_request)
        monkeypatch.setattr(doi.make_findable.__globals__['requests'], 'put', mock_datacite_request)
        monkeypatch.setattr(doi.delete.__globals__['requests'], 'delete', mock_datacite_request)

        if response_ok:
            doi.create_draft()
            assert DOI.objects(doi=doi.doi).first() is not None
            assert DOI.objects(doi=doi.doi).first().state == 'draft'

            if is_findable:
                doi.make_findable()
                assert DOI.objects(doi=doi.doi).first().state == 'findable'
            else:
                doi.delete()
                assert DOI.objects(doi=doi.doi).first() is None

        elif not response_ok and config.datacite.enabled:
            with pytest.raises(DOIException):
                doi.create_draft()
