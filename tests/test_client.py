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
from typing import List, Tuple
import pytest

from nomad.client import query_archive
from nomad.metainfo import MSection, SubSection
from nomad.datamodel import EntryArchive, User
from nomad.datamodel.metainfo.public import section_run

from tests.app.flask.conftest import client, session_client  # pylint: disable=unused-import
from tests.app.conftest import other_test_user_auth, test_user_auth  # pylint: disable=unused-import
from tests.app.flask.test_app import BlueprintClient
from tests.processing import test_data as test_processing


# TODO with the existing published_wo_user_metadata fixture there is only one entry
# that does not allow to properly test pagination and scrolling

@pytest.fixture(scope='function')
def api(client, monkeypatch):
    monkeypatch.setattr('nomad.config.client.url', '')
    api = BlueprintClient(client, '/api')
    monkeypatch.setattr('nomad.client.requests', api)
    return api


def assert_results(
        results: List[MSection],
        sub_section_defs: List[SubSection] = None,
        total=1):
    assert len(results) == total
    for result in results:
        assert result.m_def == EntryArchive.m_def
        if sub_section_defs:
            current = result
            for sub_section_def in sub_section_defs:
                for other_sub_section_def in current.m_def.all_sub_sections.values():
                    if other_sub_section_def != sub_section_def:
                        assert len(current.m_get_sub_sections(other_sub_section_def)) == 0

                sub_sections = current.m_get_sub_sections(sub_section_def)
                assert len(sub_sections) > 0
                current = sub_sections[0]


def test_query(api, published_wo_user_metadata):
    assert_results(query_archive())


def test_query_query(api, published_wo_user_metadata):
    assert_results(query_archive(query=dict(upload_id=[published_wo_user_metadata.upload_id])))


@pytest.mark.parametrize('q_schema,sub_sections', [
    ({'section_run': '*'}, [EntryArchive.section_run]),
    ({'section_run': {'section_system': '*'}}, [EntryArchive.section_run, section_run.section_system]),
    ({'section_run[0]': {'section_system': '*'}}, [EntryArchive.section_run, section_run.section_system])
])
def test_query_schema(api, published_wo_user_metadata, q_schema, sub_sections):
    assert_results(query_archive(required=q_schema), sub_section_defs=sub_sections)


def test_query_authentication(api, published, other_test_user_auth, test_user_auth, other_test_user):
    # The published test uploads uploader in calc and upload's user id do not match
    # due to testing the uploader change via publish metadata.

    assert_results(query_archive(authentication=other_test_user_auth), total=0)
    assert_results(query_archive(authentication=test_user_auth), total=1)


@pytest.fixture(scope='function')
def many_uploads(non_empty_uploaded: Tuple[str, str], test_user: User, proc_infra):
    _, upload_file = non_empty_uploaded
    for index in range(0, 4):
        upload = test_processing.run_processing(('test_upload_%d' % index, upload_file), test_user)
        upload.publish_upload()  # pylint: disable=no-member
        try:
            upload.block_until_complete(interval=.01)
        except Exception:
            pass


@pytest.fixture(scope='function', autouse=True)
def patch_multiprocessing_and_api(monkeypatch):
    class TestPool:
        ''' A fake multiprocessing pool, because multiprocessing does not work well in pytest. '''
        def __init__(self, n):
            pass

        def map(self, f, args):
            return [f(arg) for arg in args]

        def __enter__(self, *args, **kwargs):
            return self

        def __exit__(self, *args, **kwargs):
            pass

    monkeypatch.setattr('multiprocessing.Pool', TestPool)
    monkeypatch.setattr('nomad.client.get_json', lambda response: response.json)
    monkeypatch.setattr('nomad.client.get_length', lambda response: int(response.headers['Content-Length']))


def test_parallel_query(api, many_uploads, monkeypatch):
    result = query_archive(required=dict(section_run='*'), parallel=2)
    assert_results(result, total=4)
    assert result._statistics.nentries == 4
    assert result._statistics.loaded_nentries == 4
    assert result._statistics.last_response_nentries == 4
    assert result._statistics.napi_calls == 1

    result = query_archive(required=dict(section_run='*'), parallel=2, per_page=1)
    assert_results(result, total=4)
    assert result._statistics.nentries == 4
    assert result._statistics.loaded_nentries == 4
    assert result._statistics.last_response_nentries == 2
    assert result._statistics.napi_calls == 2
