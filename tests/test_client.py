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

from nomad.client import query_archive, Auth
from nomad.metainfo import MSection, SubSection
from nomad.datamodel import EntryArchive, User
from nomad.datamodel.metainfo.simulation.run import Run

from tests.app.conftest import other_test_user_auth, test_user_auth  # pylint: disable=unused-import
from tests.processing import test_data as test_processing

# TODO most nomad.client functionality is only tested indirectly via its use in nomad.cli


def test_requests_auth(api_v1):
    rv = api_v1.get('users/me', auth=Auth(from_api=True))
    assert rv.status_code == 200


# TODO with the existing published_wo_user_metadata fixture there is only one entry
# that does not allow to properly test pagination and scrolling

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


def test_query(api_v1, published_wo_user_metadata):
    assert_results(query_archive())


def test_query_query(api_v1, published_wo_user_metadata):
    assert_results(query_archive(query=dict(upload_id=[published_wo_user_metadata.upload_id])))


@pytest.mark.parametrize('q_schema,sub_sections', [
    ({'run': '*'}, [EntryArchive.run]),
    ({'run': {'system': '*'}}, [EntryArchive.run, Run.system]),
    ({'run[0]': {'system': '*'}}, [EntryArchive.run, Run.system])
])
def test_query_required(api_v1, published_wo_user_metadata, q_schema, sub_sections):
    assert_results(query_archive(required=q_schema), sub_section_defs=sub_sections)


def test_query_authentication(api_v1, published, other_test_user, test_user):
    # The published test uploads uploader in entry and upload's user id do not match
    # due to testing the uploader change via publish metadata.

    assert_results(query_archive(authentication=Auth(other_test_user.username, 'password', from_api=True)), total=0)
    assert_results(query_archive(authentication=Auth(test_user.username, 'password', from_api=True)), total=1)


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


def test_parallel_query(api_v1, many_uploads, monkeypatch):
    result = query_archive(required=dict(run='*'), parallel=2)
    assert_results(result, total=4)
    assert result._statistics.nentries == 4
    assert result._statistics.loaded_nentries == 4
    assert result._statistics.last_response_nentries == 4
    assert result._statistics.napi_calls == 1

    result = query_archive(required=dict(run='*'), parallel=2, per_page=1)
    assert_results(result, total=4)
    assert result._statistics.nentries == 4
    assert result._statistics.loaded_nentries == 4
    assert result._statistics.last_response_nentries == 2
    assert result._statistics.napi_calls == 2
