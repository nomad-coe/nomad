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

from httpx import AsyncClient
import pytest

from nomad.app.main import app
from nomad.client.archive import ArchiveQuery
from nomad.datamodel import EntryArchive, User
from nomad.datamodel.metainfo import runschema, SCHEMA_IMPORT_ERROR
from nomad.metainfo import MSection, SubSection
from tests.fixtures.users import test_users
from tests.processing import test_data as test_processing


# TODO: more tests


def assert_results(
    results: List[MSection], sub_section_defs: List[SubSection] = None, total=1
):
    assert len(results) == total
    for result in results:
        assert result.m_def == EntryArchive.m_def
        if sub_section_defs:
            current = result
            for sub_section_def in sub_section_defs:
                for other_sub_section_def in current.m_def.all_sub_sections.values():
                    if other_sub_section_def != sub_section_def:
                        assert (
                            len(current.m_get_sub_sections(other_sub_section_def)) == 0
                        )

                sub_sections = current.m_get_sub_sections(sub_section_def)
                assert len(sub_sections) > 0
                current = sub_sections[0]


@pytest.fixture(scope='function')
def many_uploads(non_empty_uploaded: Tuple[str, str], test_user: User, proc_infra):
    _, upload_file = non_empty_uploaded
    for index in range(0, 4):
        upload = test_processing.run_processing(
            ('test_upload_%d' % index, upload_file), test_user
        )
        upload.publish_upload()  # pylint: disable=no-member
        try:
            upload.block_until_complete(interval=0.01)
        except Exception:
            pass


@pytest.fixture(scope='session')
def async_api_v1(monkeysession):
    """
    This fixture provides an HTTP client with AsyncClient that accesses
    the fast api. The patch will redirect all requests to the fast api under test.
    """
    test_client = AsyncClient(app=app)

    monkeysession.setattr(
        'nomad.client.archive.ArchiveQuery._fetch_url',
        'http://testserver/api/v1/entries/query',
    )
    monkeysession.setattr(
        'nomad.client.archive.ArchiveQuery._download_url',
        'http://testserver/api/v1/entries/archive/query',
    )

    monkeysession.setattr('httpx.AsyncClient.get', getattr(test_client, 'get'))
    monkeysession.setattr('httpx.AsyncClient.put', getattr(test_client, 'put'))
    monkeysession.setattr('httpx.AsyncClient.post', getattr(test_client, 'post'))
    monkeysession.setattr('httpx.AsyncClient.delete', getattr(test_client, 'delete'))

    def mocked_auth_headers(self) -> dict:
        for user in test_users.values():
            if user['username'] == self.user or user['email'] == self.user:
                return dict(Authorization=f'Bearer {user["user_id"]}')
        return {}

    monkeysession.setattr('nomad.client.api.Auth.headers', mocked_auth_headers)

    return test_client


def test_async_query_basic(async_api_v1, published_wo_user_metadata):
    async_query = ArchiveQuery()

    assert_results(async_query.download())

    async_query = ArchiveQuery(
        query=dict(upload_id=[published_wo_user_metadata.upload_id])
    )

    assert_results(async_query.download())


@pytest.mark.skipif(runschema is None, reason=SCHEMA_IMPORT_ERROR)
@pytest.mark.parametrize(
    'q_required,sub_sections',
    [
        ({'run': '*'}, [EntryArchive.run]),
        ({'run': {'system': '*'}}, [EntryArchive.run, runschema.run.Run.system]),
        ({'run[0]': {'system': '*'}}, [EntryArchive.run, runschema.run.Run.system]),
    ],
)
def test_async_query_required(
    async_api_v1, published_wo_user_metadata, q_required, sub_sections
):
    async_query = ArchiveQuery(required=q_required)

    assert_results(async_query.download(), sub_section_defs=sub_sections)


def test_async_query_auth(async_api_v1, published, other_test_user, test_user):
    async_query = ArchiveQuery(username=other_test_user.username, password='password')

    assert_results(async_query.download(), total=0)

    async_query = ArchiveQuery(username=test_user.username, password='password')

    assert_results(async_query.download(), total=1)


def test_async_query_parallel(async_api_v1, many_uploads, monkeypatch):
    async_query = ArchiveQuery(required=dict(run='*'))

    assert_results(async_query.download(), total=4)
    assert_results(async_query.download(), total=0)

    async_query = ArchiveQuery(required=dict(run='*'), page_size=1)

    assert_results(async_query.download(), total=4)
