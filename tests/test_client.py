from typing import List
import pytest

from nomad.client import query_archive
from nomad.metainfo import MSection, SubSection
from nomad.datamodel.metainfo.public import section_run

from tests.app.test_app import BlueprintClient


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
        sub_sections: List[SubSection] = None,
        total=1):
    assert len(results) == total
    for result in results:
        if sub_sections:
            for sub_section in result.m_def.all_sub_sections.values():
                if sub_section in sub_sections:
                    assert len(result.m_get_sub_sections(sub_section)) > 0
                else:
                    assert len(result.m_get_sub_sections(sub_section)) == 0


def test_query(api, published_wo_user_metadata):
    assert_results(query_archive())


def test_query_query(api, published_wo_user_metadata):
    assert_results(query_archive(query=dict(upload_id=[published_wo_user_metadata.upload_id])))


def test_query_schema(api, published_wo_user_metadata):
    q_schema = {'section_run': {'section_system': '*'}}
    assert_results(
        query_archive(query_schema=q_schema),
        sub_sections=[section_run.section_system])


def test_query_scroll(api, published_wo_user_metadata):
    assert_results(query_archive(scroll=True))


def test_query_authentication(api, published, other_test_user_auth, test_user_auth):
    assert_results(query_archive(authentication=other_test_user_auth))
    assert_results(query_archive(authentication=test_user_auth), total=0)
