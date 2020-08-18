from typing import List
import pytest

from nomad.client import query_archive
from nomad.metainfo import MSection, SubSection
from nomad.datamodel import EntryArchive
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
    # TODO this is a workarround and should be changed in conftest, especially since we do not need this
    # feature anymore.
    published.user_id = other_test_user.user_id
    published.save()

    assert_results(query_archive(authentication=other_test_user_auth), total=1)
    assert_results(query_archive(authentication=test_user_auth), total=0)
