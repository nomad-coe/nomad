import pytest

from nomad.archive_query import ArchiveQuery, ArchiveMetainfo
from tests.app.test_app import BlueprintClient


class TestArchiveMetainfo:
    @pytest.fixture(scope='function')
    def data(self):
        return [
            {'calc_1': {'secA': {'propA': 1.0, 'propB': 'X'}}},
            {'calc_2': {'secA': {'propA': 2.0, 'propB': 'Y'}}}]

    def assert_metainfo(self, metainfo):
        for calc in metainfo.calcs:
            assert isinstance(calc.secA.propA, float)
            assert calc.secA.m_to_dict() is not None

    def test_query_from_data(self, data):
        metainfo = ArchiveMetainfo(archive_data=data)
        self.assert_metainfo(metainfo)


class TestArchiveQuery:
    @pytest.fixture(scope='function')
    def api(self, client, monkeypatch):
        monkeypatch.setattr('nomad.config.client.url', '')
        return BlueprintClient(client, '/api')

    def test_query_from_json(self, api, published_wo_user_metadata, test_user_auth, monkeypatch):
        monkeypatch.setattr('nomad.archive_query.requests', api)
        q_params = {'pagination': {'order': 1, 'per_page': 5}}
        q_schema = {'section_entry_info': '*'}
        q = ArchiveQuery(q_params, query_schema=q_schema, authentication=test_user_auth)
        q.query()
        for calc in q.metainfo:
            assert calc.section_entry_info.calc_id is not None

    def test_query_from_kwargs(self, api, published_wo_user_metadata, other_test_user_auth, monkeypatch):
        monkeypatch.setattr('nomad.archive_query.requests', api)
        q_schema = {'section_entry_info': '*'}
        q = ArchiveQuery(
            scroll=dict(scroll=True), pagination=dict(per_page=5), query_schema=q_schema,
            authentication=other_test_user_auth)
        q.query()
        for calc in q.metainfo:
            assert calc.section_entry_info.calc_id is not None
