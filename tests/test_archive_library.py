import pytest
import os

from nomad.archive import ArchiveFileDB
from nomad.archive_query import ArchiveQuery, ArchiveMetainfo
from tests.app.test_app import BlueprintClient


@pytest.fixture(scope='function')
def example_msgdb():
    def create_msgdb(payload):
        filename = 'archive_test.msg'
        msgdbo = ArchiveFileDB(filename, mode='w', entry_toc_depth=1)
        msgdbo.add_data(payload)
        msgdbo.close()
        msgdbo = ArchiveFileDB(filename, mode='r')
        return msgdbo

    filename = 'archive_test.msg'
    yield create_msgdb
    os.remove(filename)


class TestArchiveMetainfo:
    @pytest.fixture(scope='function')
    def data(self):
        return [
            {'calc_1': {'secA': {'propA': 1.0, 'propB': 'X'}}},
            {'calc_2': {'secA': {'propA': 2.0, 'propB': 'Y'}}}]

    def assert_metainfo(self, metainfo):
        for calc in metainfo.calcs:
            assert calc.secA({'propA': '*'}) is not None
            assert calc({'secA': {'propA': '*', 'propB': '*'}}) is not None

    def test_query_from_file(self, data, example_msgdb):
        _ = example_msgdb(data)
        metainfo = ArchiveMetainfo(archive_data='archive_test.msg', archive_schema={'secA': '*'})
        self.assert_metainfo(metainfo)

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
        q_params = {'Pagination': {'order': 1, 'per_page': 5}}
        q_schema = {'section_entry_info': '*'}
        q = ArchiveQuery(q_params, archive_data=q_schema, authentication=test_user_auth)
        q.query()
        for calc in q.metainfo:
            assert calc.section_entry_info.calc_id is not None

    def test_query_from_kwargs(self, api, published_wo_user_metadata, other_test_user_auth, monkeypatch):
        monkeypatch.setattr('nomad.archive_query.requests', api)
        q_schema = {'section_entry_info': '*'}
        q = ArchiveQuery(order=1, per_page=5, scroll=True, archive_data=q_schema, authentication=other_test_user_auth)
        q.query()
        for calc in q.metainfo:
            assert calc.section_entry_info.calc_id is not None
