import pytest
import os

from nomad.archive_library.filedb import ArchiveFileDB
from nomad.archive_library.metainfo import ArchiveMetainfo
from nomad.archive_library.query import ArchiveQuery
from tests.app.test_app import BlueprintClient


@pytest.fixture(scope='function')
def example_msgdb():
    def create_msgdb(payload):
        filename = 'archive_test.msg'
        msgdbo = ArchiveFileDB(filename, mode='w', max_lfragment=1)
        msgdbo.add_data(payload)
        msgdbo.close()
        msgdbo = ArchiveFileDB(filename, mode='r')
        return msgdbo

    filename = 'archive_test.msg'
    yield create_msgdb
    os.remove(filename)


class TestArchiveFileDB:
    def get_value(self, data, key):
        if key in data:
            return data[key]
        for v in data.values():
            return self.get_value(v, key)

    def get_keys(self, data, key=''):
        if data is None:
            return [key]
        keys = []
        for k, v in data.items():
            keys += self.get_keys(v, k)
        return keys

    @pytest.mark.parametrize('payload', [
        [
            'tests/data/proc/examples_archive/3Sqa0yIQnBrAautsn38YNhyZrOoE.json',
            'tests/data/proc/examples_archive/3Sqa0yIQnBrAautsn38YNhyZrOoE.log'],
        [
            {'secA': {'propA': 'X'}}]])
    def test_pack(self, example_msgdb, payload):
        fo = example_msgdb(payload)
        assert fo is not None

    @pytest.mark.parametrize('schema, dtype', [
        ({'secA': {'subsecA1[0]': {'propA1a': None}}}, {'propA1a': float}),
        ({'secB': {'propB1a': None}}, {'propB1a': list}),
        ({'secA': {'subsecA1[-1:]': {'propA1a': None}}}, {'propA1a': list})])
    def test_query(self, example_msgdb, schema, dtype):
        payload = [
            {'calc1': {
                'secA': {'subsecA1': [{'propA1a': 1.0}]}, 'secB': {'propB1a': ['a', 'b']}}},
            {'calc2': {
                'secA': {'subsecA1': [{'propA1a': 2.0}]}, 'secB': {'propB1a': ['c', 'd']}}}]
        msgdb = example_msgdb(payload)
        calc_ids = msgdb.ids.keys()
        calc_ids = [c for c in calc_ids if not os.path.dirname(c)]
        calc_ids = [c for c in calc_ids if not c.endswith('log') and c]
        qs = {calc_id: schema for calc_id in calc_ids}
        results = msgdb.query(qs)
        assert len(results) == len(calc_ids)
        for calc_id in results:
            for key in self.get_keys(schema):
                assert(isinstance(self.get_value(results[calc_id], key), dtype[key]))

    def test_error(self, example_msgdb):
        vals = [
            [{'atom_labels': ['X']}, {'atom_labels': ['X', 'X']}],
            [{'atom_labels': ['X']}], {'atom_labels': 'X'}]
        payload = [{'calc_%d' % i: {'section_run': {'section_system': vals[i]}}} for i in range(len(vals))]
        msgdb = example_msgdb(payload)
        # invalid key
        qs = {'calc_%d' % i: {'sction_run': {'section_system[:]': {'atom_labels': None}}} for i in range(len(vals))}
        results = msgdb.query(qs)
        assert results == qs
        # invalid index
        qs = {'calc_0': {'section_run': {'section_system[-3:-1]': None}}}
        results = msgdb.query(qs)
        assert results == qs
        # invalid calculation
        qs = {'calc_100': None}
        results = msgdb.query(qs)
        assert results == {}


class TestArchiveMetainfo:
    @pytest.fixture(scope='function')
    def data(self):
        return [
            {'calc_1': {'secA': {'propA': 1.0, 'propB': 'X'}}},
            {'calc_2': {'secA': {'propA': 2.0, 'propB': 'Y'}}}]

    def assert_metainfo(self, metainfo):
        for calc in metainfo.calcs:
            assert calc.secA({'propA': None}) is not None
            assert calc({'secA': {'propA': None, 'propB': None}}) is not None

    def test_query_from_file(self, data, example_msgdb):
        _ = example_msgdb(data)
        metainfo = ArchiveMetainfo(archive_data='archive_test.msg')
        self.assert_metainfo(metainfo)

    def test_query_from_data(self, data):
        metainfo = ArchiveMetainfo(archive_data=data)
        self.assert_metainfo(metainfo)


class TestArchiveQuery:
    @pytest.fixture(scope='function')
    def api(self, client, monkeypatch):
        monkeypatch.setattr('nomad.config.api_url', lambda *args, **kwargs: '')
        return BlueprintClient(client, '/api')

    @pytest.mark.parametrize('db', ['zip', 'msg'])
    def test_query_from_json(self, api, published_wo_user_metadata, other_test_user_auth, db, monkeypatch):
        monkeypatch.setattr('nomad.archive_library.query.requests', api)
        q_params = {'order': 1, 'per_page': 5, 'scroll': False, 'db': db}
        q_schema = {'section_entry_info': None}
        q = ArchiveQuery(q_params, q_schema)
        metainfo = q.query()
        for c in metainfo.calcs:
            assert c.section_entry_info({'calc_id': None}) is not None
