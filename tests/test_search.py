# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from elasticsearch_dsl import Q

from nomad import datamodel, search, processing, parsing, infrastructure, config, coe_repo
from nomad.search import Entry, aggregate_search, authors


def test_init_mapping(elastic):
    pass


def test_index_skeleton_calc(elastic):
    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test_upload', calc_id='test_calc')

    create_entry(calc_with_metadata)


def test_index_normalized_calc(elastic, normalized: parsing.LocalBackend):
    calc_with_metadata = normalized.to_calc_with_metadata()

    create_entry(calc_with_metadata)


def test_index_normalized_calc_with_metadata(
        elastic, normalized: parsing.LocalBackend, example_user_metadata: dict):

    calc_with_metadata = normalized.to_calc_with_metadata()
    calc_with_metadata.apply_user_metadata(example_user_metadata)

    create_entry(calc_with_metadata)


def test_index_upload(elastic, processed: processing.Upload):
    pass


def test_search(elastic, normalized: parsing.LocalBackend):
    calc_with_metadata = normalized.to_calc_with_metadata()
    create_entry(calc_with_metadata)
    refresh_index()

    total, hits, aggs = aggregate_search()
    assert total == 1
    assert hits[0]['calc_id'] == calc_with_metadata.calc_id
    assert 'Bulk' in aggs['system']
    assert aggs['system']['Bulk'] == 1


def test_authors(elastic, normalized: parsing.LocalBackend, test_user: coe_repo.User, other_test_user: coe_repo.User):
    calc_with_metadata = normalized.to_calc_with_metadata()
    calc_with_metadata.uploader = test_user.to_popo()
    create_entry(calc_with_metadata)
    calc_with_metadata.calc_id = 'other test calc'
    calc_with_metadata.uploader = other_test_user.to_popo()
    create_entry(calc_with_metadata)
    refresh_index()

    results, after = authors(per_page=1)
    assert len(results) == 1
    name = list(results.keys())[0]
    assert after == name


def refresh_index():
    infrastructure.elastic_client.indices.refresh(index=config.elastic.index_name)


def create_entry(calc_with_metadata: datamodel.CalcWithMetadata):
    search.Entry.from_calc_with_metadata(calc_with_metadata).save()
    assert_entry(calc_with_metadata.calc_id)


def assert_entry(calc_id):
    refresh_index()
    calc = Entry.get(calc_id)
    assert calc is not None

    search = Entry.search().query(Q('term', calc_id=calc_id))[0:10]
    assert search.count() == 1
    results = list(hit.to_dict() for hit in search)
    assert results[0]['calc_id'] == calc_id


def assert_search_upload(upload_id, published: bool = False):
    refresh_index()
    search = Entry.search().query('match_all')[0:10]
    if search.count() > 0:
        for hit in search:
            hit = hit.to_dict()
            if published:
                assert int(hit.get('pid')) > 0
                assert hit.get('published')

            for coauthor in hit.get('coauthors', []):
                assert coauthor.get('name', None) is not None


if __name__ == '__main__':
    from test_datamodel import generate_calc
    from elasticsearch.helpers import bulk
    import sys
    print('Generate index with random example calculation data. First arg is number of items')
    infrastructure.setup_elastic()
    n = 100
    if len(sys.argv) > 1:
        n = int(sys.argv[1])

    def gen_data():
        for pid in range(0, n):
            calc = generate_calc(pid)
            calc = Entry.from_calc_with_metadata(calc)
            yield calc.to_dict(include_meta=True)

    bulk(infrastructure.elastic_client, gen_data())
