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

from typing import List
from elasticsearch_dsl import Q

from nomad import datamodel, search, processing, parsing, infrastructure, config, coe_repo
from nomad.search import Entry, metrics_search, quantity_search, scroll_search, entry_search


def test_init_mapping(elastic):
    pass


def test_index_skeleton_calc(elastic):
    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test_upload', calc_id='test_calc')

    create_entry(calc_with_metadata)


def test_index_normalized_calc(elastic, normalized: parsing.LocalBackend):
    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test upload id', calc_id='test id')
    calc_with_metadata.apply_domain_metadata(normalized)

    create_entry(calc_with_metadata)


def test_index_normalized_calc_with_metadata(
        elastic, normalized: parsing.LocalBackend, example_user_metadata: dict):

    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test upload id', calc_id='test id')
    calc_with_metadata.apply_domain_metadata(normalized)
    calc_with_metadata.apply_user_metadata(example_user_metadata)

    create_entry(calc_with_metadata)


def test_index_upload(elastic, processed: processing.Upload):
    pass


def test_entry_search(elastic, normalized: parsing.LocalBackend):
    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test upload id', calc_id='test id')
    calc_with_metadata.apply_domain_metadata(normalized)
    create_entry(calc_with_metadata)
    refresh_index()

    results = entry_search()
    assert len(results['results']) > 0


def test_metrics_search(elastic, normalized: parsing.LocalBackend):
    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test upload id', calc_id='test id')
    calc_with_metadata.apply_domain_metadata(normalized)
    create_entry(calc_with_metadata)
    refresh_index()

    assert 'users' in search.metrics_names
    assert 'datasets' in search.metrics_names
    assert 'unique_code_runs' in search.metrics_names

    use_metrics = search.metrics_names

    results = metrics_search(metrics_to_use=use_metrics, with_entries=True, with_date_histogram=True)
    quantities = results['quantities']
    hits = results['results']
    assert results['pagination']['total'] == 1
    assert hits[0]['calc_id'] == calc_with_metadata.calc_id
    assert 'bulk' in quantities['system']

    example_quantity = quantities['system']['bulk']

    def assert_metrics(container, metrics_names):
        assert container['code_runs'] == 1
        for metric in metrics_names:
            assert metric in container

    assert_metrics(example_quantity, use_metrics)
    assert_metrics(quantities['total']['all'], use_metrics)

    assert 'quantities' not in hits[0]


def test_scroll_search(elastic, normalized: parsing.LocalBackend):
    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test upload id', calc_id='test id')
    calc_with_metadata.apply_domain_metadata(normalized)
    create_entry(calc_with_metadata)
    refresh_index()

    results = scroll_search()
    scroll_id = results['scroll']['scroll_id']
    assert results['scroll']['total'] == 1
    assert len(results['results']) == 1
    assert scroll_id is not None

    results = scroll_search(scroll_id=scroll_id)
    assert results['scroll']['total'] == 1
    assert len(results['results']) == 0
    assert 'scroll_id' not in results['scroll']


def test_quantity_search(elastic, normalized: parsing.LocalBackend, test_user: coe_repo.User, other_test_user: coe_repo.User):
    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test upload id', calc_id='test id')
    calc_with_metadata.apply_domain_metadata(normalized)
    calc_with_metadata.uploader = test_user.to_popo()
    create_entry(calc_with_metadata)
    calc_with_metadata.calc_id = 'other test id'
    calc_with_metadata.uploader = other_test_user.to_popo()
    create_entry(calc_with_metadata)
    refresh_index()

    results = quantity_search(quantities=dict(authors=None), size=1, with_entries=False)
    assert len(results['quantities']['authors']['values'].keys()) == 1
    name = list(results['quantities']['authors']['values'].keys())[0]
    assert results['quantities']['authors']['after'] == name


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


def assert_search_upload(upload: datamodel.UploadWithMetadata, additional_keys: List[str] = [], **kwargs):
    keys = ['calc_id', 'upload_id', 'mainfile', 'calc_hash']
    refresh_index()
    search = Entry.search().query('match_all')[0:10]
    assert search.count() == len(list(upload.calcs))
    if search.count() > 0:
        for hit in search:
            hit = hit.to_dict()
            for key, value in kwargs.items():
                assert hit.get(key, None) == value

            if 'pid' in hit:
                assert int(hit.get('pid')) > 0

            for key in keys:
                assert key in hit

            for key in additional_keys:
                assert key in hit
                assert hit[key] != config.services.unavailable_value

            for coauthor in hit.get('coauthors', []):
                assert coauthor.get('name', None) is not None


if __name__ == '__main__':
    from test_datamodel import generate_calc
    from elasticsearch.helpers import bulk
    import sys
    print('Generate index with random example calculation data. First arg is number of items')
    infrastructure.setup_logging()
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
