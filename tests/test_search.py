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
import pytest

from nomad import datamodel, search, processing, parsing, infrastructure, config
from nomad.search import Entry, SearchRequest


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


@pytest.fixture()
def example_search_data(elastic, normalized: parsing.LocalBackend):
    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test upload id', calc_id='test id')
    calc_with_metadata.apply_domain_metadata(normalized)
    create_entry(calc_with_metadata)
    refresh_index()

    return normalized


def test_search_entry(example_search_data):
    results = SearchRequest().execute()
    assert results['total'] > 0


def test_search_scan(elastic, example_search_data):
    results = list(SearchRequest().execute_scan())
    assert len(results) > 0


def test_search_paginated(elastic, example_search_data):
    results = SearchRequest().execute_paginated()
    assert results['total'] > 0
    assert len(results['results']) > 0
    pagination = results['pagination']
    assert pagination['total'] > 0
    assert 'page' in pagination
    assert 'per_page' in pagination


def test_search_scroll(elastic, example_search_data):
    request = SearchRequest()
    results = request.execute_scrolled()
    scroll_id = results['scroll']['scroll_id']
    assert results['scroll']['total'] == 1
    assert len(results['results']) == 1
    assert scroll_id is not None

    results = request.execute_scrolled(scroll_id=scroll_id)
    assert results['scroll']['total'] == 1
    assert len(results['results']) == 0
    assert 'scroll_id' not in results['scroll']


def assert_metrics(container, metrics_names):
    assert container['code_runs'] == 1
    for metric in metrics_names:
        assert metric in container


def test_search_statistics(elastic, example_search_data):
    assert 'authors' in search.metrics_names
    assert 'datasets' in search.metrics_names
    assert 'unique_entries' in search.metrics_names

    use_metrics = search.metrics_names

    request = SearchRequest().statistic('system', size=10, metrics_to_use=use_metrics).date_histogram()
    results = request.execute()

    statistics = results['statistics']
    assert 'results' not in results
    assert 'bulk' in statistics['system']

    example_statistic = statistics['system']['bulk']
    assert_metrics(example_statistic, use_metrics)
    assert_metrics(statistics['total']['all'], [])

    assert 'quantities' not in results


def test_search_totals(elastic, example_search_data):
    use_metrics = search.metrics_names

    request = SearchRequest().totals(metrics_to_use=use_metrics)
    results = request.execute()

    statistics = results['statistics']
    assert 'results' not in results
    assert len(statistics) == 1

    assert_metrics(statistics['total']['all'], [])

    assert 'quantities' not in results


def test_search_exclude(elastic, example_search_data):
    for item in SearchRequest().execute_paginated()['results']:
        assert 'atoms' in item

    for item in SearchRequest().exclude('atoms').execute_paginated()['results']:
        assert 'atoms' not in item


def test_search_include(elastic, example_search_data):
    for item in SearchRequest().execute_paginated()['results']:
        assert 'atoms' in item

    for item in SearchRequest().include('calc_id').execute_paginated()['results']:
        assert 'atoms' not in item
        assert 'calc_id' in item


@pytest.mark.parametrize("order_by", [None, 'upload_id'])
def test_search_quantity(
        elastic, normalized: parsing.LocalBackend, test_user: datamodel.User,
        other_test_user: datamodel.User, order_by: str):

    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test upload id', calc_id='test id')
    calc_with_metadata.apply_domain_metadata(normalized)
    calc_with_metadata.uploader = test_user.user_id
    create_entry(calc_with_metadata)

    calc_with_metadata.calc_id = 'other test id'
    calc_with_metadata.uploader = other_test_user.user_id
    create_entry(calc_with_metadata)
    refresh_index()

    request = SearchRequest().quantity(
        name='authors', size=1, examples=1, order_by=order_by)
    results = request.execute()
    assert len(results['quantities']['authors']['values'].keys()) == 1
    name = list(results['quantities']['authors']['values'].keys())[0]
    assert len(results['quantities']['authors']['values'][name]['examples']) == 1
    if order_by is None:
        assert results['quantities']['authors']['after'] == name
    else:
        assert results['quantities']['authors']['after'] == \
            results['quantities']['authors']['values'][name]['examples'][0][order_by]


def refresh_index():
    infrastructure.elastic_client.indices.refresh(index=config.elastic.index_name)


def create_entry(calc_with_metadata: datamodel.CalcWithMetadata):
    entry = search.Entry.from_calc_with_metadata(calc_with_metadata)
    entry.save()
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
    infrastructure.setup_mongo()
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
