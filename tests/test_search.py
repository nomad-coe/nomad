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

from typing import List, Iterable
from elasticsearch_dsl import Q
import pytest

from nomad import datamodel, search, processing, parsing, infrastructure, config
from nomad.search import Entry, SearchRequest


def test_init_mapping(elastic):
    pass


def test_index_skeleton_calc(elastic):
    entry_metadata = datamodel.EntryMetadata(
        domain='dft', upload_id='test_upload', calc_id='test_calc',
        mainfile='test/mainfile', files=['test/file1', 'test/file2'])

    create_entry(entry_metadata)


def test_index_normalized_calc(elastic, normalized: parsing.LocalBackend):
    entry_metadata = datamodel.EntryMetadata(
        domain='dft', upload_id='test upload id', calc_id='test id')
    entry_metadata.apply_domain_metadata(normalized)

    search_entry = create_entry(entry_metadata)
    entry = search.flat(search_entry.to_dict())

    assert 'calc_id' in entry
    assert 'atoms' in entry
    assert 'dft.code_name' in entry
    assert 'dft.optimade.elements_ratios' in entry


def test_index_normalized_calc_with_metadata(
        elastic, normalized: parsing.LocalBackend, internal_example_user_metadata: dict):
    entry_metadata = datamodel.EntryMetadata(
        domain='dft', upload_id='test upload id', calc_id='test id')
    entry_metadata.apply_domain_metadata(normalized)
    internal_example_user_metadata.pop('embargo_length')  # is for uploads only
    entry_metadata.apply_user_metadata(internal_example_user_metadata)

    entry = create_entry(entry_metadata)

    assert getattr(entry, 'with_embargo') == internal_example_user_metadata['with_embargo']
    assert getattr(entry, 'comment') == internal_example_user_metadata['comment']


def test_index_upload(elastic, processed: processing.Upload):
    pass


@pytest.fixture()
def example_search_data(elastic, normalized: parsing.LocalBackend):
    entry_metadata = datamodel.EntryMetadata(
        domain='dft', upload_id='test upload id', calc_id='test id')
    entry_metadata.apply_domain_metadata(normalized)
    create_entry(entry_metadata)
    refresh_index()

    return normalized


@pytest.fixture()
def example_ems_search_data(elastic, parsed_ems: parsing.LocalBackend):
    entry_metadata = datamodel.EntryMetadata(
        domain='ems', upload_id='test upload id', calc_id='test id')
    entry_metadata.apply_domain_metadata(parsed_ems)
    create_entry(entry_metadata)
    refresh_index()

    return parsed_ems


def test_search_entry(example_search_data):
    results = SearchRequest(domain='dft').execute()
    assert results['total'] > 0


def test_search_scan(elastic, example_search_data):
    results = list(SearchRequest(domain='dft').execute_scan())
    assert len(results) > 0


def test_search_paginated(elastic, example_search_data):
    results = SearchRequest(domain='dft').execute_paginated()
    assert results['total'] > 0
    assert len(results['results']) > 0
    pagination = results['pagination']
    assert pagination['total'] > 0
    assert 'page' in pagination
    assert 'per_page' in pagination


def test_search_scroll(elastic, example_search_data):
    request = SearchRequest(domain='dft')
    results = request.execute_scrolled()
    scroll_id = results['scroll']['scroll_id']
    assert results['scroll']['total'] == 1
    assert len(results['results']) == 1
    assert scroll_id is not None

    results = request.execute_scrolled(scroll_id=scroll_id)
    assert results['scroll']['total'] == 1
    assert len(results['results']) == 0
    assert 'scroll_id' not in results['scroll']


def test_domain(elastic, example_ems_search_data):
    assert len(list(SearchRequest(domain='ems').execute_scan())) > 0
    assert len(list(SearchRequest(domain='ems').domain().execute_scan())) > 0
    assert len(list(SearchRequest(domain='ems').domain('dft').execute_scan())) == 0
    assert len(list(SearchRequest(domain='dft').domain('dft').execute_scan())) == 0

    results = SearchRequest(domain='ems').statistic('ems.method', size=10).execute()
    statistics = results['statistics']
    assert 'ems.method' in statistics
    assert 'Bare eyes' in statistics['ems.method']

    results = SearchRequest(domain='ems').default_statistics().execute()
    statistics = results['statistics']
    assert 'ems.method' in statistics
    assert 'Bare eyes' in statistics['ems.method']


def assert_metrics(container, metrics_names):
    assert container['code_runs'] == 1
    for metric in metrics_names:
        assert metric in container


def test_search_statistics(elastic, example_search_data):
    assert 'authors' in search.metrics_names
    assert 'datasets' in search.metrics_names
    assert 'unique_entries' in search.metrics_names

    use_metrics = search.metrics_names

    request = SearchRequest(domain='dft').statistic(
        'dft.system', size=10, metrics_to_use=use_metrics).date_histogram()
    results = request.execute()

    statistics = results['statistics']
    assert 'results' not in results
    assert 'bulk' in statistics['dft.system']

    example_statistic = statistics['dft.system']['bulk']
    assert_metrics(example_statistic, use_metrics)
    assert_metrics(statistics['total']['all'], [])

    assert 'quantities' not in results


def test_search_totals(elastic, example_search_data):
    use_metrics = search.metrics_names

    request = SearchRequest(domain='dft').totals(metrics_to_use=use_metrics)
    results = request.execute()

    statistics = results['statistics']
    assert 'results' not in results
    assert len(statistics) == 1

    assert_metrics(statistics['total']['all'], [])

    assert 'quantities' not in results


def test_search_exclude(elastic, example_search_data):
    for item in SearchRequest().execute_paginated()['results']:
        assert 'atoms' in search.flat(item)

    for item in SearchRequest().exclude('atoms').execute_paginated()['results']:
        assert 'atoms' not in search.flat(item)


def test_search_include(elastic, example_search_data):
    for item in SearchRequest().execute_paginated()['results']:
        assert 'atoms' in search.flat(item)

    for item in SearchRequest().include('calc_id').execute_paginated()['results']:
        item = search.flat(item)
        assert 'atoms' not in item
        assert 'calc_id' in item


@pytest.mark.parametrize("order_by", [None, 'upload_id'])
def test_search_quantity(
        elastic, normalized: parsing.LocalBackend, test_user: datamodel.User,
        other_test_user: datamodel.User, order_by: str):

    entry_metadata = datamodel.EntryMetadata(
        domain='dft', upload_id='test upload id', calc_id='test id')
    entry_metadata.apply_domain_metadata(normalized)
    entry_metadata.uploader = test_user.user_id
    create_entry(entry_metadata)

    entry_metadata.calc_id = 'other test id'
    entry_metadata.uploader = other_test_user.user_id
    create_entry(entry_metadata)
    refresh_index()

    request = SearchRequest(domain='dft').quantity(
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


def create_entry(entry_metadata: datamodel.EntryMetadata):
    entry = search.create_entry(entry_metadata)
    entry.save()
    assert_entry(entry_metadata.calc_id)
    return entry


def assert_entry(calc_id):
    refresh_index()
    calc = Entry.get(calc_id)
    assert calc is not None

    search = Entry.search().query(Q('term', calc_id=calc_id))[0:10]
    assert search.count() == 1
    results = list(hit.to_dict() for hit in search)
    assert results[0]['calc_id'] == calc_id


def assert_search_upload(
        upload_entries: Iterable[datamodel.EntryMetadata],
        additional_keys: List[str] = [], **kwargs):
    keys = ['calc_id', 'upload_id', 'mainfile', 'calc_hash']
    refresh_index()
    search_results = Entry.search().query('match_all')[0:10]
    assert search_results.count() == len(list(upload_entries))
    if search_results.count() > 0:
        for hit in search_results:
            hit = search.flat(hit.to_dict())

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
            calc = Entry.from_entry_metadata(calc)
            yield calc.to_dict(include_meta=True)

    bulk(infrastructure.elastic_client, gen_data())
