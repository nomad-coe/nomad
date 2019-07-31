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

"""
This module represents calculations in elastic search.
"""

from typing import Iterable, Dict, Tuple, List, Any
from elasticsearch_dsl import Document, InnerDoc, Keyword, Text, Date, \
    Object, Boolean, Search, Q, A, analyzer, tokenizer
from elasticsearch_dsl.document import IndexMeta
import elasticsearch.helpers
from elasticsearch.exceptions import NotFoundError
from datetime import datetime

from nomad import config, datamodel, infrastructure, datamodel, coe_repo, utils

path_analyzer = analyzer(
    'path_analyzer',
    tokenizer=tokenizer('path_tokenizer', 'pattern', pattern='/'))


class AlreadyExists(Exception): pass


class ElasticSearchError(Exception): pass


class ScrollIdNotFound(Exception): pass


class User(InnerDoc):

    @classmethod
    def from_user_popo(cls, user):
        self = cls(user_id=user.id)

        if 'first_name' not in user:
            user = coe_repo.User.from_user_id(user.id).to_popo()

        name = '%s, %s' % (user['last_name'], user['first_name'])
        self.name = name

        return self

    user_id = Keyword()
    name = Text(fields={'keyword': Keyword()})


class Dataset(InnerDoc):

    @classmethod
    def from_dataset_popo(cls, dataset):
        return cls(
            id=dataset.id,
            doi=dataset.doi['value'] if dataset.doi is not None else None,
            name=dataset.name)

    id = Keyword()
    doi = Keyword()
    name = Keyword()


class WithDomain(IndexMeta):
    """ Override elasticsearch_dsl metaclass to sneak in domain specific mappings """
    def __new__(cls, name, bases, attrs):
        for quantity in datamodel.Domain.instance.quantities.values():
            attrs[quantity.name] = quantity.elastic_mapping
        return super(WithDomain, cls).__new__(cls, name, bases, attrs)


class Entry(Document, metaclass=WithDomain):

    class Index:
        name = config.elastic.index_name

    upload_id = Keyword()
    upload_time = Date()
    calc_id = Keyword()
    calc_hash = Keyword()
    pid = Keyword()
    mainfile = Keyword()
    files = Text(multi=True, analyzer=path_analyzer, fields={'keyword': Keyword()})
    uploader = Object(User)

    with_embargo = Boolean()
    published = Boolean()

    processed = Boolean()
    last_processing = Date()
    nomad_version = Keyword()
    nomad_commit = Keyword()

    authors = Object(User, multi=True)
    owners = Object(User, multi=True)
    comment = Text()
    references = Keyword()
    datasets = Object(Dataset)

    @classmethod
    def from_calc_with_metadata(cls, source: datamodel.CalcWithMetadata) -> 'Entry':
        entry = Entry(meta=dict(id=source.calc_id))
        entry.update(source)
        return entry

    def update(self, source: datamodel.CalcWithMetadata) -> None:
        self.upload_id = source.upload_id
        self.upload_time = source.upload_time
        self.calc_id = source.calc_id
        self.calc_hash = source.calc_hash
        self.pid = None if source.pid is None else str(source.pid)

        self.processed = source.processed
        self.last_processing = source.last_processing
        self.nomad_version = source.nomad_version
        self.nomad_commit = source.nomad_commit

        self.mainfile = source.mainfile
        if source.files is None:
            self.files = [self.mainfile]
        elif self.mainfile not in source.files:
            self.files = [self.mainfile] + source.files
        else:
            self.files = source.files

        self.uploader = User.from_user_popo(source.uploader) if source.uploader is not None else None

        self.with_embargo = source.with_embargo
        self.published = source.published
        self.authors = [User.from_user_popo(user) for user in source.coauthors]
        self.owners = [User.from_user_popo(user) for user in source.shared_with]
        if self.uploader is not None:
            if self.uploader not in self.authors:
                self.authors.append(self.uploader)
            if self.uploader not in self.owners:
                self.owners.append(self.uploader)
        self.comment = source.comment
        self.references = [ref.value for ref in source.references]
        self.datasets = [Dataset.from_dataset_popo(ds) for ds in source.datasets]

        for quantity in datamodel.Domain.instance.quantities.values():
            setattr(
                self, quantity.name,
                quantity.elastic_value(getattr(source, quantity.metadata_field)))


def delete_upload(upload_id):
    """ Delete all entries with given ``upload_id`` from the index. """
    index = Entry._default_index()
    Search(index=index).query('match', upload_id=upload_id).delete()


def publish(calcs: Iterable[datamodel.CalcWithMetadata]) -> None:
    """ Update all given calcs with their metadata and set ``publish = True``. """
    def elastic_updates():
        for calc in calcs:
            entry = Entry.from_calc_with_metadata(calc)
            entry.published = True
            entry = entry.to_dict(include_meta=True)
            source = entry.pop('_source')
            entry['doc'] = source
            entry['_op_type'] = 'update'
            yield entry

    elasticsearch.helpers.bulk(infrastructure.elastic_client, elastic_updates())
    refresh()


def index_all(calcs: Iterable[datamodel.CalcWithMetadata]) -> None:
    """ Adds all given calcs with their metadata to the index. """
    def elastic_updates():
        for calc in calcs:
            entry = Entry.from_calc_with_metadata(calc)
            entry = entry.to_dict(include_meta=True)
            entry['_op_type'] = 'index'
            yield entry

    elasticsearch.helpers.bulk(infrastructure.elastic_client, elastic_updates())
    refresh()


def refresh():
    infrastructure.elastic_client.indices.refresh(config.elastic.index_name)


aggregations = datamodel.Domain.instance.aggregations
""" The available aggregations in :func:`aggregate_search` and their maximum aggregation size """

search_quantities = datamodel.Domain.instance.search_quantities
"""The available search quantities """

metrics = {
    'datasets': ('cardinality', 'datasets.id'),
    'unique_code_runs': ('cardinality', 'calc_hash')
}
"""
The available search metrics. Metrics are integer values given for each entry that can
be used in aggregations, e.g. the sum of all total energy calculations or cardinality of
all unique geometries.
"""

metrics.update(**datamodel.Domain.instance.metrics)

metrics_names = list(metric for metric in metrics.keys())


order_default_quantity = None
for quantity in datamodel.Domain.instance.quantities.values():
    if quantity.order_default:
        order_default_quantity = quantity.name


def _construct_search(
        q: Q = None, time_range: Tuple[datetime, datetime] = None,
        search_parameters: Dict[str, Any] = {}, **kwargs) -> Search:

    search = Search(index=config.elastic.index_name)

    if q is not None:
        search = search.query(q)

    if time_range is not None:
        search = search.query('range', upload_time=dict(gte=time_range[0], lte=time_range[1]))

    for key, value in search_parameters.items():
        quantity = search_quantities.get(key, None)
        if quantity is None:
            if key in ['page', 'per_page', 'order', 'order_by']:
                continue
            else:
                raise KeyError('Unknown quantity %s' % key)

        value = quantity.elastic_value(value)

        if isinstance(value, list):
            values = value
        else:
            values = [value]

        for item in values:
            search = search.query(Q(quantity.elastic_search_type, **{quantity.elastic_field: item}))

    search = search.source(exclude=['quantities'])
    return search


def _execute_paginated_search(
        search: Search,
        page: int = 1, per_page: int = 10,
        order_by: str = order_default_quantity, order: int = -1,
        **kwargs) -> Tuple[Any, Dict[str, Any]]:

    if order_by not in search_quantities:
        raise KeyError('Unknown order quantity %s' % order_by)

    order_by_quantity = search_quantities[order_by]

    if order == 1:
        search = search.sort(order_by_quantity.elastic_field)
    else:
        search = search.sort('-%s' % order_by_quantity.elastic_field)
    paginated_search = search[(page - 1) * per_page: page * per_page]
    response = paginated_search.execute()  # pylint: disable=E1101

    total_results = response.hits.total
    search_results = [hit.to_dict() for hit in response.hits]

    return response, {
        'pagination': {
            'page': page,
            'per_page': per_page,
            'total': total_results
        },
        'results': search_results
    }


def scroll_search(
        scroll_id: str = None, size: int = 1000, scroll: str = u'5m',
        q: Q = None, search_parameters: Dict[str, Any] = {}) -> Dict[str, Any]:
    """
    Alternative search based on ES scroll API. Can be used similar to
    :func:`aggregate_search`, but pagination is replaced with scrolling, no ordering,
    no property, and no metrics information is available.

    he search is limited to parameters :param:`q` and :param:`search_parameters`,
    which work exactly as in :func:`entry_search`.

    Scrolling is done by calling this function again and again with the same ``scroll_id``.
    Each time, this function will return the next batch of search results. If the
    ``scroll_id`` is not available anymore, a new ``scroll_id`` is assigned and scrolling
    starts from the beginning again.

    Arguments:
        scroll_id: The scroll id to receive the next batch from. None will create a new
            scroll.
        size: The batch size in number of hits.
        scroll: The time the scroll should be kept alive (i.e. the time between requests
            to this method) in ES time units. Default is 5 minutes.

    Returns:
        A dict with keys 'scroll' and 'results'. The key 'scroll' holds a dict with
        'total', 'scroll_id', 'size'.
    """
    es = infrastructure.elastic_client

    if scroll_id is None:
        # initiate scroll
        search = _construct_search(q, search_parameters=search_parameters)
        resp = es.search(body=search.to_dict(), scroll=scroll, size=size, index=config.elastic.index_name)  # pylint: disable=E1123

        scroll_id = resp.get('_scroll_id')
        if scroll_id is None:
            # no results for search query
            return dict(scroll=dict(total=0, size=size), results=[])

    else:
        try:
            resp = es.scroll(scroll_id, scroll=scroll)  # pylint: disable=E1123
        except NotFoundError:
            raise ScrollIdNotFound()

    total = resp['hits']['total']
    results = list(hit['_source'] for hit in resp['hits']['hits'])

    # since we are using the low level api here, we should check errors
    if resp["_shards"]["successful"] < resp["_shards"]["total"]:
        utils.get_logger(__name__).error('es operation was unsuccessful on at least one shard')
        raise ElasticSearchError('es operation was unsuccessful on at least one shard')

    if len(results) == 0:
        es.clear_scroll(body={'scroll_id': [scroll_id]}, ignore=(404, ))  # pylint: disable=E1123
        scroll_id = None

    scroll_info = dict(total=total, size=size)
    if scroll_id is not None:
        scroll_info.update(scroll_id=scroll_id)

    return dict(scroll=scroll_info, results=results)


def entry_search(
        q: Q = None,
        page: int = 1, per_page: int = 10,
        order_by: str = order_default_quantity, order: int = -1,
        time_range: Tuple[datetime, datetime] = None,
        search_parameters: Dict[str, Any] = {}) -> Dict[str, Any]:
    """
    Performs a search and returns a paginated list of search results.

    The search is determimed by the given elasticsearch_dsl query param:`q`,
    param:`time_range` and additional :param:`search_parameters`.
    The search_parameters have to match general or domain specific metadata quantities.
    See module:`datamodel`.

    The search results are paginated. Pagination is controlled by the pagination parameters
    param:`page` and param:`per_page`. The results are ordered.

    Arguments:
        page: The page to return starting with page 1
        per_page: Results per page
        q: An *elasticsearch_dsl* query used to further filter the results (via ``and``)
        time_range: A tuple to filter for uploads within with start, end ``upload_time``.
        search_parameters: Adds a ``and`` search for each key, value pair. Where the key corresponds
            to a quantity and the value is the value to search for in this quantity.

    Returns:
        A dict with keys 'pagination' and 'results' (similar to pagination in the REST API).
        The pagination key holds a dict with keys 'total', 'page', 'per_page'. The
        results key holds an array with the found entries.
    """
    search = _construct_search(q, time_range, search_parameters=search_parameters)
    _, results = _execute_paginated_search(search, page, per_page, order_by, order)

    return results


def quantity_search(
        quantities: Dict[str, Any], with_entries: bool = True, size: int = 100,
        **kwargs) -> Dict[str, Any]:
    """
    Performs a search like :func:`entry_search`, but instead of entries, returns the values
    of the given quantities that are exhibited by the entries in the search results.
    In contrast to :func:`metrics_search` it allows to scroll through all values via
    elasticsearch's composite aggregations.
    Optionally, it will also return the entries.

    This can be used to implement continues scrolling through authors, datasets, or uploads
    within the searched entries.

    Arguments:
        quantities: A dict, where the keys are quantity names, and the values are either
            None, or the 'after' value. This allows to scroll over various requests, by
            providing the 'after' value of the last search. The 'after' value is
            part of the return.
        with_entries: If True, the method will also return the entry search results. See
            :func:`entry_search`.
        size: The size of the quantity lists to return with each call.
        **kwargs: Additional arguments are passed to the underlying entry search.

    Returns:
        A dictionary with key 'quantities' (and optionally the keys of the
        return of :func:`entry_search` ). The 'quantities' key will hold a dict
        of quantities, each quantity is a dictionary with 'after' and 'values' key.
        The 'values' key holds a dict with actual values as keys and their entry count
        as values (i.e. number of entries with that value).
    """

    search = _construct_search(**kwargs)
    for quantity_name, after in quantities.items():
        quantity = search_quantities[quantity_name]
        terms = A('terms', field=quantity.elastic_field)

        composite = dict(sources={quantity_name: terms}, size=size)
        if after is not None:
            composite['after'] = {quantity_name: after}

        search.aggs.bucket(quantity_name, 'composite', **composite)

    response, entry_results = _execute_paginated_search(search, **kwargs)

    def create_quantity_result(quantity):
        values = getattr(response.aggregations, quantity)
        result = dict(values={
            getattr(bucket.key, quantity): bucket.doc_count
            for bucket in values.buckets})

        if hasattr(values, 'after_key'):
            result.update(after=getattr(values.after_key, quantity))

        return result

    quantity_results = {
        quantity: create_quantity_result(quantity)
        for quantity in quantities.keys()
    }

    results = dict(quantities=quantity_results)
    if with_entries:
        results.update(**entry_results)

    return results


def metrics_search(
        quantities: Dict[str, int] = aggregations, metrics_to_use: List[str] = [],
        with_entries: bool = True, **kwargs) -> Dict[str, Any]:
    """
    Performs a search like :func:`entry_search`, but instead of entries, returns the given
    metrics aggregated for (a limited set of values) of the given quantities calculated
    from the entries in the search results.
    In contrast to :func:`property_search` the amount of values for each quantity is
    limited.
    Optionally, it will also return the entries.

    This can be used to display statistics over the searched entries and allows to
    implement faceted search on the top values for each quantity.

    The metrics contain overall and per quantity value sums of code runs (calcs), unique code runs,
    datasets, and additional domain specific metrics (e.g. total energies, and unique geometries for DFT
    calculations). The quantities that can be aggregated to metrics are defined
    in module:`datamodel`. Aggregations and respective metrics are calculated for
    aggregations given in param:`aggregations` and metrics in param:`aggregation_metrics`.
    As a pseudo aggregation param:`total_metrics` are calculation over all search results.
    The param:`aggregations` gives tuples of quantities and default aggregation sizes.

    Arguments:
        aggregations: A customized list of aggregations to perform. Keys are index fields,
            and values the amount of buckets to return. Only works on *keyword* field.
        metrics_to_use: The metrics used to aggregate over. Can be ``unique_code_runs``, ``datasets``,
            other domain specific metrics. The basic doc_count metric ``code_runs`` is always given.
        **kwargs: Additional arguments are passed to the underlying entry search.

    Returns:
        A dictionary with key 'quantities' (and optionally the keys of the
        return of :func:`entry_search`). The 'quantities' key will hold a dict with a key
        for each quantity and an extra key 'total'. Each quantity key will hold a dict
        with a key for each quantity value. Each quantity value key will hold a dict
        with a key for each metric. The values will be the actual aggregated metric values.
        The pseudo quantity 'total' contains a pseudo value 'all'. It is used to
        store the metrics aggregated over all entries in the search results.
    """

    search = _construct_search(**kwargs)

    def add_metrics(parent):
        for metric in metrics_to_use:
            metric_kind, field = metrics[metric]
            parent.metric(metric, A(metric_kind, field=field))

    for quantity_name, size in quantities.items():
        # We are using elastic searchs 'composite aggregations' here. We do not really
        # compose aggregations, but only those pseudo composites allow us to use the
        # 'after' feature that allows to scan through all aggregation values.
        terms: Dict[str, Any] = None
        quantity = search_quantities[quantity_name]
        min_doc_count = 0 if quantity.zero_aggs else 1
        terms = A(
            'terms', field=quantity.elastic_field, size=size, min_doc_count=min_doc_count,
            order=dict(_key='asc'))

        buckets = search.aggs.bucket(quantity_name, terms)
        if quantity_name not in ['authors']:
            add_metrics(buckets)

    add_metrics(search.aggs)

    response, entry_results = _execute_paginated_search(search, **kwargs)

    def get_metrics(bucket, code_runs):
        result = {
            metric: bucket[metric]['value']
            for metric in metrics_to_use
            if hasattr(bucket, metric)
        }
        result.update(code_runs=code_runs)
        return result

    metrics_results = {
        quantity_name: {
            bucket.key: get_metrics(bucket, bucket.doc_count)
            for bucket in getattr(response.aggregations, quantity_name).buckets
        }
        for quantity_name in quantities.keys()
        if quantity_name not in metrics_names  # ES aggs for total metrics, and aggs for quantities stand side by side
    }

    total_metrics_result = get_metrics(response.aggregations, entry_results['pagination']['total'])
    metrics_results['total'] = dict(all=total_metrics_result)

    results = dict(quantities=metrics_results)
    if with_entries:
        results.update(**entry_results)

    return results
