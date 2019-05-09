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

from typing import Iterable, Dict, Tuple, List
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
            doi=dataset.doi.value if dataset.doi is not None else None,
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
            self.authors.append(self.uploader)
            self.owners.append(self.uploader)
        self.comment = source.comment
        self.references = [ref.value for ref in source.references]
        self.datasets = [Dataset.from_dataset_popo(ds) for ds in source.datasets]

        for quantity in datamodel.Domain.instance.quantities.keys():
            setattr(self, quantity, getattr(source, quantity))


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


def refresh():
    infrastructure.elastic_client.indices.refresh(config.elastic.index_name)


aggregations = datamodel.Domain.instance.aggregations
""" The available aggregations in :func:`aggregate_search` and their maximum aggregation size """

search_quantities = {
    'authors': ('term', 'authors.name.keyword', (
        'Search for the given author. Exact keyword matches in the form "Lastname, Firstname".')),

    'comment': ('match', 'comment', 'Search within the comments. This is a text search ala google.'),
    'paths': ('match', 'files', (
        'Search for elements in one of the file paths. The paths are split at all "/".')),

    'files': ('term', 'files.keyword', 'Search for exact file name with full path.'),
    'quantities': ('term', 'quantities', 'Search for the existence of a certain meta-info quantity'),
    'upload_id': ('term', 'upload_id', 'Search for the upload_id.'),
    'calc_id': ('term', 'calc_id', 'Search for the calc_id.'),
    'pid': ('term', 'pid', 'Search for the pid.'),
    'mainfile': ('term', 'mainfile', 'Search for the mainfile.')
}
"""
The available search quantities in :func:`aggregate_search` as tuples with *search type*,
elastic field and description.
"""

for quantity in datamodel.Domain.instance.quantities.values():
    search_spec = ('term', quantity.name, quantity.description)
    search_quantities[quantity.name] = search_spec


metrics = {
    'datasets': ('cardinality', 'datasets.id'),
    'unique_code_runs': ('cardinality', 'calc_hash')
}
"""
The available search metrics. Metrics are integer values given for each entry that can
be used in aggregations, e.g. the sum of all total energy calculations or cardinality of
all unique geometries.
"""

for key, value in datamodel.Domain.instance.metrics.items():
    metrics[key] = value

metrics_names = list(metric for metric in metrics.keys())


order_default_quantity = None
for quantity in datamodel.Domain.instance.quantities.values():
    if quantity.order_default:
        order_default_quantity = quantity.name


def _construct_search(q: Q = None, time_range: Tuple[datetime, datetime] = None, **kwargs) -> Search:
    search = Search(index=config.elastic.index_name)

    if q is not None:
        search = search.query(q)

    if time_range is not None:
        search = search.query('range', upload_time=dict(gte=time_range[0], lte=time_range[1]))

    for key, value in kwargs.items():
        query_type, field, _ = search_quantities.get(key, (None, None, None))
        if query_type is None:
            raise KeyError('Unknown quantity %s' % key)

        if isinstance(value, list):
            values = value
        else:
            values = [value]

        for item in values:
            quantity = datamodel.Domain.instance.quantities.get(key)
            if quantity is not None and quantity.multi:
                items = item.split(',')
            else:
                items = [item]

            for item in items:
                search = search.query(Q(query_type, **{field: item}))

    search = search.source(exclude=['quantities'])

    return search


def scroll_search(
        scroll_id: str = None, size: int = 1000, scroll: str = u'5m',
        q: Q = None, **kwargs) -> Tuple[str, int, List[dict]]:
    """
    Alternative search based on ES scroll API. Can be used similar to
    :func:`aggregate_search`, but pagination is replaced with scrolling, no ordering,
    and no aggregation information is given.

    Scrolling is done by calling this function again and again with the same ``scroll_id``.
    Each time, this function will return the next batch of search results. If the
    ``scroll_id`` is not available anymore, a new ``scroll_id`` is assigned and scrolling
    starts from the beginning again.

    See see :func:`aggregate_search` for additional ``kwargs``

    Arguments:
        scroll_id: The scroll id to receive the next batch from. None will create a new
            scroll.
        size: The batch size in number of hits.
        scroll: The time the scroll should be kept alive (i.e. the time between requests
            to this method) in ES time units. Default is 5 minutes.
    Returns: A tuple with ``scroll_id``, total amount of hits, and result list.
    """
    es = infrastructure.elastic_client

    if scroll_id is None:
        # initiate scroll
        search = _construct_search(q, **kwargs)
        resp = es.search(body=search.to_dict(), scroll=scroll, size=size, index=config.elastic.index_name)  # pylint: disable=E1123

        scroll_id = resp.get('_scroll_id')
        if scroll_id is None:
            # no results for search query
            return None, 0, []
    else:
        try:
            resp = es.scroll(scroll_id, scroll=scroll)  # pylint: disable=E1123
        except NotFoundError:
            raise ScrollIdNotFound()

    total = resp['hits']['total']
    results = [hit['_source'] for hit in resp['hits']['hits']]

    # since we are using the low level api here, we should check errors
    if resp["_shards"]["successful"] < resp["_shards"]["total"]:
        utils.get_logger(__name__).error('es operation was unsuccessful on at least one shard')
        raise ElasticSearchError('es operation was unsuccessful on at least one shard')

    if len(results) == 0:
        es.clear_scroll(body={'scroll_id': [scroll_id]}, ignore=(404, ))  # pylint: disable=E1123
        return None, total, []

    return scroll_id, total, results


def aggregate_search(
        page: int = 1, per_page: int = 10, order_by: str = order_default_quantity, order: int = -1,
        q: Q = None,
        time_range: Tuple[datetime, datetime] = None,
        aggregations: Dict[str, int] = aggregations,
        aggregation_metrics: List[str] = [],
        total_metrics: List[str] = [],
        **kwargs) -> Tuple[int, List[dict], Dict[str, Dict[str, Dict[str, int]]], Dict[str, int]]:
    """
    Performs a search and returns paginated search results and aggregations. The aggregations
    contain overall and per quantity value sums of code runs (calcs), unique code runs, datasets,
    and additional domain specific metrics (e.g. total energies, and unique geometries for DFT
    calculations).

    Arguments:
        page: The page to return starting with page 1
        per_page: Results per page
        q: An *elasticsearch_dsl* query used to further filter the results (via ``and``)
        time_range: A tuple to filter for uploads within with start, end ``upload_time``.
        aggregations: A customized list of aggregations to perform. Keys are index fields,
            and values the amount of buckets to return. Only works on *keyword* field.
        aggregation_metrics: The metrics used to aggregate over. Can be ``unique_code_runs``, ``datasets``,
            other domain specific metrics. The basic doc_count metric ``code_runs`` is always given.
        total_metrics: The metrics used to for total numbers (see ``aggregation_metrics``).
        **kwargs: Quantity, value pairs to search for.

    Returns: A tuple with the total hits, an array with the results, an dictionary with
        the aggregation data, and a dictionary with the overall metrics.
    """

    search = _construct_search(q, time_range, **kwargs)

    def add_metrics(parent, metrics_to_add):
        for metric in metrics_to_add:
            metric_kind, field = metrics[metric]
            parent.metric(metric, A(metric_kind, field=field))

    for aggregation, size in aggregations.items():

        if aggregation == 'authors':
            a = A('terms', field='authors.name_keyword', size=size)
        else:
            a = A('terms', field=aggregation, size=size, min_doc_count=0, order=dict(_key='asc'))

        buckets = search.aggs.bucket(aggregation, a)
        add_metrics(buckets, aggregation_metrics)

    add_metrics(search.aggs, total_metrics)

    if order_by not in search_quantities:
        raise KeyError('Unknown order quantity %s' % order_by)
    search = search.sort(order_by if order == 1 else '-%s' % order_by)

    response = search[(page - 1) * per_page: page * per_page].execute()  # pylint: disable=E1101

    total_results = response.hits.total
    search_results = [hit.to_dict() for hit in response.hits]

    def get_metrics(bucket, metrics_to_get, code_runs):
        result = {
            metric: bucket[metric]['value']
            for metric in metrics_to_get
        }
        result.update(code_runs=code_runs)
        return result

    aggregation_results = {
        aggregation: {
            bucket.key: get_metrics(bucket, aggregation_metrics, bucket.doc_count)
            for bucket in getattr(response.aggregations, aggregation).buckets
        }
        for aggregation in aggregations.keys()
        if aggregation not in metrics_names
    }

    total_metrics_result = get_metrics(response.aggregations, total_metrics, total_results)

    return total_results, search_results, aggregation_results, total_metrics_result


def authors(per_page: int = 10, after: str = None, prefix: str = None) -> Tuple[Dict[str, int], str]:
    """
    Returns the name field for all authors with the number of their calculations in
    their natural order.

    The author buckets of :func:`search` is limit to the top 10 author. This function
    in contrast, allows to paginate through all authors.

    Arguments:
        per_page: The number of authors to return
        after: Only return the authors after the given ``name_keyword`` field value
            (for pagination).
        prefix: Used to do a prefix search on authors. Be aware the return authors also
            contain the coauthors of the actual authors with prefix.

    Returns: A tuple with an ordered dict containing author ``name_keyword`` field value
        and number of calculations and the ``name_keyword`` value of the last author
        (to be used as the next ``after`` value).
    """
    composite = dict(
        size=per_page,
        sources=dict(authors=dict(terms=dict(field='authors.name.keyword'))))
    if after is not None:
        composite.update(after=dict(authors=after))

    body = dict(size=0, aggs=dict(authors=dict(composite=composite)))
    if prefix is not None:
        body.update(query=dict(prefix={'authors.name': prefix}))

    response = infrastructure.elastic_client.search(index=config.elastic.index_name, body=body)
    response = response['aggregations']['authors']
    return {
        bucket['key']['authors']: bucket['doc_count']
        for bucket in response['buckets']}, response.get('after_key', {'authors': None})['authors']
