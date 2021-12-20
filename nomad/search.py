#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

'''
This module represents calculations in elastic search.
'''

from typing import cast, Iterable, Dict, List, Any
from elasticsearch_dsl import Search, Q, A, analyzer, tokenizer
import elasticsearch.helpers
from elasticsearch.exceptions import NotFoundError, RequestError, TransportError
from datetime import datetime
import json

from nomad.datamodel.material import Material
from nomad import config, datamodel, infrastructure, utils
from nomad.metainfo.search_extension import (  # pylint: disable=unused-import
    search_quantities, metrics, order_default_quantities, groups)
from nomad.app.v1 import models as api_models
from nomad.app.v1.models import (
    EntryPagination, PaginationResponse, Query, MetadataRequired, SearchResponse, Aggregation,
    Statistic, StatisticResponse, AggregationOrderType, AggregationResponse, AggregationDataItem)


_entry_metadata_defaults = {
    quantity.name: quantity.default
    for quantity in datamodel.EntryMetadata.m_def.quantities  # pylint: disable=not-an-iterable
    if quantity.default not in [None, [], False, 0]
}

_all_author_quantities = [
    quantity.name
    for quantity in datamodel.EntryMetadata.m_def.all_quantities.values()
    if quantity.type in [datamodel.user_reference, datamodel.author_reference]]


def _es_to_entry_dict(hit, required: MetadataRequired) -> Dict[str, Any]:
    '''
    Elasticsearch entry metadata does not contain default values, if a metadata is not
    set. This will add default values to entry metadata in dict form obtained from
    elasticsearch.
    '''
    entry_dict = hit.to_dict()
    for key, value in _entry_metadata_defaults.items():
        if key not in entry_dict:
            if required is not None:
                if required.exclude and key in required.exclude:
                    continue
                if required.include and key not in required.include:
                    continue

            entry_dict[key] = value

    for author_quantity in _all_author_quantities:
        authors = entry_dict.get(author_quantity)
        if authors is None:
            continue
        if isinstance(authors, dict):
            authors = [authors]
        for author in authors:
            if 'email' in author:
                del(author['email'])

    return entry_dict


path_analyzer = analyzer(
    'path_analyzer',
    tokenizer=tokenizer('path_tokenizer', 'pattern', pattern='/'))


class AlreadyExists(Exception): pass


class ElasticSearchError(Exception): pass


class AuthenticationRequiredError(Exception): pass


class ScrollIdNotFound(Exception): pass


class InvalidQuery(Exception): pass


entry_document = datamodel.EntryMetadata.m_def.a_elastic.document
material_document = Material.m_def.a_elastic.document

for domain in datamodel.domains:
    order_default_quantities.setdefault(domain, order_default_quantities.get('__all__'))


def delete_upload(upload_id):
    ''' Delete all entries with given ``upload_id`` from the index. '''
    index = entry_document._default_index()
    Search(index=index).query('match', upload_id=upload_id).delete()


def delete_entry(calc_id):
    ''' Delete the entry with the given ``calc_id`` from the index. '''
    index = entry_document._default_index()
    Search(index=index).query('match', calc_id=calc_id).delete()


def publish(calcs: Iterable[datamodel.EntryMetadata]) -> None:
    ''' Update all given calcs with their metadata and set ``publish = True``. '''
    def elastic_updates():
        for calc in calcs:
            entry = calc.a_elastic.create_index_entry()
            entry.published = True
            entry = entry.to_dict(include_meta=True)
            source = entry.pop('_source')
            entry['doc'] = source
            entry['_op_type'] = 'update'
            yield entry

    elasticsearch.helpers.bulk(infrastructure.elastic_client, elastic_updates())
    refresh()


def index_all(calcs: Iterable[datamodel.EntryMetadata], do_refresh=True) -> None:
    '''
    Adds all given calcs with their metadata to the index.

    Returns:
        Number of failed entries.
    '''
    def elastic_updates():
        for calc in calcs:
            try:
                entry = calc.a_elastic.create_index_entry()
                entry = entry.to_dict(include_meta=True)
                entry['_op_type'] = 'index'

                yield entry

            except Exception as e:
                utils.get_logger(__name__).error(
                    'could not create index doc', exc_info=e,
                    upload_id=calc.upload_id, calc_id=calc.calc_id)

    _, failed = elasticsearch.helpers.bulk(infrastructure.elastic_client, elastic_updates(), stats_only=True)

    if do_refresh:
        refresh()

    return failed


def refresh():
    infrastructure.elastic_client.indices.refresh(config.elastic.index_name)


def _owner_es_query(owner: str, user_id: str = None):
    if owner == 'all':
        q = Q('term', published=True)
        if user_id is not None:
            q = q | Q('term', owners__user_id=user_id)
    elif owner == 'public':
        q = Q('term', published=True) & Q('term', with_embargo=False)
    elif owner == 'visible':
        q = Q('term', published=True) & Q('term', with_embargo=False)
        if user_id is not None:
            q = q | Q('term', owners__user_id=user_id)
    elif owner == 'shared':
        if user_id is None:
            raise AuthenticationRequiredError('Authentication required for owner value shared.')

        q = Q('term', owners__user_id=user_id)
    elif owner == 'user':
        if user_id is None:
            raise AuthenticationRequiredError('Authentication required for owner value user.')

        q = Q('term', uploader__user_id=user_id)
    elif owner == 'staging':
        if user_id is None:
            raise AuthenticationRequiredError('Authentication required for owner value user')
        q = Q('term', published=False) & Q('term', owners__user_id=user_id)
    elif owner == 'admin':
        if user_id is None or not datamodel.User.get(user_id=user_id).is_admin:
            raise AuthenticationRequiredError('This can only be used by the admin user.')
        q = None
    else:
        raise KeyError('Unsupported owner value')

    if q is not None:
        return q
    return Q()


class SearchRequest:
    '''
    Represents a search request and allows to execute that request.
    It allows to compose the following features: a query;
    statistics (metrics and aggregations); quantity values; scrolling, pagination for entries;
    scrolling for quantity values.

    The query part filters NOMAD data before the other features come into effect. There
    are specialized methods for configuring the :func:`owner` and :func:`time_range` queries.
    Quantity's can be search for by setting them as attributes.

    The aggregations for statistics can be requested for pre-configured quantities. These
    bucket aggregations come with a metric calculated for each each possible
    quantity value.

    The other possible form of aggregations, allows to get quantity values as results
    (e.g. get all datasets, get all users, etc.). Each value can be accompanied by metrics
    (over all entries with that value) and an example value.

    Of course, searches can return a set of search results. Search objects can be
    configured with pagination or scrolling for these results. Pagination is the default
    and also allows ordering of results. Scrolling can be used if all entries need to be
    'scrolled through'. This might be necessary, since elastic search has limits on
    possible pages (e.g. 'from' must by smaller than 10000). On the downside, there is no
    ordering on scrolling.

    There is also scrolling for quantities to go through all quantity values. There is no
    paging for aggregations.
    '''
    def __init__(self, domain: str = config.meta.default_domain, query=None):
        self._domain = domain
        self._query = query
        self._search = Search(index=config.elastic.index_name)
        self._required = None

    def domain(self, domain: str = None):
        '''
        Applies the domain of this request to the query. Allows to optionally update
        the domain of this request.
        '''
        if domain is not None:
            self._domain = domain

        self.q = self.q & Q('term', domain=self._domain)
        return self

    def owner(self, owner_type: str = 'all', user_id: str = None):
        '''
        Uses the query part of the search to restrict the results based on the owner.
        The possible types are: ``all`` for all calculations; ``public`` for
        calculations visible by everyone, excluding embargo-ed entries and entries only visible
        to the given user; ``visible`` all data that is visible by the user, excluding
        embargo-ed entries from other users; ``user`` for all calculations of to the given
        user; ``staging`` for all calculations in staging of the given user.

        Arguments:
            owner_type: The type of the owner query, see above.
            user_id: The 'owner' given as the user's unique id.

        Raises:
            KeyError: If the given owner_type is not supported
            ValueError: If the owner_type requires a user but none is given, or the
                given user is not allowed to use the given owner_type.
        '''
        self.q &= _owner_es_query(owner=owner_type, user_id=user_id)
        return self

    def search_parameters(self, **kwargs):
        '''
        Configures the existing query with additional search parameters. Kwargs are
        interpreted as key value pairs. Keys have to coresspond to valid entry quantities
        in the domain's (DFT calculations) datamodel. Alternatively search parameters
        can be set via attributes.
        '''
        for name, value in kwargs.items():
            self.search_parameter(name, value)

        return self

    def _search_parameter_to_es(self, name, value):
        quantity = search_quantities[name]

        if quantity.many and not isinstance(value, list):
            value = [value]

        if quantity.many_or and isinstance(value, List):
            return Q('terms', **{quantity.search_field: value})

        if quantity.derived:
            if quantity.many and not isinstance(value, list):
                value = [value]
            value = quantity.derived(value)

        if isinstance(value, list):
            values = value
        else:
            values = [value]

        return Q('bool', must=[
            Q('match', **{quantity.search_field: item})
            for item in values])

    def search_parameter(self, name, value):
        self.q &= self._search_parameter_to_es(name, value)
        return self

    def query(self, query):
        ''' Adds the given query as a 'and' (i.e. 'must') clause to the request. '''
        self._query &= query

        return self

    def query_expression(self, expression) -> 'SearchRequest':

        bool_operators = {'$and': 'must', '$or': 'should', '$not': 'must_not'}
        comp_operators = {'$%s' % op: op for op in ['gt', 'gte', 'lt', 'lte']}

        def _to_es(key, values):
            if key in bool_operators:
                if isinstance(values, dict):
                    values = [values]
                assert isinstance(values, list), 'bool operator requires a list of dicts or dict'
                child_es_queries = [
                    _to_es(child_key, child_value)
                    for child_query in values
                    for child_key, child_value in child_query.items()]
                return Q('bool', **{bool_operators[key]: child_es_queries})

            if key in comp_operators:
                assert isinstance(values, dict), 'comparison operator requires a dict'
                assert len(values) == 1, 'comparison operator requires exactly one quantity'
                quantity_name, value = next(iter(values.items()))
                quantity = search_quantities.get(quantity_name)
                assert quantity is not None, 'quantity %s does not exist' % quantity_name
                return Q('range', **{quantity.search_field: {comp_operators[key]: value}})

            try:
                return self._search_parameter_to_es(key, values)
            except KeyError:
                assert False, 'quantity %s does not exist' % key

        if len(expression) == 0:
            self.q &= Q()
        else:
            self.q &= Q('bool', must=[_to_es(key, value) for key, value in expression.items()])

        return self

    def time_range(self, start: datetime, end: datetime):
        ''' Adds a time range to the query. '''
        if start is None and end is None:
            return self

        if start is None:
            start = datetime.fromtimestamp(0)
        if end is None:
            end = datetime.utcnow()

        self.q &= Q('range', upload_time=dict(gte=start, lte=end))

        return self

    @property
    def q(self):
        ''' The underlying elasticsearch_dsl query object '''
        if self._query is None:
            return Q('match_all')
        else:
            return self._query

    @q.setter
    def q(self, q):
        self._query = q

    def totals(self, metrics_to_use: List[str] = []):
        '''
        Configure the request to return overall totals for the given metrics.

        The statics are returned with the other quantity statistics under the pseudo
        quantity name 'total'. 'total' contains the pseudo value 'all'. It is used to
        store the metrics aggregated over all entries in the search results.
        '''
        self._add_metrics(self._search.aggs, metrics_to_use)
        return self

    def statistics(self, statistics: List[str], metrics_to_use: List[str] = []):
        '''
        Configures the domain's default statistics.
        '''
        for statistic in statistics:
            search_quantity = search_quantities[statistic]
            statistic_order = search_quantity.statistic_order
            self.statistic(
                search_quantity.qualified_name,
                search_quantity.statistic_size,
                metrics_to_use=metrics_to_use,
                order={statistic_order: 'asc' if statistic_order == '_key' else 'desc'})

        return self

    def statistic(
            self, quantity_name: str, size: int, metrics_to_use: List[str] = [],
            order: Dict[str, str] = dict(_key='asc'), include: str = None):
        '''
        This can be used to display statistics over the searched entries and allows to
        implement faceted search on the top values for each quantity.

        The metrics contain overall and per quantity value sums of code runs (calcs),
        unique code runs, datasets, and additional domain specific metrics
        (e.g. total energies, and unique geometries for DFTcalculations). The quantities
        that can be aggregated to metrics are defined in module:`datamodel`. Aggregations
        and respective metrics are calculated for aggregations given in ``aggregations``
        and metrics in ``aggregation_metrics``. As a pseudo aggregation ``total_metrics``
        are calculation over all search results. The ``aggregations`` gives tuples of
        quantities and default aggregation sizes.

        The search results will contain a dictionary ``statistics``. This has a key
        for each configured quantity. Each quantity key will hold a dict
        with a key for each quantity value. Each quantity value key will hold a dict
        with a key for each metric. The values will be the actual aggregated metric values.

        Arguments:
            quantity_name: The quantity to aggregate statistics for. Only works on *keyword* field.
            metrics_to_use: The metrics calculated over the aggregations. Can be
                ``unique_code_runs``, ``datasets``, other domain specific metrics.
                The basic doc_count metric ``code_runs`` is always given.
            order: The order dictionary is passed to the elastic search aggregation.
            include:
                Uses an regular expression in ES to only return values that include
                the given substring.
        '''
        quantity = search_quantities[quantity_name]
        terms_kwargs = {}
        if include is not None:
            terms_kwargs['include'] = '.*%s.*' % include
        terms = A('terms', field=quantity.search_field, size=size, order=order, **terms_kwargs)

        buckets = self._search.aggs.bucket('statistics:%s' % quantity_name, terms)
        self._add_metrics(buckets, metrics_to_use)

        return self

    def _add_metrics(self, parent=None, metrics_to_use: List[str] = []):
        if parent is None:
            parent = self._search.aggs

        for metric in metrics_to_use:
            metric_quantity = metrics[metric]
            field = metric_quantity.search_field
            parent.metric(
                'metric:%s' % metric_quantity.metric_name,
                A(metric_quantity.metric, field=field))

    def date_histogram(self, metrics_to_use: List[str] = [], interval: str = '1M'):
        '''
        Adds a date histogram on the given metrics to the statistics part.
        '''
        histogram = A('date_histogram', field='upload_time', interval=interval, format='yyyy-MM-dd')
        self._add_metrics(self._search.aggs.bucket('statistics:date_histogram', histogram), metrics_to_use)

        return self

    def quantities(self, **kwargs):
        '''
        Shorthand for adding multiple quantities. See :func:`quantity`. Keywork argument
        keys are quantity name, values are tuples of size and after value.
        '''
        for name, spec in kwargs:
            size, after = spec
            self.quantity(name, after=after, size=size)

        return self

    def quantity(
            self, name, size=100, after=None, examples=0, examples_source=None,
            order_by: str = None, order: str = 'desc'):
        '''
        Adds a requests for values of the given quantity.
        It allows to scroll through all values via elasticsearch's
        composite aggregations. The response will contain the quantity values and
        an example entry for each value.

        This can be used to implement continues scrolling through authors, datasets,
        or uploads within the searched entries.

        If one or more quantities are specified,
        the search results will contain a dictionary ``quantities``. The keys are quantity
        name the values dictionary with 'after' and 'values' key.
        The 'values' key holds a dict with all the values as keys and their entry count
        as values (i.e. number of entries with that value).

        Arguments:
            name: The quantity name. Must be in :data:`quantities`.
            after: The 'after' value allows to scroll over various requests, by providing
                the 'after' value of the last search. The 'after' value is part of the
                response. Use ``None`` in the first request.
            size:
                The size gives the ammount of maximum values in the next scroll window.
                If the size is None, a maximum of 100 quantity values will be requested.
            examples:
                Number of results to return that has each value
            order_by:
                A sortable quantity that should be used to order. By default, the max of each
                value bucket is used.
            order:
                "desc" or "asc"
        '''
        if size is None:
            size = 100

        quantity = search_quantities[name]
        terms = A('terms', field=quantity.search_field)

        # We are using elastic searchs 'composite aggregations' here. We do not really
        # compose aggregations, but only those pseudo composites allow us to use the
        # 'after' feature that allows to scan through all aggregation values.
        if order_by is None:
            composite = dict(sources={name: terms}, size=size)
        else:
            sort_terms = A('terms', field=order_by, order=order)
            composite = dict(sources=[{order_by: sort_terms}, {name: terms}], size=size)
        if after is not None:
            if order_by is None:
                composite['after'] = {name: after}
            else:
                composite['after'] = {order_by: after, name: ''}

        composite_agg = self._search.aggs.bucket('quantity:%s' % name, 'composite', **composite)

        if examples > 0:
            kwargs: Dict[str, Any] = {}
            if examples_source is not None:
                kwargs.update(_source=dict(includes=examples_source))

            composite_agg.metric('examples', A('top_hits', size=examples, **kwargs))

        return self

    def global_statistics(self):
        '''
        Adds general statistics to the request. The results will have a key called
        global_statistics.
        '''
        self.owner('public')
        self._search.aggs.metric(
            'global_statistics:n_entries', A('value_count', field='calc_id'))
        self._search.aggs.metric(
            'global_statistics:n_uploads', A('cardinality', field='upload_id'))
        self._search.aggs.metric(
            'global_statistics:n_calculations', A('sum', field='dft.n_calculations'))
        self._search.aggs.metric(
            'global_statistics:n_quantities', A('sum', field='dft.n_quantities'))
        self._search.aggs.metric(
            'global_statistics:n_materials', A('cardinality', field='encyclopedia.material.material_id'))

        return self

    def exclude(self, *args):
        ''' Exclude certain elastic fields from the search results. '''
        self._search = self._search.source(excludes=args)
        self._required = MetadataRequired(exclude=args)
        return self

    def include(self, *args):
        ''' Include only the given fields in the search results. '''
        self._search = self._search.source(includes=args)
        return self

    def execute(self):
        '''
        Executes without returning actual results. Only makes sense if the request
        was configured for statistics or quantity values.
        '''
        search = self._search.query(self.q)[0:0]
        response = search.execute()
        return self._response(response)

    def execute_scan(self, order_by: str = None, order: int = -1, **kwargs):
        '''
        This execute the search as scan. The result will be a generator over the found
        entries. Everything but the query part of this object, will be ignored.
        '''
        search = self._search.query(self.q)

        if order_by is not None:
            order_by_quantity = search_quantities[order_by]

            if order == 1:
                search = search.sort(order_by_quantity.search_field)
            else:
                search = search.sort('-%s' % order_by_quantity.search_field)  # pylint: disable=no-member

            search = search.params(preserve_order=True)

        for hit in search.params(**kwargs).scan():
            yield _es_to_entry_dict(hit, self._required)

    def execute_paginated(
            self,
            page: int = 1, per_page=10, page_offset: int = None,
            order_by: str = None, order: int = -1):
        '''
        Executes the search and returns paginated results. Those are sorted.

        Arguments:
            page: The requested page, starts with 1.
            per_page: The number of entries per page.
            page_offset: Instead of a page number, use this absolute offset.
            order_by: The quantity to order by.
            order: -1 or 1 for descending or ascending order.
        '''
        if order_by is None:
            order_by_quantity = order_default_quantities[self._domain]
        else:
            order_by_quantity = search_quantities[order_by]

        search = self._search.query(self.q)

        if order == 1:
            search = search.sort(order_by_quantity.search_field)
        else:
            search = search.sort('-%s' % order_by_quantity.search_field)  # pylint: disable=no-member

        if page_offset is not None:
            search = search[page_offset: page_offset + per_page]  # pylint: disable=unsubscriptable-object
        else:
            search = search[(page - 1) * per_page: page * per_page]  # pylint: disable=unsubscriptable-object

        es_result = search.execute()

        result = self._response(es_result, with_hits=True)

        if page_offset is not None:
            result.update(pagination=dict(
                total=result['total'],
                page_offset=page_offset, per_page=per_page))
        else:
            result.update(
                pagination=dict(total=result['total'], page=page, per_page=per_page))
        return result

    def execute_scrolled(
            self, scroll_id: str = None, size: int = 1000, scroll: str = u'5m',
            order_by: str = None, order: int = -1):
        '''
        Executes a scrolling search. based on ES scroll API. Pagination is replaced with
        scrolling, no ordering is available, no statistics, no quantities will be provided.

        Scrolling is done by calling this function again and again with the same ``scroll_id``.
        Each time, this function will return the next batch of search results. If the
        ``scroll_id`` is not available anymore, a new ``scroll_id`` is assigned and scrolling
        starts from the beginning again.

        The response will contain a 'scroll' part with attributes 'total', 'scroll_id',
        and 'size'.

        Arguments:
            scroll_id: The scroll id to receive the next batch from. None will create a new
                scroll.
            size: The batch size in number of hits.
            scroll: The time the scroll should be kept alive (i.e. the time between requests
                to this method) in ES time units. Default is 5 minutes.

        TODO support order and order_by
        '''
        es = infrastructure.elastic_client

        if scroll_id is None:
            # initiate scroll
            resp = es.search(  # pylint: disable=E1123
                body=self._search.query(self.q).to_dict(), scroll=scroll, size=size,
                index=config.elastic.index_name)

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

        scroll_info = dict(total=total, size=size, scroll=True)
        if scroll_id is not None:
            scroll_info.update(scroll_id=scroll_id)

        return dict(scroll=scroll_info, results=results)

    def execute_aggregated(
            self, after: str = None, per_page: int = 1000, includes: List[str] = None):
        '''
        Uses a composite aggregation on top of the search to go through the result
        set. This allows to go arbirarely deep without using scroll. But, it will
        only return results with ``upload_id``, ``calc_id`` and the given
        quantities. The results will be 'ordered' by ``upload_id``.

        Arguments:
            after: The key that determines the start of the current page. This after
                key is returned with each response. Use None (default) for the first
                request.
            per_page: The size of each page.
            includes: A list of quantity names that should be returned in addition to
                ``upload_id`` and ``calc_id``.
        '''
        upload_id_agg = A('terms', field="upload_id")
        calc_id_agg = A('terms', field="calc_id")

        composite = dict(
            sources=[dict(upload_id=upload_id_agg), dict(calc_id=calc_id_agg)],
            size=per_page)

        if after is not None:
            upload_id, calc_id = after.split(':')
            composite['after'] = dict(upload_id=upload_id, calc_id=calc_id)

        composite_agg = self._search.aggs.bucket('ids', 'composite', **composite)
        if includes is not None:
            composite_agg.metric('examples', A('top_hits', size=1, _source=dict(includes=includes)))

        search = self._search.query(self.q)[0:0]
        response = search.execute()

        ids = response['aggregations']['ids']
        if 'after_key' in ids:
            after_dict = ids['after_key']
            after = '%s:%s' % (after_dict['upload_id'], after_dict['calc_id'])
        else:
            after = None

        id_agg_info = dict(total=response['hits']['total'], after=after, per_page=per_page)

        def transform_result(es_result):
            result = dict(
                upload_id=es_result['key']['upload_id'],
                calc_id=es_result['key']['calc_id'])

            if includes is not None:
                source = es_result['examples']['hits']['hits'][0]['_source']
                for key in source:
                    result[key] = source[key]

            return result

        results = [
            transform_result(item) for item in ids['buckets']]

        return dict(aggregation=id_agg_info, results=results)

    def _response(self, response, with_hits: bool = False) -> Dict[str, Any]:
        '''
        Prepares a response object covering the total number of results, hits, statistics,
        and quantities. Other aspects like pagination and scrolling have to be added
        elsewhere.
        '''
        result: Dict[str, Any] = dict()
        aggs = response.aggregations.to_dict()

        # total
        total = response.hits.total if hasattr(response, 'hits') else 0
        result.update(total=total)

        # hits
        if len(response.hits) > 0 or with_hits:
            result.update(results=[_es_to_entry_dict(hit, self._required) for hit in response.hits])

        # statistics
        def get_metrics(bucket, code_runs):
            result = {}
            # TODO optimize ... go through the buckets not the metrics
            for metric in metrics:
                agg_name = 'metric:%s' % metric
                if agg_name in bucket:
                    result[metric] = bucket[agg_name]['value']
            result.update(code_runs=code_runs)
            return result

        statistics_results = {
            quantity_name[11:]: {
                str(bucket['key']): get_metrics(bucket, bucket['doc_count'])
                for bucket in quantity['buckets']
            }
            for quantity_name, quantity in aggs.items()
            if quantity_name.startswith('statistics:')
        }

        # global statistics
        global_statistics_results = {
            agg_name[18:]: agg.get('value')
            for agg_name, agg in aggs.items()
            if agg_name.startswith('global_statistics:')
        }
        if len(global_statistics_results) > 0:
            result.update(global_statistics=global_statistics_results)

        # totals
        totals_result = get_metrics(aggs, total)
        statistics_results['total'] = dict(all=totals_result)

        if len(statistics_results) > 0:
            result.update(statistics=statistics_results)

        # quantities
        def create_quantity_result(quantity_name, quantity):
            values = {}
            for bucket in quantity['buckets']:
                value = dict(
                    total=bucket['doc_count'])
                if 'examples' in bucket:
                    examples = [hit['_source'] for hit in bucket['examples']['hits']['hits']]
                    value.update(examples=examples)

                values[bucket['key'][quantity_name]] = value

            result = dict(values=values)
            if 'after_key' in quantity:
                after = quantity['after_key']
                if len(after) == 1:
                    result.update(after=after[quantity_name])
                else:
                    for key in after:
                        if key != quantity_name:
                            result.update(after=after[key])
                            break

            return result

        quantity_results = {
            quantity_name[9:]: create_quantity_result(quantity_name[9:], quantity)
            for quantity_name, quantity in aggs.items()
            if quantity_name.startswith('quantity:')
        }

        if len(quantity_results) > 0:
            result.update(quantities=quantity_results)

        return result

    def __str__(self):
        return json.dumps(self._search.to_dict(), indent=2)


def flat(obj, prefix=None):
    '''
    Helper that translates nested result objects into flattened dicts with
    ``domain.quantity`` as keys.
    '''
    if isinstance(obj, dict):
        result = {}
        for key, value in obj.items():
            if isinstance(value, dict):
                value = flat(value)
                for child_key, child_value in value.items():
                    result['%s.%s' % (key, child_key)] = child_value

            else:
                result[key] = value

        return result
    else:
        return obj


def _api_to_es_query(query: api_models.Query) -> Q:
    '''
    Creates an ES query based on the API's query model. This needs to be a normalized
    query expression with explicit objects for logical, set, and comparison operators.
    Shorthand notations ala ``quantity:operator`` are not supported here; this
    needs to be resolved via the respective pydantic validator. There is also no
    validation of quantities and types.
    '''
    def quantity_to_es(name: str, value: api_models.Value) -> Q:
        # TODO depends on keyword or not, value might need normalization, etc.
        quantity = search_quantities[name]
        return Q('match', **{quantity.search_field: value})

    def parameter_to_es(name: str, value: api_models.QueryParameterValue) -> Q:

        if isinstance(value, api_models.All):
            return Q('bool', must=[
                quantity_to_es(name, item)
                for item in value.op])

        if isinstance(value, api_models.Any_):
            return Q('bool', should=[
                quantity_to_es(name, item)
                for item in value.op])

        if isinstance(value, api_models.None_):
            return Q('bool', must_not=[
                quantity_to_es(name, item)
                for item in value.op])

        if isinstance(value, api_models.ComparisonOperator):
            quantity = search_quantities[name]
            return Q('range', **{quantity.search_field: {
                type(value).__name__.lower(): value.op}})

        # list of values is treated as an "all" over the items
        if isinstance(value, list):
            return Q('bool', must=[
                quantity_to_es(name, item)
                for item in value])

        return quantity_to_es(name, value)

    def query_to_es(query: api_models.Query) -> Q:
        if isinstance(query, api_models.LogicalOperator):
            if isinstance(query, api_models.And):
                return Q('bool', must=[query_to_es(operand) for operand in query.op])

            if isinstance(query, api_models.Or):
                return Q('bool', should=[query_to_es(operand) for operand in query.op])

            if isinstance(query, api_models.Not):
                return Q('bool', must_not=query_to_es(query.op))

            raise NotImplementedError()

        if not isinstance(query, dict):
            raise NotImplementedError()

        # dictionary is like an "and" of all items in the dict
        if len(query) == 0:
            return Q()

        if len(query) == 1:
            key = next(iter(query))
            return parameter_to_es(key, query[key])

        return Q('bool', must=[
            parameter_to_es(name, value) for name, value in query.items()])

    return query_to_es(query)


def _api_to_es_statistic(es_search: Search, name: str, statistic: Statistic) -> A:
    '''
    Creates an ES aggregation based on the API's statistic model.
    '''

    quantity = search_quantities[statistic.quantity.value]
    if quantity.statistic_values is not None:
        statistic.size = len(quantity.statistic_values)

    terms_kwargs = {}
    if statistic.value_filter is not None:
        terms_kwargs['include'] = '.*%s.*' % statistic.value_filter

    order_type = '_count' if statistic.order.type_ == AggregationOrderType.entries else '_key'
    statistic_agg = es_search.aggs.bucket('statistic:%s' % name, A(
        'terms',
        field=quantity.search_field,
        size=statistic.size,
        order={order_type: statistic.order.direction.value},
        **terms_kwargs))

    for metric in statistic.metrics:
        metric_quantity = metrics[metric.value]
        statistic_agg.metric('metric:%s' % metric_quantity.metric_name, A(
            metric_quantity.metric,
            field=metric_quantity.search_field))


def _es_to_api_statistics(es_response, name: str, statistic: Statistic) -> StatisticResponse:
    '''
    Creates a StatisticResponse from elasticsearch response on a request executed with
    the given statistics.
    '''
    quantity = search_quantities[statistic.quantity.value]

    es_statistic = es_response.aggs['statistic:' + name]
    statistic_data = {}
    for bucket in es_statistic.buckets:
        value_data = dict(entries=bucket.doc_count)
        for metric in statistic.metrics:
            value_data[metric.value] = bucket['metric:' + metric.value].value
        statistic_data[bucket.key] = value_data

    if quantity.statistic_values is not None:
        for value in quantity.statistic_values:
            if value not in statistic_data:
                statistic_data[value] = dict(entries=0, **{
                    metric.value: 0 for metric in statistic.metrics})

    return StatisticResponse(data=statistic_data, **statistic.dict(by_alias=True))


def _api_to_es_aggregation(es_search: Search, name: str, agg: Aggregation) -> A:
    '''
    Creates an ES aggregation based on the API's aggregation model.
    '''
    quantity = search_quantities[agg.quantity.value]
    terms = A('terms', field=quantity.search_field, order=agg.pagination.order.value)

    # We are using elastic searchs 'composite aggregations' here. We do not really
    # compose aggregations, but only those pseudo composites allow us to use the
    # 'after' feature that allows to scan through all aggregation values.
    order_by = agg.pagination.order_by
    if order_by is None:
        composite = dict(sources={name: terms}, size=agg.pagination.page_size)
    else:
        order_quantity = search_quantities[order_by]
        sort_terms = A('terms', field=order_quantity.search_field, order=agg.pagination.order.value)
        composite = dict(sources=[{order_by: sort_terms}, {quantity.name: terms}], size=agg.pagination.page_size)

    if agg.pagination.page_after_value is not None:
        if order_by is None:
            composite['after'] = {name: agg.pagination.page_after_value}
        else:
            order_value, quantity_value = agg.pagination.page_after_value.split(':')
            composite['after'] = {quantity.name: quantity_value, order_quantity.name: order_value}

    composite_agg = es_search.aggs.bucket('agg:%s' % name, 'composite', **composite)

    if agg.entries is not None and agg.entries.size > 0:
        kwargs: Dict[str, Any] = {}
        if agg.entries.required is not None:
            if agg.entries.required.include is not None:
                kwargs.update(_source=dict(includes=agg.entries.required.include))
            else:
                kwargs.update(_source=dict(excludes=agg.entries.required.exclude))

        composite_agg.metric('entries', A('top_hits', size=agg.entries.size, **kwargs))

    # additional cardinality to get total
    es_search.aggs.metric('agg:%s:total' % name, 'cardinality', field=quantity.search_field)


def _es_to_api_aggregation(es_response, name: str, agg: Aggregation) -> AggregationResponse:
    '''
    Creates a AggregationResponse from elasticsearch response on a request executed with
    the given aggregation.
    '''
    order_by = agg.pagination.order_by
    quantity = search_quantities[agg.quantity.value]
    es_agg = es_response.aggs['agg:' + name]

    def get_entries(agg):
        if 'entries' in agg:
            return [item['_source'] for item in agg.entries.hits.hits]
        else:
            return None

    if agg.pagination.order_by is None:
        agg_data = {
            bucket.key[name]: AggregationDataItem(size=bucket.doc_count, data=get_entries(bucket))
            for bucket in es_agg.buckets}
    else:
        agg_data = {
            bucket.key[quantity.search_field]: AggregationDataItem(size=bucket.doc_count, data=get_entries(bucket))
            for bucket in es_agg.buckets}

    aggregation_dict = agg.dict(by_alias=True)
    pagination = PaginationResponse(
        total=es_response.aggs['agg:%s:total' % name]['value'],
        **aggregation_dict.pop('pagination'))

    if 'after_key' in es_agg:
        after_key = es_agg['after_key']
        if order_by is None:
            pagination.next_page_after_value = after_key[name]
        else:
            str_values = [str(v) for v in after_key.to_dict().values()]
            pagination.next_page_after_value = ':'.join(str_values)

    return AggregationResponse(data=agg_data, pagination=pagination, **aggregation_dict)


def search(
        owner: str = 'public',
        query: Query = None,
        pagination: EntryPagination = None,
        required: MetadataRequired = None,
        aggregations: Dict[str, Aggregation] = {},
        statistics: Dict[str, Statistic] = {},
        user_id: str = None) -> SearchResponse:

    # The first half of this method creates the ES query. Then the query is run on ES.
    # The second half is about transforming the ES response to a SearchResponse.

    # query and owner
    if query is None:
        query = {}
    es_query = _api_to_es_query(query)
    if owner is not None:
        es_query &= _owner_es_query(owner=owner, user_id=user_id)

    # pagination
    if pagination is None:
        pagination = EntryPagination()

    search = Search(index=config.elastic.index_name)

    search = search.query(es_query)
    order_field = search_quantities[pagination.order_by].search_field
    sort = {order_field: pagination.order.value}
    if order_field != 'calc_id':
        sort['calc_id'] = pagination.order.value
    search = search.sort(sort)
    search = search.extra(size=pagination.page_size)
    if pagination.page_after_value:
        search = search.extra(search_after=pagination.page_after_value.rsplit(':', 1))

    if pagination.page:
        search = search[(pagination.page - 1) * pagination.page_size: pagination.page * pagination.page_size]

    # required
    if required:
        if required.include is not None and pagination.order_by not in required.include:
            required.include.append(pagination.order_by)
        if required.exclude is not None and pagination.order_by in required.exclude:
            required.exclude.remove(pagination.order_by)
        search = search.source(includes=required.include, excludes=required.exclude)

    # statistics
    for name, statistic in statistics.items():
        _api_to_es_statistic(search, name, statistic)

    # aggregations
    for name, agg in aggregations.items():
        _api_to_es_aggregation(search, name, agg)

    # execute
    try:
        es_response = search.execute()
    except RequestError as e:
        raise ElasticSearchError(e)

    more_response_data = {}

    # pagination
    next_page_after_value = None
    if 0 < len(es_response.hits) < es_response.hits.total:
        last = es_response.hits[-1]
        if order_field == 'calc_id':
            next_page_after_value = last['calc_id']
        else:
            after_value = last
            for order_field_segment in order_field.split('.'):
                after_value = after_value[order_field_segment]
            next_page_after_value = '%s:%s' % (after_value, last['calc_id'])
    pagination_response = PaginationResponse(
        total=es_response.hits.total,
        next_page_after_value=next_page_after_value,
        **pagination.dict())

    # statistics
    if len(statistics) > 0:
        more_response_data['statistics'] = cast(Dict[str, Any], {
            name: _es_to_api_statistics(es_response, name, statistic)
            for name, statistic in statistics.items()})

    # aggregations
    if len(aggregations) > 0:
        more_response_data['aggregations'] = cast(Dict[str, Any], {
            name: _es_to_api_aggregation(es_response, name, aggregation)
            for name, aggregation in aggregations.items()})

    more_response_data['es_query'] = es_query.to_dict()

    return SearchResponse(
        owner=owner,
        query=query,
        pagination=pagination_response,
        required=required,
        data=[_es_to_entry_dict(hit, required) for hit in es_response.hits],
        **more_response_data)


def update_by_query(update_script: str, owner: str = 'public', query: Query = None, user_id: str = None, **kwargs):
    '''
    Uses the given painless script to update the entries by given query.

    In most cases, the elasticsearch entry index should not be updated field by field;
    you should run `index_all` instead and fully replace documents from mongodb and
    archive files.

    This method provides a faster direct method to update individual fiels, e.g. to quickly
    update fields for editing operations.
    '''

    if query is None:
        query = {}
    es_query = _api_to_es_query(query)
    es_query &= _owner_es_query(owner=owner, user_id=user_id)

    body = {
        'script': {
            'source': update_script,
            'lang': 'painless'
        },
        'query': es_query.to_dict()
    }

    body['script'].update(**kwargs)

    try:
        result = infrastructure.elastic_client.update_by_query(
            body=body, index=config.elastic.index_name)
    except TransportError as e:
        utils.get_logger(__name__).error(
            'es update_by_query script error', exc_info=e,
            es_info=json.dumps(e.info, indent=2))
        raise e

    return result
