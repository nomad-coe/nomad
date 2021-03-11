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

from typing import cast, Any, Dict
from elasticsearch.exceptions import RequestError, TransportError
from elasticsearch_dsl import Search, Q, A, analyzer, tokenizer
import json

from nomad import config, infrastructure, utils, datamodel
from nomad.search import _owner_es_query
from nomad.metainfo.elasticsearch_extension import entry_type

from . import models as api_models
from .models import (
    Pagination, PaginationResponse, Query, MetadataRequired, SearchResponse, Aggregation,
    Statistic, StatisticResponse, AggregationOrderType, AggregationResponse, AggregationDataItem)


class SearchError(Exception): pass


_entry_metadata_defaults = {
    quantity.name: quantity.default
    for quantity in datamodel.EntryMetadata.m_def.quantities  # pylint: disable=not-an-iterable
    if quantity.default not in [None, [], False, 0]
}


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

    return entry_dict


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
        quantity = entry_type.quantities[name]
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
            quantity = entry_type.quantities[name]
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

    quantity = entry_type.quantities[statistic.quantity.value]
    if quantity.values is not None:
        statistic.size = len(quantity.values)

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
        metric_name = metric.value
        metric_aggregation, metric_quantity = entry_type.metrics[metric_name]
        statistic_agg.metric('metric:%s' % metric_name, A(
            metric_aggregation,
            field=metric_quantity.qualified_field))


def _es_to_api_statistics(es_response, name: str, statistic: Statistic) -> StatisticResponse:
    '''
    Creates a StatisticResponse from elasticsearch response on a request executed with
    the given statistics.
    '''
    quantity = entry_type.quantities[statistic.quantity.value]

    es_statistic = es_response.aggs['statistic:' + name]
    statistic_data = {}
    for bucket in es_statistic.buckets:
        value_data = dict(entries=bucket.doc_count)
        for metric in statistic.metrics:
            value_data[metric.value] = bucket['metric:' + metric.value].value
        statistic_data[bucket.key] = value_data

    if quantity.values is not None:
        for value in quantity.values:
            if value not in statistic_data:
                statistic_data[value] = dict(entries=0, **{
                    metric.value: 0 for metric in statistic.metrics})

    return StatisticResponse(data=statistic_data, **statistic.dict(by_alias=True))


def _api_to_es_aggregation(es_search: Search, name: str, agg: Aggregation) -> A:
    '''
    Creates an ES aggregation based on the API's aggregation model.
    '''
    quantity = entry_type.quantities[agg.quantity.value]
    terms = A('terms', field=quantity.search_field, order=agg.pagination.order.value)

    # We are using elastic searchs 'composite aggregations' here. We do not really
    # compose aggregations, but only those pseudo composites allow us to use the
    # 'after' feature that allows to scan through all aggregation values.
    order_by = agg.pagination.order_by
    if order_by is None:
        composite = dict(sources={name: terms}, size=agg.pagination.size)
    else:
        order_quantity = entry_type.quantities[order_by]
        sort_terms = A('terms', field=order_quantity.search_field, order=agg.pagination.order.value)
        composite = dict(sources=[{order_by: sort_terms}, {quantity.name: terms}], size=agg.pagination.size)

    if agg.pagination.after is not None:
        if order_by is None:
            composite['after'] = {name: agg.pagination.after}
        else:
            order_value, quantity_value = agg.pagination.after.split(':')
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
    quantity = entry_type.quantities[agg.quantity.value]
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
            pagination.next_after = after_key[name]
        else:
            str_values = [str(v) for v in after_key.to_dict().values()]
            pagination.next_after = ':'.join(str_values)

    return AggregationResponse(data=agg_data, pagination=pagination, **aggregation_dict)


def search(
        owner: str = 'public',
        query: Query = None,
        pagination: Pagination = None,
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
    es_query &= _owner_es_query(owner=owner, user_id=user_id)

    # pagination
    if pagination is None:
        pagination = Pagination()

    search = Search(index=config.elastic.index_name)

    search = search.query(es_query)
    order_field = entry_type.quantities[pagination.order_by].search_field
    sort = {order_field: pagination.order.value}
    if order_field != 'calc_id':
        sort['calc_id'] = pagination.order.value
    search = search.sort(sort)
    search = search.extra(size=pagination.size)
    if pagination.after:
        search = search.extra(search_after=pagination.after.rsplit(':', 1))

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
        raise SearchError(e)

    more_response_data = {}

    # pagination
    next_after = None
    if 0 < len(es_response.hits) < es_response.hits.total:
        last = es_response.hits[-1]
        if order_field == 'calc_id':
            next_after = last['calc_id']
        else:
            after_value = last
            for order_field_segment in order_field.split('.'):
                after_value = after_value[order_field_segment]
            next_after = '%s:%s' % (after_value, last['calc_id'])
    pagination_response = PaginationResponse(
        total=es_response.hits.total,
        next_after=next_after,
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
        raise SearchError(e)

    return result
