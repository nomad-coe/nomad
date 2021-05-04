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

from typing import cast, Any, Dict, Union
from elasticsearch.exceptions import RequestError
from elasticsearch_dsl import Search, A

from nomad.metainfo.elasticsearch_extension import entry_type, entry_index, Index, index_entries
from nomad.app.v1.models import (
    Pagination, PaginationResponse, Query, MetadataRequired, SearchResponse, Aggregation,
    Statistic, StatisticResponse, AggregationOrderType, AggregationResponse, AggregationDataItem)

from .common import SearchError, _api_to_es_query, _es_to_entry_dict, _owner_es_query


def _api_to_es_statistic(es_search: Search, name: str, statistic: Statistic) -> A:
    '''
    Creates an ES aggregation based on the API's statistic model.
    '''

    quantity = entry_type.quantities[statistic.quantity]
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

    for metric_name in statistic.metrics:
        metric_aggregation, metric_quantity = entry_type.metrics[metric_name]
        statistic_agg.metric('metric:%s' % metric_name, A(
            metric_aggregation,
            field=metric_quantity.qualified_field))


def _es_to_api_statistics(es_response, name: str, statistic: Statistic) -> StatisticResponse:
    '''
    Creates a StatisticResponse from elasticsearch response on a request executed with
    the given statistics.
    '''
    quantity = entry_type.quantities[statistic.quantity]

    es_statistic = es_response.aggs['statistic:' + name]
    statistic_data = {}
    for bucket in es_statistic.buckets:
        value_data = dict(entries=bucket.doc_count)
        for metric in statistic.metrics:
            value_data[metric] = bucket['metric:' + metric].value
        statistic_data[bucket.key] = value_data

    if quantity.values is not None:
        for value in quantity.values:
            if value not in statistic_data:
                statistic_data[value] = dict(entries=0, **{
                    metric: 0 for metric in statistic.metrics})

    return StatisticResponse(data=statistic_data, **statistic.dict(by_alias=True))


def _api_to_es_aggregation(es_search: Search, name: str, agg: Aggregation) -> A:
    '''
    Creates an ES aggregation based on the API's aggregation model.
    '''
    quantity = entry_type.quantities[agg.quantity]
    terms = A('terms', field=quantity.search_field, order=agg.pagination.order.value)

    # We are using elastic searchs 'composite aggregations' here. We do not really
    # compose aggregations, but only those pseudo composites allow us to use the
    # 'after' feature that allows to scan through all aggregation values.
    order_by = agg.pagination.order_by
    if order_by is None:
        composite = dict(sources={name: terms}, size=agg.pagination.page_size)
    else:
        order_quantity = entry_type.quantities[order_by]
        sort_terms = A('terms', field=order_quantity.search_field, order=agg.pagination.order.value)
        composite = {
            'sources': [
                {order_quantity.search_field: sort_terms},
                {quantity.search_field: terms}
            ],
            'size': agg.pagination.page_size
        }

    if agg.pagination.page_after_value is not None:
        if order_by is None:
            composite['after'] = {name: agg.pagination.page_after_value}
        else:
            order_value, quantity_value = agg.pagination.page_after_value.split(':')
            composite['after'] = {quantity.search_field: quantity_value, order_quantity.search_field: order_value}

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
    quantity = entry_type.quantities[agg.quantity]
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
        pagination: Pagination = None,
        required: MetadataRequired = None,
        aggregations: Dict[str, Aggregation] = {},
        statistics: Dict[str, Statistic] = {},
        user_id: str = None,
        index: Union[Index, str] = entry_index) -> SearchResponse:

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

    if pagination.order_by is None:
        pagination.order_by = 'entry_id'

    if isinstance(index, Index):
        index = index.index_name
    search = Search(index=index)

    search = search.query(es_query)
    order_field = entry_type.quantities[pagination.order_by].search_field
    sort = {order_field: pagination.order.value}
    if order_field != 'entry_id':
        sort['entry_id'] = pagination.order.value
    search = search.sort(sort)
    search = search.extra(size=pagination.page_size)
    if pagination.page_after_value:
        search = search.extra(search_after=pagination.page_after_value.rsplit(':', 1))

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
    next_page_after_value = None
    if 0 < len(es_response.hits) < es_response.hits.total:
        last = es_response.hits[-1]
        if order_field == 'entry_id':
            next_page_after_value = last['entry_id']
        else:
            after_value = last
            for order_field_segment in order_field.split('.'):
                after_value = after_value[order_field_segment]
            next_page_after_value = '%s:%s' % (after_value, last['entry_id'])
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

    result = SearchResponse(
        owner='all' if owner is None else owner,
        query=query,
        pagination=pagination_response,
        required=required,
        data=[_es_to_entry_dict(hit, required) for hit in es_response.hits],
        **more_response_data)

    return result


def _index(entries, **kwargs):
    index_entries(entries, **kwargs)
