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

from typing import cast, Any, Dict, Union, List
from elasticsearch.exceptions import RequestError
from elasticsearch_dsl import Search, A, Q
from elasticsearch_dsl.query import Query as EsQuery
from pydantic.error_wrappers import ErrorWrapper

from nomad.metainfo.elasticsearch_extension import (
    material_type, entry_type, material_entry_type,
    entry_index, Index, index_entries, DocumentType, SearchQuantity)
from nomad.app.v1 import models as api_models
from nomad.app.v1.models import (
    Pagination, PaginationResponse, Query, MetadataRequired, MetadataResponse, Aggregation,
    Statistic, StatisticResponse, AggregationOrderType, AggregationResponse, AggregationDataItem,
    Value)

from .common import SearchError, _es_to_entry_dict, _owner_es_query


class QueryValidationError(Exception):
    def __init__(self, error, loc):
        self.errors = [ErrorWrapper(Exception(error), loc=loc)]


def validate_quantity(
        quantity_name: str, value: Value = None, doc_type: DocumentType = None,
        loc: List[str] = None) -> SearchQuantity:
    '''
    Validates the given quantity name and value against the given document type.

    Returns:
        A metainfo elasticsearch extension SearchQuantity object.

    Raises: QueryValidationError
    '''
    assert quantity_name is not None

    if doc_type == material_entry_type and not quantity_name.startswith('entries'):
        quantity_name = f'entries.{quantity_name}'

    if doc_type == material_type and quantity_name.startswith('entries'):
        doc_type = material_entry_type

    if doc_type is None:
        doc_type = entry_type

    quantity = doc_type.quantities.get(quantity_name)
    if quantity is None:
        raise QueryValidationError(
            f'{quantity_name} is not a {doc_type} quantity',
            loc=[quantity_name] if loc is None else loc)

    return quantity


def validate_api_query(
        query: Query, doc_type: DocumentType, owner_query: EsQuery,
        prefix: str = None) -> EsQuery:
    '''
    Creates an ES query based on the API's query model. This needs to be a normalized
    query expression with explicit objects for logical, set, and comparison operators.
    Shorthand notations ala ``quantity:operator`` are not supported here; this
    needs to be resolved via the respective pydantic validator.

    However, this function performs validation of quantities and types and raises
    a QueryValidationError accordingly. This exception is populated with pydantic
    errors.

    Arguments:
        query: The api query object.
        doc_type:
            The elasticsearch metainfo extension document type that this query needs to
            be verified against.
        owner_query:
            A prebuild ES query that is added to nested entries query. Only for
            materials queries.
        prefix:
            An optional prefix that is added to all quantity names. Used for recursion.

    Returns:
        A elasticsearch dsl query object.

    Raises: QueryValidationError
    '''

    def match(name: str, value: Value) -> EsQuery:
        # TODO non keyword quantities, quantities with value transformation, type checks
        quantity = validate_quantity(name, value, doc_type=doc_type)
        return Q('match', **{quantity.search_field: value})

    def validate_query(query: Query) -> EsQuery:
        return validate_api_query(
            query, doc_type=doc_type, owner_query=owner_query, prefix=prefix)

    def validate_criteria(name: str, value: Any):
        if prefix is not None:
            name = f'{prefix}.{name}'

        # handle prefix and nested queries
        for nested_key in doc_type.nested_object_keys:
            if len(name) < len(nested_key):
                break
            if not name.startswith(nested_key):
                continue
            if prefix is not None and prefix.startswith(nested_key):
                continue
            if nested_key == name and isinstance(value, api_models.Nested):
                continue

            value = api_models.Nested(query={name[len(nested_key) + 1:]: value})
            name = nested_key
            break

        if isinstance(value, api_models.All):
            return Q('bool', must=[match(name, item) for item in value.op])

        elif isinstance(value, api_models.Any_):
            return Q('bool', should=[match(name, item)for item in value.op])

        elif isinstance(value, api_models.None_):
            return Q('bool', must_not=[match(name, item) for item in value.op])

        elif isinstance(value, api_models.ComparisonOperator):
            # TODO typecheck?
            quantity = validate_quantity(name, None, doc_type=doc_type)
            field = quantity.search_field
            return Q('range', **{field: {type(value).__name__.lower(): value.op}})

        elif isinstance(value, (api_models.And, api_models.Or, api_models.Not)):
            return validate_query(value)

        elif isinstance(value, api_models.Nested):
            sub_doc_type = material_entry_type if name == 'entries' else doc_type

            sub_query = validate_api_query(
                value.query, doc_type=sub_doc_type, prefix=name, owner_query=owner_query)

            if name in doc_type.nested_object_keys:
                if name == 'entries':
                    sub_query &= owner_query
                return Q('nested', path=name, query=sub_query)
            else:
                return sub_query

        # list of values is treated as an "all" over the items
        elif isinstance(value, list):
            return Q('bool', must=[match(name, item) for item in value])

        elif isinstance(value, dict):
            assert False, (
                'Using dictionaries as criteria values directly is not supported. Use the '
                'Nested model.')

        else:
            return match(name, value)

    if isinstance(query, api_models.And):
        return Q('bool', must=[validate_query(operand) for operand in query.op])

    if isinstance(query, api_models.Or):
        return Q('bool', should=[validate_query(operand) for operand in query.op])

    if isinstance(query, api_models.Not):
        return Q('bool', must_not=validate_query(query.op))

    if isinstance(query, dict):
        # dictionary is like an "and" of all items in the dict
        if len(query) == 0:
            return Q()

        if len(query) == 1:
            key = next(iter(query))
            return validate_criteria(key, query[key])

        return Q('bool', must=[
            validate_criteria(name, value) for name, value in query.items()])

    raise NotImplementedError()


def validate_pagination(pagination: Pagination, doc_type: DocumentType, loc: List[str] = None):
    order_quantity = None
    if pagination.order_by is not None:
        order_quantity = validate_quantity(
            pagination.order_by, doc_type=doc_type, loc=['pagination', 'order_by'])
        if not order_quantity.definition.is_scalar:
            raise QueryValidationError(
                'the order_by quantity must be a scalar',
                loc=(loc if loc else []) + ['pagination', 'order_by'])

    page_after_value = pagination.page_after_value
    if page_after_value is not None and \
            pagination.order_by is not None and \
            pagination.order_by != doc_type.id_field and \
            ':' not in page_after_value:

        pagination.page_after_value = '%s:' % page_after_value

    return order_quantity, page_after_value


def _api_to_es_statistic(
        es_search: Search, name: str, statistic: Statistic, doc_type: DocumentType) -> A:
    '''
    Creates an ES aggregation based on the API's statistic model.
    '''

    quantity = validate_quantity(statistic.quantity, loc=['statistic', 'quantity'], doc_type=doc_type)

    if not quantity.aggregateable:
        raise QueryValidationError(
            'the statistic quantity cannot be aggregated',
            loc=['statistic', 'quantity'])

    if statistic.size is None:
        statistic.size = quantity.statistics_size

    if quantity.values is not None:
        statistic.size = len(quantity.values)

    terms_kwargs = {}
    if statistic.value_filter is not None:
        terms_kwargs['include'] = '.*%s.*' % statistic.value_filter

    aggs = es_search.aggs
    for nested_key in doc_type.nested_object_keys:
        if statistic.quantity.startswith(nested_key):
            aggs = es_search.aggs.bucket('nested_statistic:%s' % name, 'nested', path=nested_key)

    order_type = '_count' if statistic.order.type_ == AggregationOrderType.entries else '_key'
    statistic_agg = aggs.bucket('statistic:%s' % name, A(
        'terms',
        field=quantity.search_field,
        size=statistic.size,
        order={order_type: statistic.order.direction.value},
        **terms_kwargs))

    for metric_name in statistic.metrics:
        metrics = doc_type.metrics
        if nested_key == 'entries':
            metrics = material_entry_type.metrics
        if metric_name not in metrics:
            raise QueryValidationError(
                'metric must be the qualified name of a suitable search quantity',
                loc=['statistic', 'metrics'])
        metric_aggregation, metric_quantity = metrics[metric_name]
        statistic_agg.metric('metric:%s' % metric_name, A(
            metric_aggregation,
            field=metric_quantity.qualified_field))


def _es_to_api_statistics(
        es_response, name: str, statistic: Statistic, doc_type: DocumentType) -> StatisticResponse:
    '''
    Creates a StatisticResponse from elasticsearch response on a request executed with
    the given statistics.
    '''
    quantity = validate_quantity(statistic.quantity, doc_type=doc_type)

    es_aggs = es_response.aggs
    for nested_key in doc_type.nested_object_keys:
        if statistic.quantity.startswith(nested_key):
            es_aggs = es_response.aggs[f'nested_statistic:{name}']

    es_statistic = es_aggs['statistic:' + name]
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


def _api_to_es_aggregation(
        es_search: Search, name: str, agg: Aggregation, doc_type: DocumentType) -> A:
    '''
    Creates an ES aggregation based on the API's aggregation model.
    '''
    order_quantity, page_after_value = validate_pagination(
        agg.pagination, doc_type=doc_type, loc=['aggration'])

    quantity = validate_quantity(agg.quantity, doc_type=doc_type, loc=['aggregation', 'quantity'])
    if not quantity.aggregateable:
        raise QueryValidationError(
            'the aggregation quantity cannot be aggregated',

            loc=['aggregation', 'quantity'])

    terms = A('terms', field=quantity.search_field, order=agg.pagination.order.value)

    # We are using elastic searchs 'composite aggregations' here. We do not really
    # compose aggregations, but only those pseudo composites allow us to use the
    # 'after' feature that allows to scan through all aggregation values.
    if order_quantity is None:
        composite = {
            'sources': {
                name: terms
            },
            'size': agg.pagination.page_size
        }
    else:
        sort_terms = A('terms', field=order_quantity.search_field, order=agg.pagination.order.value)
        composite = {
            'sources': [
                {order_quantity.search_field: sort_terms},
                {quantity.search_field: terms}
            ],
            'size': agg.pagination.page_size
        }

    if page_after_value is not None:
        if order_quantity is None:
            composite['after'] = {name: page_after_value}
        else:
            order_value, quantity_value = page_after_value.split(':')
            composite['after'] = {quantity.search_field: quantity_value, order_quantity.search_field: order_value}

    aggs = es_search.aggs
    for nested_key in doc_type.nested_object_keys:
        if agg.quantity.startswith(nested_key):
            aggs = es_search.aggs.bucket('nested_agg:%s' % name, 'nested', path=nested_key)

    composite_agg = aggs.bucket('agg:%s' % name, 'composite', **composite)

    if agg.entries is not None and agg.entries.size > 0:
        kwargs: Dict[str, Any] = {}
        if agg.entries.required is not None:
            if agg.entries.required.include is not None:
                kwargs.update(_source=dict(includes=agg.entries.required.include))
            else:
                kwargs.update(_source=dict(excludes=agg.entries.required.exclude))

        composite_agg.metric('entries', A('top_hits', size=agg.entries.size, **kwargs))

    # additional cardinality to get total
    aggs.metric('agg:%s:total' % name, 'cardinality', field=quantity.search_field)


def _es_to_api_aggregation(
        es_response, name: str, agg: Aggregation, doc_type: DocumentType) -> AggregationResponse:
    '''
    Creates a AggregationResponse from elasticsearch response on a request executed with
    the given aggregation.
    '''
    order_by = agg.pagination.order_by
    quantity = validate_quantity(agg.quantity, doc_type=doc_type)

    nested = False
    es_aggs = es_response.aggs
    for nested_key in doc_type.nested_object_keys:
        if agg.quantity.startswith(nested_key):
            es_aggs = es_response.aggs[f'nested_agg:{name}']
            nested = True

    es_agg = es_aggs['agg:' + name]

    def get_entries(agg):
        if 'entries' in agg:
            if nested:
                return [{nested_key: item['_source']} for item in agg.entries.hits.hits]
            else:
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
        total=es_aggs['agg:%s:total' % name]['value'],
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
        query: Union[Query, EsQuery] = None,
        pagination: Pagination = None,
        required: MetadataRequired = None,
        aggregations: Dict[str, Aggregation] = {},
        statistics: Dict[str, Statistic] = {},
        user_id: str = None,
        index: Index = entry_index) -> MetadataResponse:

    # The first half of this method creates the ES query. Then the query is run on ES.
    # The second half is about transforming the ES response to a MetadataResponse.

    doc_type = index.doc_type

    # owner and query
    owner_query = _owner_es_query(owner=owner, user_id=user_id, doc_type=doc_type)

    if query is None:
        query = {}

    if isinstance(query, EsQuery):
        es_query = cast(EsQuery, query)
    else:
        es_query = validate_api_query(
            cast(Query, query), doc_type=doc_type, owner_query=owner_query)

    if doc_type != entry_type:
        es_query &= Q('nested', path='entries', query=owner_query)
    else:
        es_query &= owner_query

    # pagination
    if pagination is None:
        pagination = Pagination()

    if pagination.order_by is None:
        pagination.order_by = doc_type.id_field

    search = Search(index=index.index_name)

    search = search.query(es_query)
    # TODO this depends on doc_type
    if pagination.order_by is None:
        pagination.order_by = doc_type.id_field
    order_quantity, page_after_value = validate_pagination(pagination, doc_type=doc_type)
    order_field = order_quantity.search_field
    sort = {order_field: pagination.order.value}
    if order_field != doc_type.id_field:
        sort[doc_type.id_field] = pagination.order.value
    search = search.sort(sort)
    search = search.extra(size=pagination.page_size)
    if page_after_value:
        search = search.extra(search_after=page_after_value.rsplit(':', 1))

    # required
    if required:
        for list_ in [required.include, required.exclude]:
            for quantity in [] if list_ is None else list_:
                # TODO validate quantities with wildcards
                if '*' not in quantity:
                    validate_quantity(quantity, doc_type=doc_type, loc=['required'])

        if required.include is not None and pagination.order_by not in required.include:
            required.include.append(pagination.order_by)
        if required.exclude is not None and pagination.order_by in required.exclude:
            required.exclude.remove(pagination.order_by)

        if required.include is not None and doc_type.id_field not in required.include:
            required.include.append(doc_type.id_field)

        if required.exclude is not None and doc_type.id_field in required.exclude:
            required.exclude.remove(doc_type.id_field)

        search = search.source(includes=required.include, excludes=required.exclude)

    # statistics
    for name, statistic in statistics.items():
        _api_to_es_statistic(search, name, statistic, doc_type=doc_type)

    # aggregations
    for name, agg in aggregations.items():
        _api_to_es_aggregation(search, name, agg, doc_type=doc_type)

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
        if order_field == doc_type.id_field:
            next_page_after_value = last[doc_type.id_field]
        else:
            after_value = last
            for order_field_segment in order_field.split('.'):
                after_value = after_value[order_field_segment]
            next_page_after_value = '%s:%s' % (after_value, last[doc_type.id_field])
    pagination_response = PaginationResponse(
        total=es_response.hits.total,
        next_page_after_value=next_page_after_value,
        **pagination.dict())

    # statistics
    if len(statistics) > 0:
        more_response_data['statistics'] = cast(Dict[str, Any], {
            name: _es_to_api_statistics(es_response, name, statistic, doc_type=doc_type)
            for name, statistic in statistics.items()})

    # aggregations
    if len(aggregations) > 0:
        more_response_data['aggregations'] = cast(Dict[str, Any], {
            name: _es_to_api_aggregation(es_response, name, aggregation, doc_type=doc_type)
            for name, aggregation in aggregations.items()})

    more_response_data['es_query'] = es_query.to_dict()

    result = MetadataResponse(
        owner='all' if owner is None else owner,
        query=query,
        pagination=pagination_response,
        required=required,
        data=[_es_to_entry_dict(hit, required) for hit in es_response.hits],
        **more_response_data)

    return result


def _index(entries, **kwargs):
    index_entries(entries, **kwargs)
