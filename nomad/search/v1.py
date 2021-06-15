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
    Value, AggregationBase, TermsAggregation, BucketAggregation, HistogramAggregation,
    DateHistogramAggregation, MinMaxAggregation, Bucket,
    MixMaxAggregationResponse, TermsAggregationResponse, HistogramAggregationResponse,
    DateHistogramAggregationResponse, AggregationResponse)

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

        elif isinstance(value, api_models.RangeNumber):
            quantity = validate_quantity(name, None, doc_type=doc_type)
            field = quantity.search_field
            return Q('range', **{field: value.dict(
                include={'lt', 'lte', 'gt', 'gte'},
                exclude_unset=True,
            )})

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


def _api_to_es_aggregation(
        es_search: Search, name: str, agg: AggregationBase, doc_type: DocumentType) -> A:
    '''
    Creates an ES aggregation based on the API's aggregation model.
    '''

    agg_name = f'agg:{name}'
    quantity = validate_quantity(agg.quantity, doc_type=doc_type, loc=['aggregation', 'quantity'])
    es_aggs = es_search.aggs
    for nested_key in doc_type.nested_object_keys:
        if agg.quantity.startswith(nested_key):
            es_aggs = es_search.aggs.bucket('nested_agg:%s' % name, 'nested', path=nested_key)

    es_agg = None
    if isinstance(agg, TermsAggregation):
        if not quantity.aggregateable:
            raise QueryValidationError(
                'The aggregation quantity cannot be used in a terms aggregation.',
                loc=['aggregation', name, 'terms', 'quantity'])

        if agg.pagination is not None:
            if agg.size is not None:
                raise QueryValidationError(
                    f'You cannot paginate and provide an extra size parameter.',
                    loc=['aggregations', name, 'terms', 'pagination'])

            order_quantity, page_after_value = validate_pagination(
                agg.pagination, doc_type=doc_type, loc=['aggregation'])

            # We are using elastic searchs 'composite aggregations' here. We do not really
            # compose aggregations, but only those pseudo composites allow us to use the
            # 'after' feature that allows to scan through all aggregation values.
            terms = A('terms', field=quantity.search_field, order=agg.pagination.order.value)

            if order_quantity is None:
                composite = {
                    'sources': {
                        name: terms
                    },
                    'size': agg.pagination.page_size
                }

            else:
                sort_terms = A(
                    'terms',
                    field=order_quantity.search_field,
                    order=agg.pagination.order.value)

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
                    try:
                        order_value, quantity_value = page_after_value.split(':')
                        composite['after'] = {quantity.search_field: quantity_value, order_quantity.search_field: order_value}
                    except Exception:
                        raise QueryValidationError(
                            f'The pager_after_value has not the right format.',
                            loc=['aggregations', name, 'terms', 'pagination', 'page_after_value'])

            es_agg = es_aggs.bucket(agg_name, 'composite', **composite)

            # additional cardinality to get total
            es_aggs.metric('agg:%s:total' % name, 'cardinality', field=quantity.search_field)
        else:
            if agg.size is None:
                if quantity.default_aggregation_size is not None:
                    agg.size = quantity.default_aggregation_size

                elif quantity.values is not None:
                    agg.size = len(quantity.values)

                else:
                    agg.size = 10

            terms_kwargs = {}
            if agg.value_filter is not None:
                terms_kwargs['include'] = '.*%s.*' % agg.value_filter

            terms = A('terms', field=quantity.search_field, size=agg.size, **terms_kwargs)
            es_agg = es_aggs.bucket(agg_name, terms)

        if agg.entries is not None and agg.entries.size > 0:
            kwargs: Dict[str, Any] = {}
            if agg.entries.required is not None:
                if agg.entries.required.include is not None:
                    kwargs.update(_source=dict(includes=agg.entries.required.include))
                else:
                    kwargs.update(_source=dict(excludes=agg.entries.required.exclude))

            es_agg.metric('entries', A('top_hits', size=agg.entries.size, **kwargs))

    elif isinstance(agg, DateHistogramAggregation):
        if not quantity.annotation.mapping['type'] in ['date']:
            raise QueryValidationError(
                f'The quantity {quantity} cannot be used in a date histogram aggregation',
                loc=['aggregations', name, 'histogram', 'quantity'])

        es_agg = es_aggs.bucket(agg_name, A(
            'date_histogram', field=quantity.search_field, interval=agg.interval,
            format='yyyy-MM-dd'))

    elif isinstance(agg, HistogramAggregation):
        if not quantity.annotation.mapping['type'] in ['integer', 'float', 'double', 'long']:
            raise QueryValidationError(
                f'The quantity {quantity} cannot be used in a histogram aggregation',
                loc=['aggregations', name, 'histogram', 'quantity'])

        es_agg = es_aggs.bucket(agg_name, A(
            'histogram', field=quantity.search_field, interval=agg.interval))

    elif isinstance(agg, MinMaxAggregation):
        if not quantity.annotation.mapping['type'] in ['integer', 'float', 'double', 'long']:
            raise QueryValidationError(
                f'The quantity {quantity} cannot be used in a mix-max aggregation',
                loc=['aggregations', name, 'min_max', 'quantity'])

        es_aggs.metric(agg_name + ':min', A('min', field=quantity.search_field))
        es_aggs.metric(agg_name + ':max', A('max', field=quantity.search_field))

    else:
        raise NotImplementedError()

    if isinstance(agg, BucketAggregation):
        for metric_name in agg.metrics:
            metrics = doc_type.metrics
            if nested_key == 'entries':
                metrics = material_entry_type.metrics
            if metric_name not in metrics:
                raise QueryValidationError(
                    'metric must be the qualified name of a suitable search quantity',
                    loc=['statistic', 'metrics'])
            metric_aggregation, metric_quantity = metrics[metric_name]
            es_agg.metric('metric:%s' % metric_name, A(
                metric_aggregation,
                field=metric_quantity.qualified_field))


def _es_to_api_aggregation(
        es_response, name: str, agg: AggregationBase, doc_type: DocumentType):
    '''
    Creates a AggregationResponse from elasticsearch response on a request executed with
    the given aggregation.
    '''
    quantity = validate_quantity(agg.quantity, doc_type=doc_type)

    nested = False
    es_aggs = es_response.aggs
    for nested_key in doc_type.nested_object_keys:
        if agg.quantity.startswith(nested_key):
            es_aggs = es_response.aggs[f'nested_agg:{name}']
            nested = True

    aggregation_dict = agg.dict(by_alias=True)
    has_no_pagination = getattr(agg, 'pagination', None) is None

    if isinstance(agg, BucketAggregation):
        es_agg = es_aggs['agg:' + name]
        values = set()

        def get_bucket(es_bucket) -> Bucket:
            if has_no_pagination:
                if isinstance(agg, DateHistogramAggregation):
                    value = es_bucket['key_as_string']
                else:
                    value = es_bucket['key']
            elif agg.pagination.order_by is None:  # type: ignore
                value = es_bucket.key[name]
            else:
                value = es_bucket.key[quantity.search_field]

            count = es_bucket.doc_count
            metrics = {}
            for metric in agg.metrics:  # type: ignore
                metrics[metric] = es_bucket['metric:' + metric].value

            entries = None
            if 'entries' in es_bucket:
                if nested:
                    entries = [{nested_key: item['_source']} for item in es_bucket.entries.hits.hits]
                else:
                    entries = [item['_source'] for item in es_bucket.entries.hits.hits]

            values.add(value)
            if len(metrics) == 0:
                metrics = None
            return Bucket(value=value, entries=entries, count=count, metrics=metrics)

        data = [get_bucket(es_bucket) for es_bucket in es_agg.buckets]

        if has_no_pagination:
            # fill "empty" values
            if quantity.values is not None:
                for value in quantity.values:
                    if value not in values:
                        metrics = {metric: 0 for metric in agg.metrics}
                        if len(metrics) == 0:
                            metrics = None
                        data.append(Bucket(value=value, count=0, metrics=metrics))

        else:
            total = es_aggs['agg:%s:total' % name]['value']
            pagination = PaginationResponse(total=total, **aggregation_dict['pagination'])
            if pagination.page_after_value is not None and pagination.page_after_value.endswith(':'):
                pagination.page_after_value = pagination.page_after_value[0:-1]

            if 'after_key' in es_agg:
                after_key = es_agg['after_key']
                if pagination.order_by is None:
                    pagination.next_page_after_value = after_key[name]
                else:
                    str_values = [str(v) for v in after_key.to_dict().values()]
                    pagination.next_page_after_value = ':'.join(str_values)
            else:
                pagination.next_page_after_value = None

            aggregation_dict['pagination'] = pagination

        if isinstance(agg, TermsAggregation):
            return AggregationResponse(
                terms=TermsAggregationResponse(data=data, **aggregation_dict))
        elif isinstance(agg, HistogramAggregation):
            return AggregationResponse(
                histogram=HistogramAggregationResponse(data=data, **aggregation_dict))
        elif isinstance(agg, DateHistogramAggregation):
            return AggregationResponse(
                date_histogram=DateHistogramAggregationResponse(data=data, **aggregation_dict))
        else:
            raise NotImplementedError()

    elif isinstance(agg, MinMaxAggregation):
        min_value = es_aggs['agg:%s:min' % name]['value']
        max_value = es_aggs['agg:%s:max' % name]['value']

        return AggregationResponse(
            min_max=MixMaxAggregationResponse(data=[min_value, max_value], **aggregation_dict))

    else:
        raise NotImplementedError()


def _specific_agg(agg: Aggregation) -> Union[TermsAggregation, DateHistogramAggregation, HistogramAggregation, MinMaxAggregation]:
    if agg.terms is not None:
        return agg.terms

    if agg.histogram is not None:
        return agg.histogram

    if agg.date_histogram is not None:
        return agg.date_histogram

    if agg.min_max is not None:
        return agg.min_max

    raise NotImplementedError()


def search(
        owner: str = 'public',
        query: Union[Query, EsQuery] = None,
        pagination: Pagination = None,
        required: MetadataRequired = None,
        aggregations: Dict[str, Aggregation] = {},
        user_id: str = None,
        index: Index = entry_index) -> MetadataResponse:

    # The first half of this method creates the ES query. Then the query is run on ES.
    # The second half is about transforming the ES response to a MetadataResponse.

    print(query)
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

    print(es_query)

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

    # aggregations
    for name, agg in aggregations.items():
        _api_to_es_aggregation(search, name, _specific_agg(agg), doc_type=doc_type)

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

    # aggregations
    if len(aggregations) > 0:
        more_response_data['aggregations'] = cast(Dict[str, Any], {
            name: _es_to_api_aggregation(es_response, name, _specific_agg(agg), doc_type=doc_type)
            for name, agg in aggregations.items()})

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
