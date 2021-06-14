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

from typing import List, Dict, Optional, Union, Any, Mapping
import enum
from fastapi import Body, Request, HTTPException, Query as FastApiQuery
import pydantic
from pydantic import BaseModel, Field, validator, root_validator
import datetime
import numpy as np
import re
import fnmatch

from nomad import datamodel  # pylint: disable=unused-import
from nomad.utils import strip
from nomad.metainfo import Datetime, MEnum
from nomad.metainfo.search_extension import metrics, search_quantities, search_sub_sections

from .utils import parameter_dependency_from_model, update_url_query_arguments


User = datamodel.User.m_def.a_pydantic.model


calc_id = 'calc_id'
Metric = enum.Enum('Metric', {name: name for name in metrics})  # type: ignore
AggregateableQuantity = enum.Enum('AggregateableQuantity', {  # type: ignore
    name: name for name in search_quantities
    if search_quantities[name].aggregateable})

AggregateableQuantity.__doc__ = '''
    Statistics and aggregations can only be computed for those search quantities that have
    discrete values. For example a statistics aggregates a certain metric (e.g. the number of entries)
    over all entries were this quantity has the same value (bucket aggregation, think historgam here).
'''

Value = Union[str, int, float, bool, datetime.datetime]
ComparableValue = Union[str, int, float, datetime.datetime]


class AggregationOrderType(str, enum.Enum):
    '''
    Allows to order statistics or aggregations by either quantity values (`values`) or number
    of entries (`entries`).
    '''
    values = 'values'
    entries = 'entries'


class HTTPExceptionModel(BaseModel):
    detail: str


class NoneEmptyBaseModel(BaseModel):
    @root_validator
    def check_exists(cls, values):  # pylint: disable=no-self-argument
        assert any(value is not None for value in values.values())
        return values


class All(NoneEmptyBaseModel):
    op: List[Value] = Field(None, alias='all')


class None_(NoneEmptyBaseModel):
    op: List[Value] = Field(None, alias='none')


class Any_(NoneEmptyBaseModel):
    op: List[Value] = Field(None, alias='any')


class ComparisonOperator(NoneEmptyBaseModel): pass


class Lte(ComparisonOperator):
    op: ComparableValue = Field(None, alias='lte')


class Lt(ComparisonOperator):
    op: ComparableValue = Field(None, alias='lt')


class Gte(ComparisonOperator):
    op: ComparableValue = Field(None, alias='gte')


class Gt(ComparisonOperator):
    op: ComparableValue = Field(None, alias='gt')


class LogicalOperator(NoneEmptyBaseModel):

    @validator('op', check_fields=False)
    def validate_query(cls, query):  # pylint: disable=no-self-argument
        if isinstance(query, list):
            return [_validate_query(item) for item in query]

        return _validate_query(query)


class And(LogicalOperator):
    op: List['Query'] = Field(None, alias='and')


class Or(LogicalOperator):
    op: List['Query'] = Field(None, alias='or')


class Not(LogicalOperator):
    op: 'Query' = Field(None, alias='not')


ops = {
    'lte': Lte,
    'lt': Lt,
    'gte': Gte,
    'gt': Gt,
    'all': All,
    'none': None_,
    'any': Any_
}


QueryParameterValue = Union[Value, List[Value], Lte, Lt, Gte, Gt, Any_, All, None_]

Query = Union[
    Mapping[str, QueryParameterValue], And, Or, Not]


And.update_forward_refs()
Or.update_forward_refs()
Not.update_forward_refs()


class Owner(str, enum.Enum):
    '''
    The `owner` allows to limit the scope of the searched based on entry ownership.
    This is useful, if you only want to search among all publically downloadable
    entries, or only among your own entries, etc.

    These are the possible owner values and their meaning:
    * `all`: Consider all entries.
    * `public` (default): Consider all entries that can be publically downloaded,
        i.e. only published entries without embargo
    * `user`: Only consider entries that belong to you.
    * `shared`: Only consider entries that belong to you or are shared with you.
    * `visible`: Consider all entries that are visible to you. This includes
        entries with embargo or unpublished entries that belong to you or are
        shared with you.
    * `staging`: Only search through unpublished entries.
    '''

    # There seems to be a slight bug in fast API. When it creates the example in OpenAPI
    # it will ignore any given default or example and simply take the first enum value.
    # Therefore, we put public first, which is the most default and save in most contexts.
    public = 'public'
    all_ = 'all'
    visible = 'visible'
    shared = 'shared'
    user = 'user'
    staging = 'staging'
    admin = 'admin'


class WithQuery(BaseModel):
    owner: Optional[Owner] = Body('public')
    query: Optional[Query] = Body(
        None,
        embed=True,
        description=strip('''
            A query can be very simple list of parameters. Different parameters are combined
            with a logical **and**, values of the same parameter with also with a logical **and**.
            The following would search for all entries that are VASP calculations,
            contain *Na* **and** *Cl*, **and** are authored by *Stefano Curtarolo*
            **and** *Chris Wolverton*.
            ```
            {
                "atoms": ["Na", "Cl"],
                "dft.code_name": "VASP",
                "authors": ["Stefano Curtarolo", "Chris Wolverton"]
            }
            ```

            A short cut to change the logical combination of values in a list, is to
            add a suffix to the quantity `:any`:
            ```
            {
                "atoms": ["Na", "Cl"],
                "dft.code_name": "VASP",
                "authors:any": ["Stefano Curtarolo", "Chris Wolverton"]
            }
            ```

            Otherwise, you can also write complex logical combinations of parameters like this:
            ```
            {
                "and": [
                    {
                        "or": [
                            {
                                "atoms": ["Cl", "Na"]
                            },
                            {
                                "atoms": ["H", "O"]
                            }
                        ]
                    },
                    {
                        "not": {
                            "dft.crystal": "cubic"
                        }
                    }
                ]
            }
            ```
            Other short-cut prefixes are `none:` and `any:` (the default).

            By default all quantity values have to **equal** the given values to match. For
            some values you can also use comparison operators like this:
            ```
            {
                "upload_time": {
                    "gt": "2020-01-01",
                    "lt": "2020-08-01"
                },
                "dft.workflow.section_geometry_optimization.final_energy_difference": {
                    "lte": 1.23e-18
                }
            }
            ```

            or shorter with suffixes:
            ```
            {
                "upload_time:gt": "2020-01-01",
                "upload_time:lt": "2020-08-01",
                "dft.workflow.section_geometry_optimization.final_energy_difference:lte" 1.23e-18
            }
            ```

            The searchable quantities are a subset of the NOMAD Archive quantities defined
            in the NOMAD Metainfo. The most common quantities are: %s.
        ''' % ', '.join(reversed([
            '`%s`' % name
            for name in search_quantities
            if (name.startswith('dft') or '.' not in name) and len(name) < 20
        ]))),
        example={
            'upload_time:gt': '2020-01-01',
            'atoms': ['Ti', 'O'],
            'dft.code_name': 'VASP',
            'dft.workflow.section_geometry_optimization.final_energy_difference:lte': 1.23e-18,
            'dft.quantities': 'section_dos',
            'dft.system:any': ['bulk', '2d']
        })

    @validator('query')
    def validate_query(cls, query):  # pylint: disable=no-self-argument
        return _validate_query(query)


def _validate_query(query: Query):
    if isinstance(query, dict):
        for key, value in list(query.items()):
            # Note, we loop over a list of items, not query.items(). This is because we
            # may modify the query in the loop.
            if ':' in key:
                quantity, qualifier = key.split(':')
            else:
                quantity, qualifier = key, None

            assert quantity in search_quantities, '%s is not a searchable quantity' % key
            if qualifier is not None:
                assert quantity not in query, 'a quantity can only appear once in a query'
                assert qualifier in ops, 'unknown quantity qualifier %s' % qualifier
                del(query[key])
                query[quantity] = ops[qualifier](**{qualifier: value})  # type: ignore
            elif isinstance(value, list):
                query[quantity] = All(all=value)

    return query


def query_parameters(
    request: Request,
    owner: Optional[Owner] = FastApiQuery(
        'public', description=strip(Owner.__doc__)),
    q: Optional[List[str]] = FastApiQuery(
        [], description=strip('''
            Since we cannot properly offer forms for all parameters in the OpenAPI dashboard,
            you can use the parameter `q` and encode a query parameter like this
            `atoms__H` or `n_atoms__gt__3`. Multiple usage of `q` will combine parameters with
            logical *and*.
        '''))) -> WithQuery:

    # copy parameters from request
    query_params = {
        key: request.query_params.getlist(key)
        for key in request.query_params}

    # add the encoded parameters
    for parameter in q:
        fragments = parameter.split('__')
        if len(fragments) == 1 or len(fragments) > 3:
            raise HTTPException(422, detail=[{
                'loc': ['query', 'q'],
                'msg': 'wrong format, use <quantity>[__<op>]__<value>'}])
        name_op, value = '__'.join(fragments[:-1]), fragments[-1]
        quantity_name = name_op.split('__')[0]

        if quantity_name not in search_quantities:
            raise HTTPException(422, detail=[{
                'loc': ['query', parameter],
                'msg': '%s is not a search quantity' % quantity_name}])

        query_params.setdefault(name_op, []).append(value)

    # transform query parameters to query
    query: Dict[str, Any] = {}
    for key in query_params:
        op = None
        if '__' in key:
            quantity_name, op = key.split('__')
        else:
            quantity_name = key

        if quantity_name not in search_quantities:
            continue

        quantity = search_quantities[quantity_name]
        type_ = quantity.definition.type
        if type_ is Datetime:
            type_ = datetime.datetime.fromisoformat
        elif isinstance(type_, MEnum):
            type_ = str
        elif isinstance(type_, np.dtype):
            type_ = float
        elif type_ not in [int, float, bool]:
            type_ = str
        values = query_params[key]
        values = [type_(value) for value in values]

        if op is None:
            if quantity.many_and:
                op = 'all'
            if quantity.many_or:
                op = 'any'

        if op is None:
            if len(values) > 1:
                raise HTTPException(
                    status_code=422,
                    detail=[{
                        'loc': ['query', key],
                        'msg':'search parameter %s does not support multiple values' % key}])
            query[quantity_name] = values[0]

        elif op == 'all':
            query[quantity_name] = All(all=values)
        elif op == 'any':
            query[quantity_name] = Any_(any=values)
        elif op in ops:
            if len(values) > 1:
                raise HTTPException(
                    status_code=422,
                    detail=[{
                        'loc': ['query', key],
                        'msg': 'operator %s does not support multiple values' % op}])
            query[quantity_name] = ops[op](**{op: values[0]})
        else:
            raise HTTPException(
                422, detail=[{'loc': ['query', key], 'msg': 'operator %s is unknown' % op}])

    return WithQuery(query=query, owner=owner)


class Direction(str, enum.Enum):
    '''
    Order direction, either ascending (`asc`) or descending (`desc`)
    '''
    asc = 'asc'
    desc = 'desc'


class MetadataRequired(BaseModel):
    ''' Defines which metadata quantities are included or excluded in the response. '''

    include: Optional[List[str]] = Field(
        None, description=strip('''
            Quantities to include for each result. Only those quantities will be
            returned. The entry id quantity `calc_id` will always be included.
        '''))
    exclude: Optional[List[str]] = Field(
        None, description=strip('''
            Quantities to exclude for each result. Only all other quantities will
            be returned. The quantity `calc_id` cannot be excluded.
        '''))

    @validator('include', 'exclude')
    def validate_include(cls, value, values, field):  # pylint: disable=no-self-argument
        if value is None:
            return None

        for item in value:
            assert item in search_quantities or item in search_sub_sections or item[-1] == '*', \
                f'required fields ({item}) must be valid search quantities or contain wildcards'

        if field.name == 'include' and 'calc_id' not in value:
            value.append('calc_id')

        if field.name == 'exclude':
            if 'calc_id' in value:
                value.remove('calc_id')

        return value


metadata_required_parameters = parameter_dependency_from_model(
    'metadata_required_parameters', MetadataRequired)


class Pagination(BaseModel):
    ''' Defines the order, size, and page of results. '''

    page_size: Optional[int] = Field(
        10, description=strip('''
            The page size, e.g. the maximum number of items contained in one response.
            A `page_size` of 0 will return no results.
        '''))
    order_by: Optional[str] = Field(
        None,  # type: ignore
        description=strip('''
            The results are ordered by the values of this field. If omitted, default
            ordering is applied.
        '''))
    order: Optional[Direction] = Field(
        Direction.asc, description=strip('''
            The ordering direction of the results based on `order_by`. Its either
            ascending `asc` or decending `desc`. Default is `asc`.
        '''))
    page_after_value: Optional[str] = Field(
        None, description=strip('''
            This attribute defines the position after which the page begins, and is used
            to navigate through the total list of results.

            When requesting the first page, no value should be provided for
            `page_after_value`. Each response will contain a value `next_page_after_value`,
            which can be used to obtain the next page (by setting `page_after_value` in
            your next request to this value).

            The field is encoded as a string, and the format of `page_after_value` and
            `next_page_after_value` depends on which API method is used.

            Some API functions additionally allows a simplified navigation, by specifying
            the page number in the key `page`. It is however always possible to use
            `page_after_value` and `next_page_after_value` to iterate through the results.
            '''))
    page: Optional[int] = Field(
        None, description=strip('''
            The number of the page (1-based). When provided in a request, this attribute
            can be used instead of `page_after_value` to jump to a particular results page.

            **NOTE #1**: the option to request pages by submitting the `page` number is
            limited. There are api calls where this attribute cannot be used for indexing,
            or where it can only be used partially. **If you want to just iterate through
            all the results, aways use the `page_after_value` and `next_page_after_value`!**

            **NOTE #2**: In a request, you should never specify *both* `page` and
            `page_after_value` (at most one index can be provided).
        '''))

    @validator('page_size')
    def validate_page_size(cls, page_size):  # pylint: disable=no-self-argument
        assert page_size >= 0, 'page_size must be >= 0'
        return page_size

    @validator('order_by')
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        '''
        Override this in your Pagination class to ensure that a valid attribute is selected.
        This method has to be implemented!
        '''
        raise NotImplementedError('Validation of `order_by` not implemented!')

    @validator('page_after_value')
    def validate_page_after_value(cls, page_after_value, values):  # pylint: disable=no-self-argument
        '''
        Override this in your Pagination class to implement validation of the
        `page_after_value` value.
        This method has to be implemented!
        '''
        raise NotImplementedError('Validation of `page_after_value` not implemented!')

    @validator('page')
    def validate_page(cls, page, values):  # pylint: disable=no-self-argument
        if page is not None:
            assert page >= 1, 'page must be >= 1'
        return page

    @root_validator(skip_on_failure=True)
    def validate_values(cls, values):  # pylint: disable=no-self-argument
        # Because of a bug in pydantic (#2670), root validators can't be overridden, so
        # we invoke a class method, which *can* be overridden.
        return cls._root_validation(values)

    @classmethod
    def _root_validation(cls, values):
        page = values.get('page')
        page_after_value = values.get('page_after_value')
        page_size = values.get('page_size')
        assert page is None or page_after_value is None, 'Cannot specify both `page` and `page_after_value'
        if page_size == 0:
            assert page is None, 'Cannot specify `page` when `page_size` is set to 0'
            assert page_after_value is None, 'Cannot specify `page_after_value` when `page_size` is set to 0'
        return values

    def get_simple_index(self):
        '''
        If simple, index-based pagination is used, this method can be used to get the
        corresponding index (0-based). It will look on either `page` or `page_after_value`.
        If neither index is provided, we return 0 (i.e. the first index).
        '''
        if self.page is None and self.page_after_value is None:
            return 0
        if self.page is not None:
            return (self.page - 1) * self.page_size
        rv = int(self.page_after_value) + 1
        assert rv >= 0
        return rv


class PaginationResponse(Pagination):
    total: int = Field(
        ..., description=strip('''
        The total number of results that fit the given query. This is independent of
        any pagination and aggregations.
        '''))
    next_page_after_value: Optional[str] = Field(
        None, description=strip('''
        The *next* value to be used as `page_after_value` in a follow up requests, to get
        the next page of results. If no more results are available, `next_page_after_value`
        will not be set.
        '''))
    page_url: Optional[str] = Field(
        None, description=strip('''
        The url of the current page. Only applicable for GET requests.
        '''))
    next_page_url: Optional[str] = Field(
        None, description=strip('''
        The url to get the next page. Only applicable for GET requests.
        '''))
    prev_page_url: Optional[str] = Field(
        None, description=strip('''
        The url to get the previous page. **NOTE:** Only applicable for some API methods,
        (namely, where indexing by `page` is possible), and only for GET requests.
        '''))
    first_page_url: Optional[str] = Field(
        None, description=strip('''
        The url to get the first page. Only applicable for GET requests.
        '''))

    @validator('order_by')
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        # No validation - behaviour of this field depends on api method
        return order_by

    @validator('page_after_value')
    def validate_page_after_value(cls, page_after_value, values):  # pylint: disable=no-self-argument
        # No validation - behaviour of this field depends on api method
        return page_after_value

    @classmethod
    def _root_validation(cls, values):  # pylint: disable=no-self-argument
        # No validation
        return values

    def populate_urls(self, request: Request):
        '''
        Populates the urls (`page_url`, `next_page_url`, `first_page_url` from the
        request and `next_page_after_value`. Only applicable for GET requests.
        '''
        assert request.method.upper() == 'GET', 'Trying to populate urls, but method is not GET.'
        original_url = str(request.url)
        self.page_url = original_url
        if self.page_size:
            self.first_page_url = update_url_query_arguments(
                original_url, page=None, page_after_value=None)
        if self.next_page_after_value:
            self.next_page_url = update_url_query_arguments(
                original_url, page=None, page_after_value=self.next_page_after_value)
        if self.page and self.page > 1:
            self.prev_page_url = update_url_query_arguments(
                original_url, page=self.page - 1, page_after_value=None)

    def populate_simple_index_and_urls(self, request: Request):
        '''
        If simple, index-based pagination is used, this method can be used to populate
        the `page`, `page_after_value` and urls (if it is a GET request) automatically.
        Assumes that the field `total` is populated.
        '''
        if not self.page_size:
            self.page = 1
            self.page_after_value = None
            self.next_page_after_value = None
        else:
            ind = self.get_simple_index()
            self.page = ind // self.page_size + 1
            self.page_after_value = None if self.page == 1 else str(ind - 1)
            if self.page_size * self.page >= self.total:
                self.next_page_after_value = None
            else:
                self.next_page_after_value = str(ind + self.page_size - 1)

            if self.page < 1 or (
                    self.total == 0 and self.page != 1) or (
                    self.total > 0 and (self.page - 1) * self.page_size >= self.total):
                raise HTTPException(400, detail='Page out of range requested.')
        if request.method.upper() == 'GET':
            self.populate_urls(request)


class EntryBasedPagination(Pagination):
    order_by: Optional[str] = Field(
        calc_id,  # type: ignore
        description=strip('''
            The results are ordered by the values of this field. If omitted, default
            ordering is applied.
        '''))

    @validator('order_by')
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        if order_by is None:
            return order_by

        assert order_by in search_quantities, 'order_by must be a valid search quantity'
        quantity = search_quantities[order_by]
        assert quantity.definition.is_scalar, 'the order_by quantity must be a scalar'
        return order_by

    @validator('page_after_value')
    def validate_page_after_value(cls, page_after_value, values):  # pylint: disable=no-self-argument
        order_by = values.get('order_by', calc_id)
        if page_after_value is not None and order_by is not None and order_by != calc_id:
            if ':' not in page_after_value:
                page_after_value = '%s:' % page_after_value
        return page_after_value


class EntryPagination(EntryBasedPagination):
    page: Optional[int] = Field(
        None, description=strip('''
            For simple, index-based pagination, this should contain the number of the
            requested page (1-based). When provided in a request, this attribute can be
            used instead of `page_after_value` to jump to a particular results page.

            However, you can only retreive up to the 10.000th entry with a page number.
            Only one, `page_after_value` *or* `page` can be provided.
        '''))

    @validator('page')
    def validate_page(cls, page, values):  # pylint: disable=no-self-argument
        if page is not None:
            assert page > 0, 'Page has to be larger than 1.'
            assert page * values.get('page_size', 10) < 10000, 'Pagination by `page` is limited to 10.000 entries.'

        return page


entry_pagination_parameters = parameter_dependency_from_model(
    'entry_pagination_parameters', EntryPagination)


class AggregationPagination(EntryBasedPagination):
    order_by: Optional[str] = Field(
        None,  # type: ignore
        description=strip('''
            The results are ordered by the values of this field. If omitted, default
            ordering is applied.
        '''))

    @validator('page')
    def validate_page(cls, page, values):  # pylint: disable=no-self-argument
        assert page is None, 'Pagination by `page` is not possible for aggregations, use `page_after_value`'
        return page


class AggregatedEntities(BaseModel):
    size: Optional[pydantic.conint(gt=0)] = Field(  # type: ignore
        1, description=strip('''
        The maximum number of entries that should be returned for each value in the
        aggregation.
        '''))
    required: Optional[MetadataRequired] = Field(
        None, description=strip('''
        This allows to determined what fields should be returned for each entry.
        '''))


class Aggregation(BaseModel):
    quantity: AggregateableQuantity = Field(
        ..., description=strip('''
        The manatory name of the quantity for the aggregation. Aggregations
        can only be computed for those search metadata that have discrete values;
        an aggregation buckets entries that have the same value for this quantity.'''))
    pagination: Optional[AggregationPagination] = Field(
        AggregationPagination(), description=strip('''
        Only the data few values are returned for each API call. Pagination allows to
        get the next set of values based on the last value in the last call.
        '''))
    entries: Optional[AggregatedEntities] = Field(
        None, description=strip('''
        Optionally, a set of entries can be returned for each value.
        '''))


class StatisticsOrder(BaseModel):
    type_: Optional[AggregationOrderType] = Field(AggregationOrderType.entries, alias='type')
    direction: Optional[Direction] = Field(Direction.desc)


class Statistic(BaseModel):
    quantity: AggregateableQuantity = Field(
        ..., description=strip('''
        The manatory name of the quantity that the statistic is calculated for. Statistics
        can only be computed for those search metadata that have discrete values; a statistics
        aggregates a certain metric (e.g. the number of entries) over all entries were
        this quantity has the same value (bucket aggregation, think historgam here).

        There is one except and these are date/time values quantities (most notably `upload_time`).
        Here each statistic value represents an time interval. The interval can
        be determined via `datetime_interval`.'''))
    metrics: Optional[List[Metric]] = Field(
        [], description=strip('''
        By default the returned statistics will provide the number of entries for each
        value. You can add more metrics. For each metric an additional number will be
        provided for each value. Metrics are also based on search metadata. Depending on
        the metric the number will represent either a sum (`calculations` for the number
        of individual calculation in each code run) or an amount of different values
        (i.e. `materials` for the amount of different material hashes).'''))
    datetime_interval: Optional[pydantic.conint(gt=0)] = Field(  # type: ignore
        None, description=strip('''
        While statistics in general are only possible for quantities with discrete values,
        these is one exception. These are date/time values quantities (most notably `upload_time`).
        Here each statistic value represents an time interval.

        A date/time interval is a number of seconds greater than 0. This will only be used for
        date/time valued quantities (e.g. `upload_time`).
        '''))
    value_filter: Optional[pydantic.constr(regex=r'^[a-zA-Z0-9_\-\s]+$')] = Field(  # type: ignore
        None, description=strip('''
        An optional filter for values. Only values that contain the filter as substring
        will be part of the statistics.
        '''))
    size: Optional[pydantic.conint(gt=0)] = Field(  # type: ignore
        None, description=strip('''
        An optional maximum size of values in the statistics. The default depends on the
        quantity.
        '''))
    order: Optional[StatisticsOrder] = Field(
        StatisticsOrder(), description=strip('''
        The values in the statistics are either ordered by the entry count or by the
        natural ordering of the values.
        '''))

    @root_validator(skip_on_failure=True)
    def fill_default_size(cls, values):  # pylint: disable=no-self-argument
        if 'size' not in values or values['size'] is None:
            values['size'] = search_quantities[values['quantity'].value].statistic_size

        return values


class WithQueryAndPagination(WithQuery):
    pagination: Optional[EntryPagination] = Body(
        None,
        example={
            'page_size': 5,
            'order_by': 'upload_time'
        })


class EntriesMetadata(WithQueryAndPagination):
    required: Optional[MetadataRequired] = Body(
        None,
        example={
            'include': ['calc_id', 'mainfile', 'upload_id', 'authors', 'upload_time']
        })
    statistics: Optional[Dict[str, Statistic]] = Body(
        {},
        description=strip('''
            This allows to define additional statistics that should be returned.
            Statistics aggregate entries that show the same quantity values for a given quantity.
            A simple example is the number of entries for each `dft.code_name`. These statistics
            will be computed only over the query results. This allows to get an overview about
            query results.
        '''),
        example={
            'by_code_name': {
                'metrics': ['uploads', 'datasets'],
                'quantity': 'dft.code_name'
            }
        })
    aggregations: Optional[Dict[str, Aggregation]] = Body(
        {},
        example={
            'uploads': {
                'quantity': 'upload_id',
                'pagination': {
                    'page_size': 10,
                    'order_by': 'upload_time'
                },
                'entries': {
                    'size': 1,
                    'required': {
                        'include': ['mainfile']
                    }
                }
            }
        },
        description=strip('''
            Defines additional aggregations to return. An aggregation list entries
            for the values of a quantity, e.g. to get all uploads and their entries.
        '''))


class Files(BaseModel):
    ''' Configures the download of files. '''
    compress: Optional[bool] = Field(
        False, description=strip('''
        By default the returned zip file is not compressed. This allows to enable compression.
        Compression will reduce the rate at which data is provided, often below
        the rate of the compression. Therefore, compression is only sensible if the
        network connection is limited.'''))
    glob_pattern: Optional[str] = Field(
        None, description=strip('''
        An optional *glob* (or unix style path) pattern that is used to filter the
        returned files. Only files matching the pattern are returned. The pattern is only
        applied to the end of the full path. Internally
        [fnmatch](https://docs.python.org/3/library/fnmatch.html) is used.'''))
    re_pattern: Optional[str] = Field(
        None, description=strip('''
        An optional regexp that is used to filter the returned files. Only files matching
        the pattern are returned. The pattern is applied in search mode to the full
        path of the files. With `$` and `^` you can control if you want to match the
        whole path.

        A re pattern will replace a given glob pattern.'''))

    @validator('glob_pattern')
    def validate_glob_pattern(cls, glob_pattern):  # pylint: disable=no-self-argument
        # compile the glob pattern into re
        if glob_pattern is None:
            return None

        return re.compile(fnmatch.translate(glob_pattern) + r'$')

    @validator('re_pattern')
    def validate_re_pattern(cls, re_pattern):  # pylint: disable=no-self-argument
        # compile an re
        if re_pattern is None:
            return None
        try:
            return re.compile(re_pattern)
        except re.error as e:
            assert False, 'could not parse the re pattern: %s' % e

    @root_validator()
    def vaildate(cls, values):  # pylint: disable=no-self-argument
        # use the compiled glob pattern as re
        if values.get('re_pattern') is None:
            values['re_pattern'] = values.get('glob_pattern')
        return values


files_parameters = parameter_dependency_from_model(
    'files_parameters', Files)


ArchiveRequired = Union[str, Dict[str, Any]]


_archive_required_field = Body(
    '*',
    embed=True,
    description=strip('''
        The `required` part allows you to specify what parts of the requested archives
        should be returned. The NOMAD Archive is a hierarchical data format and
        you can *require* certain branches (i.e. *sections*) in the hierarchy.
        By specifing certain sections with specific contents or all contents (via
        the directive `"*"`), you can determine what sections and what quantities should
        be returned. The default is the whole archive, i.e., `"*"`.

        For example to specify that you are only interested in the `section_metadata`
        use:

        ```
        {
            "section_metadata": "*"
        }
        ```

        Or to only get the `energy_total` from each individual calculations, use:
        ```
        {
            "section_run": {
                "section_single_configuration_calculation": {
                    "energy_total": "*"
                }
            }
        }
        ```

        You can also request certain parts of a list, e.g. the last calculation:
        ```
        {
            "section_run": {
                "section_single_configuration_calculation[-1]": "*"
            }
        }
        ```

        These required specifications are also very useful to get workflow results.
        This works because we can use references (e.g. workflow to final result calculation)
        and the API will resolve these references and return the respective data.
        For example just the total energy value and reduced formula from the resulting
        calculation:
        ```
        {
            "section_workflow": {
                "calculation_result_ref": {
                    "energy_total": "*",
                    "single_configuration_calculation_to_system_ref": {
                        "chemical_composition_reduced": "*"
                    }
                }
            }
        }
        ```

        You can also resolve all references in a branch with the `include-resolved`
        directive. This will resolve all references in the branch, and also all references
        in referenced sections:
        ```
        {
            "section_workflow":
                "calculation_result_ref": "include-resolved"
            }
        }
        ```

        By default, the targets of "resolved" references are added to the archive at
        their original hierarchy positions.
        This means, all references are still references, but they are resolvable within
        the returned data, since they targets are now part of the data. Another option
        is to add
        `"resolve-inplace": true` to the root of required. Here, the reference targets will
        replace the references:
        ```
        {
            "resolve-inplace": true,
            "section_workflow":
                "calculation_result_ref": "include-resolved"
            }
        }
        ```
    '''),
    example={
        'section_run': {
            'section_single_configuration_calculation[-1]': {
                'energy_total': '*'
            },
            'section_system[-1]': '*'
        },
        'section_metadata': '*'
    })


class EntriesArchive(WithQueryAndPagination):
    required: Optional[ArchiveRequired] = _archive_required_field


class EntryArchiveRequest(BaseModel):
    required: Optional[ArchiveRequired] = _archive_required_field


class EntriesArchiveDownload(WithQuery):
    files: Optional[Files] = Body(None)


class EntriesRaw(WithQuery):
    pagination: Optional[EntryPagination] = Body(None)


class EntriesRawDownload(WithQuery):
    files: Optional[Files] = Body(
        None,
        example={
            'glob_pattern': 'vasp*.xml*'
        })


class StatisticResponse(Statistic):
    data: Dict[str, Dict[str, int]] = Field(
        None, description=strip('''
        The returned statistics data as dictionary. The key is a string representation of the values.
        The concrete type depends on the quantity that was used to create the statistics.
        Each dictionary value is a dictionary itself. The keys are the metric names the
        values the metric values. The key `entries` that gives the amount of entries with
        this value is always returned.'''))


class AggregationDataItem(BaseModel):
    data: Optional[List[Dict[str, Any]]] = Field(
        None, description=strip('''The entries that were requested for each value.'''))
    size: int = Field(
        None, description=strip('''The amount of entries with this value.'''))


class AggregationResponse(Aggregation):
    pagination: PaginationResponse  # type: ignore
    data: Dict[str, AggregationDataItem] = Field(
        None, description=strip('''
        The aggregation data as a dictionary. The key is a string representation of the values.
        The dictionary values contain the aggregated data depending if `entries` where
        requested.'''))


class CodeResponse(BaseModel):
    curl: str
    requests: str
    nomad_lab: Optional[str]


class EntriesMetadataResponse(EntriesMetadata):
    pagination: PaginationResponse  # type: ignore
    statistics: Optional[Dict[str, StatisticResponse]]  # type: ignore
    aggregations: Optional[Dict[str, AggregationResponse]]  # type: ignore
    data: List[Dict[str, Any]] = Field(
        None, description=strip('''
        The entries data as a list. Each item is a dictionary with the metadata for each
        entry.'''))
    code: Optional[CodeResponse]


class EntryRawFile(BaseModel):
    path: str = Field(None)
    size: int = Field(None)


class EntryRaw(BaseModel):
    calc_id: str = Field(None)
    upload_id: str = Field(None)
    mainfile: str = Field(None)
    files: List[EntryRawFile] = Field(None)


class EntriesRawResponse(EntriesRaw):
    pagination: PaginationResponse = Field(None)  # type: ignore
    data: List[EntryRaw] = Field(None)


class EntryMetadataResponse(BaseModel):
    entry_id: str = Field(None)
    required: MetadataRequired = Field(None)
    data: Dict[str, Any] = Field(
        None, description=strip('''A dictionary with the metadata of the requested entry.'''))


class EntryRawResponse(BaseModel):
    entry_id: str = Field(...)
    data: EntryRaw = Field(...)


class EntryArchive(BaseModel):
    calc_id: str = Field(None)
    upload_id: str = Field(None)
    parser_name: str = Field(None)
    archive: Dict[str, Any] = Field(None)


class EntriesArchiveResponse(EntriesArchive):
    pagination: PaginationResponse = Field(None)  # type: ignore
    data: List[EntryArchive] = Field(None)


class EntryArchiveResponse(EntryArchiveRequest):
    entry_id: str = Field(...)
    data: EntryArchive = Field(None)


class SearchResponse(EntriesMetadataResponse):
    es_query: Any = Field(
        None, description=strip('''The elasticsearch query that was used to retrieve the results.'''))
