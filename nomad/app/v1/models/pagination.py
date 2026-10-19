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

from fastapi import HTTPException, Request
from pydantic import Field, field_validator, model_validator
from pydantic_core import PydanticCustomError

from nomad.app.v1.utils import update_url_query_arguments
from nomad.config.models.pagination import Direction, PaginationBaseModel  # noqa: F401
from nomad.utils import strip


class Pagination(PaginationBaseModel):
    """Defines the order, size, and page of results."""

    page_after_value: str | None = Field(
        None,
        description=strip("""
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
        """),
    )
    page: int | None = Field(
        None,
        description=strip("""
            The number of the page (1-based). When provided in a request, this attribute
            can be used instead of `page_after_value` to jump to a particular results page.

            **NOTE #1**: the option to request pages by submitting the `page` number is
            limited. There are api calls where this attribute cannot be used for indexing,
            or where it can only be used partially. **If you want to just iterate through
            all the results, always use the `page_after_value` and `next_page_after_value`!**

            **NOTE #2**: Only one, `page`, `page_offset` or `page_after_value`, can be used.
        """),
    )
    page_offset: int | None = Field(
        None,
        description=strip("""
            The number of skipped entries. When provided in a request, this attribute
            can be used instead of `page_after_value` to jump to a particular results page.

            **NOTE #1**: the option to request pages by submitting the `page_offset` number is
            limited. There are api calls where this attribute cannot be used for indexing,
            or where it can only be used partially. **If you want to just iterate through
            all the results, always use the `page_after_value` and `next_page_after_value`!**

            **NOTE #2**: Only one, `page`, `page_offset` or `page_after_value`, can be used.
        """),
    )

    @field_validator('order_by')
    @classmethod
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        """
        Override this in your Pagination class to ensure that a valid attribute is selected.
        This method has to be implemented!
        """
        raise NotImplementedError('Validation of `order_by` not implemented!')

    @field_validator('page')
    @classmethod
    def validate_page(cls, page, values):  # pylint: disable=no-self-argument
        if page is not None and page < 1:
            raise PydanticCustomError('invalid_page', 'page must be >= 1')
        return page

    @field_validator('page_offset')
    @classmethod
    def validate_page_offset(cls, page_offset, values):  # pylint: disable=no-self-argument
        if page_offset is not None and page_offset < 0:
            raise PydanticCustomError('invalid_page_offset', 'page_offset must be >= 0')
        return page_offset

    @model_validator(mode='after')
    @classmethod
    def validate_values(cls, values):  # pylint: disable=no-self-argument
        # Because of a bug in pydantic (#2670), root validators can't be overridden, so
        # we invoke a class method, which *can* be overridden.
        return cls._root_validation(values)

    @classmethod
    def _root_validation(cls, values):
        page_offset = values.page_offset
        page = values.page
        page_after_value = values.page_after_value
        page_size = values.page_size

        n_offset_criteria = (
            (1 if page_offset else 0)
            + (1 if page else 0)
            + (1 if page_after_value else 0)
        )
        if n_offset_criteria > 1:
            raise PydanticCustomError(
                'multiple_pagination_criteria',
                'Can only specify one of: page_offset, page, or page_after_value',
            )

        if page_size == 0:
            if page_offset is not None:
                raise PydanticCustomError(
                    'invalid_pagination',
                    'Cannot specify page_offset when page_size is set to 0',
                )
            if page is not None:
                raise PydanticCustomError(
                    'invalid_pagination',
                    'Cannot specify page when page_size is set to 0',
                )
            if page_after_value is not None:
                raise PydanticCustomError(
                    'invalid_pagination',
                    'Cannot specify page_after_value when page_size is set to 0',
                )

        return values

    def get_simple_index(self):
        """
        If simple, index-based pagination is used, this method can be used to get the
        corresponding index (0-based). It will look on either `page` or `page_after_value`.
        If neither index is provided, we return 0 (i.e. the first index).
        """
        if self.page_offset is not None:
            return self.page_offset
        if self.page is not None:
            return (self.page - 1) * self.page_size
        if self.page_after_value is not None:
            rv = int(self.page_after_value) + 1
            assert rv >= 0
            return rv
        return 0

    def order_result(self, result):
        """
        Override this in your Pagination class to implement ordering of the results.
        This method has to be implemented!
        """
        raise NotImplementedError('Ordering of results not implemented!')

    def paginate_result(self, result, pick_value):
        """
        Override this in your Pagination class to implement pagination of the results.
        This method has to be implemented!
        """
        if self.page is not None:
            start = (self.page - 1) * self.page_size
            end = start + self.page_size
        elif self.page_offset is not None:
            start = self.page_offset
            end = start + self.page_size
        elif self.page_after_value is not None:
            start = 0
            for index, item in enumerate(result):
                if pick_value(item) == self.page_after_value:
                    start = index + 1
                    break
            end = start + self.page_size
        else:
            start, end = 0, self.page_size

        total_size = result.count()
        first, last = min(start, total_size), min(end, total_size)

        return [] if first == last else result[first:last]


class PaginationResponse(Pagination):
    total: int = Field(
        ...,
        description=strip(
            """
        The total number of results that fit the given query. This is independent of
        any pagination and aggregations.
        """
        ),
    )
    next_page_after_value: str | None = Field(
        None,
        description=strip(
            """
        The *next* value to be used as `page_after_value` in a follow up requests, to get
        the next page of results. If no more results are available, `next_page_after_value`
        will not be set.
        """
        ),
    )
    page_url: str | None = Field(
        None,
        description=strip(
            """
        The url of the current page. Only applicable for GET requests.
        """
        ),
    )
    next_page_url: str | None = Field(
        None,
        description=strip(
            """
        The url to get the next page. Only applicable for GET requests.
        """
        ),
    )
    prev_page_url: str | None = Field(
        None,
        description=strip(
            """
        The url to get the previous page. **NOTE:** Only applicable for some API methods,
        (namely, where indexing by `page` is possible), and only for GET requests.
        """
        ),
    )
    first_page_url: str | None = Field(
        None,
        description=strip(
            """
        The url to get the first page. Only applicable for GET requests.
        """
        ),
    )

    @field_validator('order_by')
    @classmethod
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        # No validation - behaviour of this field depends on api method
        return order_by

    @classmethod
    def _root_validation(cls, values):  # pylint: disable=no-self-argument
        # No validation
        return values

    def populate_urls(self, request: Request):
        """
        Populates the urls (`page_url`, `next_page_url`, `first_page_url` from the
        request and `next_page_after_value`. Only applicable for GET requests.
        """
        assert request.method.upper() == 'GET', (
            'Trying to populate urls, but method is not GET.'
        )
        original_url = str(request.url)
        self.page_url = original_url
        if self.page_size:
            self.first_page_url = update_url_query_arguments(
                original_url, page=None, page_after_value=None
            )
        if self.next_page_after_value:
            self.next_page_url = update_url_query_arguments(
                original_url, page=None, page_after_value=self.next_page_after_value
            )
        if self.page and self.page > 1:
            self.prev_page_url = update_url_query_arguments(
                original_url, page=self.page - 1, page_after_value=None
            )

    def populate_simple_index_and_urls(self, request: Request):
        """
        If simple, index-based pagination is used, this method can be used to populate
        the `page`, `page_after_value` and urls (if it is a GET request) automatically.
        Assumes that the field `total` is populated.
        """
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

            if (
                self.page < 1
                or (self.total == 0 and self.page != 1)
                or (0 < self.total <= (self.page - 1) * self.page_size)
            ):
                raise HTTPException(400, detail='Page out of range requested.')
        if request.method.upper() == 'GET':
            self.populate_urls(request)
