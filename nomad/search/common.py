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
from elasticsearch_dsl import Q

from nomad import datamodel
from nomad.metainfo.elasticsearch_extension import entry_type, DocumentType

from nomad.app.v1 import models as api_models
from nomad.app.v1.models import MetadataRequired


class SearchError(Exception): pass


class AuthenticationRequiredError(Exception): pass


_entry_metadata_defaults = {
    quantity.name: quantity.default
    for quantity in datamodel.EntryMetadata.m_def.quantities  # pylint: disable=not-an-iterable
    if quantity.default not in [None, [], False, 0]
}


def _es_to_entry_dict(hit, required: MetadataRequired = None) -> Dict[str, Any]:
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

        return quantity_to_es(name, cast(api_models.Value, value))

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


def _owner_es_query(owner: str, user_id: str = None, doc_type: DocumentType = entry_type):
    def term_query(**kwargs):
        prefix = '' if doc_type == entry_type else 'entries.'
        return Q('term', **{
            (prefix + field): value for field, value in kwargs.items()})

    if owner == 'all':
        q = term_query(published=True)
        if user_id is not None:
            q = q | term_query(owners__user_id=user_id)
    elif owner == 'public':
        q = term_query(published=True) & term_query(with_embargo=False)
    elif owner == 'visible':
        q = term_query(published=True) & term_query(with_embargo=False)
        if user_id is not None:
            q = q | term_query(owners__user_id=user_id)
    elif owner == 'shared':
        if user_id is None:
            raise AuthenticationRequiredError('Authentication required for owner value shared.')

        q = term_query(owners__user_id=user_id)
    elif owner == 'user':
        if user_id is None:
            raise AuthenticationRequiredError('Authentication required for owner value user.')

        q = term_query(uploader__user_id=user_id)
    elif owner == 'staging':
        if user_id is None:
            raise AuthenticationRequiredError('Authentication required for owner value user')
        q = term_query(published=False) & term_query(owners__user_id=user_id)
    elif owner == 'admin':
        if user_id is None or not datamodel.User.get(user_id=user_id).is_admin:
            raise AuthenticationRequiredError('This can only be used by the admin user.')
        q = None
    elif owner is None:
        q = None
    else:
        raise KeyError('Unsupported owner value')

    if q is not None:
        return q
    return Q()
