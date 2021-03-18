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
This module provides an interface to elasticsearch. Other parts of NOMAD must not
interact with elasticsearch to maintain a clear coherent interface and allow for change.

Currently NOMAD uses two distinct pairs of elasticsearch indices. One based on the
old datamodel and metainfo layout based on metadata, domain, and encyclopedia sections (v0).
The other one based on the new section results (v1). Both provide an entry and derived materials
index. Write operations (index, publish, edit, lift embargo, delete) are common; defined
here in the module ``__init__.py``. There might be specific versions
in v0 and v1, but with a common interface here. Read operations are different and
should be used as per use-case directly from the ``v0`` and ``v1`` submodules.

Most common functions take an ``index`` keyword arg (either None, ``v0_index``, or
``v1_index``) and runs the function on either both, v0, or v1 entries index.

Most common functions also take an ``update_materials`` keyword arg with allows to
update the v1 materials index according to the performed changes.
'''

from typing import Union, List, Iterable
from enum import Enum
import json
from elasticsearch.exceptions import TransportError

from nomad import config, infrastructure, utils
from nomad.datamodel import EntryArchive, EntryMetadata

from . import v0
from . import v1
from .common import (
    _api_to_es_query, _owner_es_query,
    SearchError, AuthenticationRequiredError)


v0_index = 'v0'
v1_index = 'v1'


def run_on_both_indexes(func):
    '''
    A decorator that takes an ``index`` keyword arg (either None, ``v0_index``, or ``v1_index``)
    and runs the decorated function on either both, v0, or v1 entries index. The decorated
    function must take an elasticsearch index name as first argument.
    '''
    def wrapper(*args, index: str = None, **kwargs):
        if index is None:
            indices = [v0_index, v1_index]
        elif isinstance(index, List):
            indices = index
        else:
            indices = [index]

        for index in indices:
            if index == v0_index:
                index_name = config.elastic.index_name
            elif index == v1_index:
                index_name = config.elastic.entries_index
            else:
                index_name = index

            return func(index_name, *args, **kwargs)

    return wrapper


@run_on_both_indexes
def update_by_query(
        update_script: str,
        query: dict = None,
        owner: str = None,
        user_id: str = None,
        index: str = None,
        refresh: bool = False,
        **kwargs):
    '''
    Uses the given painless script to update the entries by given query.

    In most cases, the elasticsearch entry index should not be updated field by field;
    you should run `index_all` instead and fully replace documents from mongodb and
    archive files.

    This method provides a faster direct method to update individual fields, e.g. to quickly
    update fields for editing operations.
    '''
    if query is None:
        query = {}
    es_query = _api_to_es_query(query)
    if owner is not None:
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
            body=body, index=index)
    except TransportError as e:
        utils.get_logger(__name__).error(
            'es update_by_query script error', exc_info=e,
            es_info=json.dumps(e.info, indent=2))
        raise SearchError(e)

    if refresh:
        infrastructure.elastic_client.indices.refresh(index=index)

    return result


@run_on_both_indexes
def delete_by_query(
        query: dict,
        owner: str = None,
        user_id: str = None,
        index: str = None,
        refresh: bool = False):
    '''
    Deletes all entries that match the given query.
    '''
    if query is None:
        query = {}
    es_query = _api_to_es_query(query)
    es_query &= _owner_es_query(owner=owner, user_id=user_id)

    body = {
        'query': es_query.to_dict()
    }

    try:
        result = infrastructure.elastic_client.delete_by_query(
            body=body, index=index)
    except TransportError as e:
        utils.get_logger(__name__).error(
            'es delete_by_query error', exc_info=e,
            es_info=json.dumps(e.info, indent=2))
        raise SearchError(e)

    if refresh:
        infrastructure.elastic_client.indices.refresh(index=index)

    return result


@run_on_both_indexes
def refresh(index: str = None):
    '''
    Refreshes the specified indices.
    '''

    try:
        infrastructure.elastic_client.indices.refresh(index=index)
    except TransportError as e:
        utils.get_logger(__name__).error(
            'es delete_by_query error', exc_info=e,
            es_info=json.dumps(e.info, indent=2))
        raise SearchError(e)


def index(
        entries: Union[EntryArchive, List[EntryArchive]],
        index: str = None,
        update_materials: bool = False,
        refresh: bool = False):
    '''
    Index the given entries based on their archive. Either creates or updates the underlying
    elasticsearch documents.
    '''
    if index is None:
        indices = [v0_index, v1_index]
    else:
        indices = [index]

    if not isinstance(entries, list):
        entries = [entries]

    if v0_index in indices:
        assert not update_materials, 'update materials not supported for v0 index'
        v0.index_all([entry.section_metadata for entry in entries], do_refresh=refresh)

    if v1_index in indices:
        v1.index(entries=entries, update_materials=update_materials, refresh=refresh)


def publish(entries: Iterable[EntryMetadata], index: str = None):
    '''
    Publishes the given entries based on their entry metadata. Sets publishes to true,
    and updates most user provided metadata with a partial update.
    '''
    if index == v0_index or index is None:
        v0.publish(entries)
    if index == v1_index or index is None:
        v1.publish(entries)


def lift_embargo(query: dict, **kwargs):
    '''
    Removes the embargo flag based on a query.
    '''
    update_by_query(
        index=None,
        update_script='ctx._source.with_embargo = false;',
        query=query,
        **kwargs)


def delete_upload(upload_id: str, **kwargs):
    '''
    Deletes the given upload.
    '''
    delete_by_query(index=None, query=dict(upload_id=upload_id), **kwargs)


def delete_entry(entry_id: str, index: str = None, **kwargs):
    '''
    Deletes the given entry.
    '''
    if index == v0_index or index is None:
        delete_by_query(index=v0_index, query=dict(calc_id=entry_id), **kwargs)
    if index == v1_index or index is None:
        delete_by_query(index=v1_index, query=dict(entry_id=entry_id), **kwargs)
