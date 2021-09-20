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

Currently NOMAD uses one entry index and two distinct materials indices. The entries
index is based on two different mappings, once used by the old flask api (v0) and one
used by the new fastapi api (v1). The mappings are used at the same time and the documents
are merged. Write operations (index, publish, edit, lift embargo, delete) are common; defined
here in the module ``__init__.py``. Read operations are different and
should be used as per use-case directly from the ``v0`` and ``v1`` submodules.

Most common functions also take an ``update_materials`` keyword arg with allows to
update the v1 materials index according to the performed changes. TODO this is only
partially implemented.
'''

from typing import Union, List, Iterable, Any
from enum import Enum
import json
import elasticsearch
from elasticsearch.exceptions import TransportError

from nomad import config, infrastructure, utils
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.metainfo.elasticsearch_extension import index_entries, entry_type, entry_index

from .common import (
    _api_to_es_query, _owner_es_query,
    SearchError, AuthenticationRequiredError)


def update_by_query(
        update_script: str,
        query: Any = None,
        owner: str = None,
        user_id: str = None,
        index: str = None,
        refresh: bool = False,
        **kwargs):
    '''
    Uses the given painless script to update the entries by given query.

    In most cases, the elasticsearch entry index should not be updated field by field;
    you should run `index` instead and fully replace documents from mongodb and
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
            body=body, index=config.elastic.entries_index)
    except TransportError as e:
        utils.get_logger(__name__).error(
            'es update_by_query script error', exc_info=e,
            es_info=json.dumps(e.info, indent=2))
        raise SearchError(e)

    if refresh:
        _refresh()

    return result


def delete_by_query(
        query: dict,
        owner: str = None,
        user_id: str = None,
        update_materials: bool = False,
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
            body=body, index=config.elastic.entries_index)
    except TransportError as e:
        utils.get_logger(__name__).error(
            'es delete_by_query error', exc_info=e,
            es_info=json.dumps(e.info, indent=2))
        raise SearchError(e)

    if refresh:
        _refresh()

    if update_materials:
        # TODO update the matrials index at least for v1
        pass

    return result


def refresh():
    '''
    Refreshes the specified indices.
    '''

    try:
        infrastructure.elastic_client.indices.refresh(index=config.elastic.entries_index)
    except TransportError as e:
        utils.get_logger(__name__).error(
            'es delete_by_query error', exc_info=e,
            es_info=json.dumps(e.info, indent=2))
        raise SearchError(e)


_refresh = refresh


def index(
        entries: Union[EntryArchive, List[EntryArchive]],
        update_materials: bool = False,
        refresh: bool = True):
    '''
    Index the given entries based on their archive. Either creates or updates the underlying
    elasticsearch documents. If an underlying elasticsearch document already exists it
    will be fully replaced.
    '''
    if not isinstance(entries, list):
        entries = [entries]

    index_entries(entries=entries, update_materials=update_materials)

    if refresh:
        _refresh()


# TODO this depends on how we merge section metadata
def publish(entries: Iterable[EntryMetadata], index: str = None) -> int:
    '''
    Publishes the given entries based on their entry metadata. Sets publishes to true,
    and updates most user provided metadata with a partial update. Returns the number
    of failed updates.
    '''
    return update_metadata(
        entries, index=index, published=True, update_materials=True, refresh=True)


def update_metadata(
        entries: Iterable[EntryMetadata], index: str = None,
        update_materials: bool = False, refresh: bool = False,
        **kwargs) -> int:
    '''
    Update all given entries with their given metadata. Additionally apply kwargs.
    Returns the number of failed updates. This is doing a partial update on the underlying
    elasticsearch documents.
    '''

    def elastic_updates():
        for entry_metadata in entries:
            entry_archive = entry_metadata.m_parent
            if entry_archive is None:
                entry_archive = EntryArchive(metadata=entry_metadata)
            entry_doc = entry_type.create_index_doc(entry_archive)

            entry_doc.update(**kwargs)

            yield dict(
                doc=entry_doc,
                _id=entry_metadata.calc_id,
                _type=entry_index.doc_type.name,
                _index=entry_index.index_name,
                _op_type='update')

    updates = list(elastic_updates())
    _, failed = elasticsearch.helpers.bulk(
        infrastructure.elastic_client, updates, stats_only=True)

    if update_materials:
        # TODO update the matrials index at least for v1
        pass

    if refresh:
        _refresh()

    return failed


def delete_upload(upload_id: str, refresh: bool = False, **kwargs):
    '''
    Deletes the given upload.
    '''
    delete_by_query(query=dict(upload_id=upload_id), **kwargs)

    if refresh:
        _refresh()


def delete_entry(entry_id: str, index: str = None, refresh: bool = False, **kwargs):
    '''
    Deletes the given entry.
    '''
    delete_by_query(query=dict(calc_id=entry_id), **kwargs)

    if refresh:
        _refresh()
