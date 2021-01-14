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

from typing import Dict, Iterator, Any, List, Set, cast
from fastapi import APIRouter, Depends, Path, status, HTTPException
from fastapi.responses import StreamingResponse
import os.path
import io
import json
import orjson

from nomad import search, files, config, utils
from nomad.utils import strip
from nomad.archive import (
    query_archive, ArchiveQueryError, compute_required_with_referenced,
    read_partial_archives_from_mongo, filter_archive)

from .auth import get_optional_user
from ..utils import create_streamed_zipfile, File, create_responses
from ..models import (
    Pagination, WithQuery, MetadataRequired, EntriesMetadataResponse, EntriesMetadata,
    EntryMetadataResponse, query_parameters, metadata_required_parameters, Files, Query,
    pagination_parameters, files_parameters, User, Owner, HTTPExceptionModel, EntriesRaw,
    EntriesRawResponse, EntriesRawDownload, EntryRaw, EntryRawFile, EntryRawResponse,
    EntriesArchiveDownload, EntryArchiveResponse, EntriesArchive, EntriesArchiveResponse,
    ArchiveRequired)


router = APIRouter()
default_tag = 'entries'
metadata_tag = 'entries/metadata'
raw_tag = 'entries/raw'
archive_tag = 'entries/archive'

logger = utils.get_logger(__name__)


_bad_owner_response = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Unauthorized. The given owner requires authorization,
        but no or bad authentication credentials are given.''')}

_bad_id_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Entry not found. The given id does not match any entry.''')}

_raw_download_response = 200, {
    'content': {'application/zip': {}},
    'description': strip('''
        A zip file with the requested raw files. The file is streamed.
        The content length is not known in advance.
    ''')}

_archive_download_response = 200, {
    'content': {'application/zip': {}},
    'description': strip('''
        A zip file with the requested archive files. The file is streamed.
        The content length is not known in advance.
    ''')}


_bad_archive_required_response = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The given required specification could not be understood.''')}


def perform_search(*args, **kwargs):
    try:
        return search.search(*args, **kwargs)
    except search.AuthenticationRequiredError as e:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=str(e))
    except search.ElasticSearchError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Elasticsearch could not process your query: %s' % str(e))


@router.post(
    '/query', tags=['entries/metadata'],
    summary='Search entries and retrieve their metadata',
    response_model=EntriesMetadataResponse,
    responses=create_responses(_bad_owner_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_entries_metadata_query(
        data: EntriesMetadata,
        user: User = Depends(get_optional_user)):

    '''
    Executes a *query* and returns a *page* of the results with *required* result data
    as well as *statistics* and *aggregated* data.

    This is the basic search operation to retrieve metadata for entries that match
    certain search criteria (`query` and `owner`). All parameters (including `query`, `owner`)
    are optional. Look at the body schema or parameter documentation for more details.

    By default the *empty* search (that returns everything) is performed. Only a small
    page of the search results are returned at a time; use `pagination` in subsequent
    requests to retrive more data. Each entry has a lot of different *metadata*, use
    `required` to limit the data that is returned.

    The `statistics` and `aggregations` keys will further allow to return statistics
    and aggregated data over all search results.
    '''

    return perform_search(
        owner=data.owner,
        query=data.query,
        pagination=data.pagination,
        required=data.required,
        statistics=data.statistics,
        aggregations=data.aggregations,
        user_id=user.user_id if user is not None else None)


@router.get(
    '', tags=[metadata_tag],
    summary='Search entries and retrieve their metadata',
    response_model=EntriesMetadataResponse,
    responses=create_responses(_bad_owner_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_entries_metadata(
        with_query: WithQuery = Depends(query_parameters),
        pagination: Pagination = Depends(pagination_parameters),
        required: MetadataRequired = Depends(metadata_required_parameters),
        user: User = Depends(get_optional_user)):
    '''
    Executes a *query* and returns a *page* of the results with *required* result data.
    This is a version of `/entries/query`. Queries work a little different, because
    we cannot put complex queries into URL parameters.

    In addition to the `q` parameter (see parameter documentation for details), you can use all NOMAD
    search quantities as parameters, e.g. `?atoms=H&atoms=O`. Those quantities can be
    used with additional operators attached to their names, e.g. `?n_atoms__gte=3` for
    all entries with more than 3 atoms. Operators are `all`, `any`, `none`, `gte`,
    `gt`, `lt`, `lte`.
    '''

    return perform_search(
        owner=with_query.owner, query=with_query.query,
        pagination=pagination, required=required,
        user_id=user.user_id if user is not None else None)


def _do_exaustive_search(owner: Owner, query: Query, include: List[str], user: User) -> Iterator[Dict[str, Any]]:
    after = None
    while True:
        response = perform_search(
            owner=owner, query=query,
            pagination=Pagination(size=100, after=after, order_by='upload_id'),
            required=MetadataRequired(include=include),
            user_id=user.user_id if user is not None else None)

        after = response.pagination.next_after

        for result in response.data:
            yield result

        if after is None or len(response.data) == 0:
            break


class _Uploads():
    '''
    A helper class that caches subsequent access to upload files the same upload.
    '''
    def __init__(self):
        self._upload_files = None

    def get_upload_files(self, upload_id: str) -> files.UploadFiles:
        if self._upload_files is not None and self._upload_files.upload_id != upload_id:
            self._upload_files.close()

        if self._upload_files is None or self._upload_files.upload_id != upload_id:
            self._upload_files = files.UploadFiles.get(
                upload_id, is_authorized=lambda *args, **kwargs: True)

        return self._upload_files

    def close(self):
        if self._upload_files is not None:
            self._upload_files.close()


def _create_entry_raw(entry_metadata: Dict[str, Any], uploads: _Uploads):
    calc_id = entry_metadata['calc_id']
    upload_id = entry_metadata['upload_id']
    mainfile = entry_metadata['mainfile']

    upload_files = uploads.get_upload_files(upload_id)
    mainfile_dir = os.path.dirname(mainfile)

    files = []
    for file_name, file_size in upload_files.raw_file_list(directory=mainfile_dir):
        path = os.path.join(mainfile_dir, file_name)
        files.append(EntryRawFile(path=path, size=file_size))

    return EntryRaw(calc_id=calc_id, upload_id=upload_id, mainfile=mainfile, files=files)


def _answer_entries_raw_request(
        owner: Owner, query: Query, pagination: Pagination, user: User):

    if owner == Owner.all_:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            The owner=all is not allowed for this operation as it will search for entries
            that you might now be allowed to access.
            '''))

    search_response = perform_search(
        owner=owner, query=query,
        pagination=pagination,
        required=MetadataRequired(include=['calc_id', 'upload_id', 'mainfile']),
        user_id=user.user_id if user is not None else None)

    uploads = _Uploads()
    try:
        response_data = [
            _create_entry_raw(entry_metadata, uploads)
            for entry_metadata in search_response.data]
    finally:
        uploads.close()

    return EntriesRawResponse(
        owner=search_response.owner,
        query=search_response.query,
        pagination=search_response.pagination,
        data=response_data)


def _answer_entries_raw_download_request(owner: Owner, query: Query, files: Files, user: User):
    if owner == Owner.all_:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            The owner=all is not allowed for this operation as it will search for entries
            that you might now be allowed to access.
            '''))

    response = perform_search(
        owner=owner, query=query,
        pagination=Pagination(size=0),
        required=MetadataRequired(include=[]),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total > config.max_entry_download:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='The limit of maximum number of entries in a single download (%d) has been exeeded (%d).' % (
                config.max_entry_download, response.pagination.total))

    uploads = _Uploads()
    files_params = Files() if files is None else files
    manifest = []
    search_includes = ['calc_id', 'upload_id', 'mainfile']
    streamed_paths: Set[str] = set()

    try:
        # a generator of File objects to create the streamed zip from
        def raw_file_generator():
            # go through all entries that match the query
            for entry_metadata in _do_exaustive_search(owner, query, include=search_includes, user=user):
                upload_id = entry_metadata['upload_id']
                mainfile = entry_metadata['mainfile']

                upload_files = uploads.get_upload_files(upload_id)
                mainfile_dir = os.path.dirname(mainfile)

                # go through all files that belong to this entry
                all_filtered = True
                files = upload_files.raw_file_list(directory=mainfile_dir)
                for file_name, file_size in files:
                    path = os.path.join(mainfile_dir, file_name)

                    # apply the filter
                    if files_params.re_pattern is not None and not files_params.re_pattern.search(path):
                        continue
                    all_filtered = False

                    # add upload_id to path used in streamed zip
                    streamed_path = os.path.join(upload_id, path)

                    # check if already streamed
                    if streamed_path in streamed_paths:
                        continue
                    streamed_paths.add(streamed_path)

                    # yield the file
                    with upload_files.raw_file(path, 'rb') as f:
                        yield File(path=streamed_path, f=f, size=file_size)

                if not all_filtered or len(files) == 0:
                    entry_metadata['mainfile'] = os.path.join(upload_id, mainfile)
                    manifest.append(entry_metadata)

            # add the manifest at the end
            manifest_content = json.dumps(manifest).encode()
            yield File(path='manifest.json', f=io.BytesIO(manifest_content), size=len(manifest_content))

        # create the streaming response with zip file contents
        content = create_streamed_zipfile(raw_file_generator(), compress=files_params.compress)
        return StreamingResponse(content, media_type='application/zip')
    except Exception as e:
        logger.error('exception while streaming download', exc_info=e)
    finally:
        uploads.close()


_entries_raw_query_docstring = strip('''
    Will perform a search and return a *page* of raw file metadata for entries fulfilling
    the query. This allows you to get a complete list of all rawfiles with their full
    path in their respective upload and their sizes. The first returned files for each
    entry, is their respective *mainfile*.

    Each entry on NOMAD represents a set of raw files. These are the input and output
    files (as well as additional auxiliary files) in their original form, i.e. as
    provided by the uploader. More specifically, an entry represents a code-run identified
    by a certain *mainfile*. This is usually the main output file of the code. All other
    files in the same directory are considered the entries *auxiliary* no matter their role
    or if they were actually parsed by NOMAD.

    This operation supports the usual `owner`, `query`, and `pagination` parameters.
    ''')


@router.post(
    '/raw/query',
    tags=[raw_tag],
    summary='Search entries and get their raw files metadata',
    description=_entries_raw_query_docstring,
    response_model=EntriesRawResponse,
    responses=create_responses(_bad_owner_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_entries_raw_query(data: EntriesRaw, user: User = Depends(get_optional_user)):

    return _answer_entries_raw_request(
        owner=data.owner, query=data.query, pagination=data.pagination, user=user)


@router.get(
    '/raw',
    tags=[raw_tag],
    summary='Search entries and get raw their raw files metadata',
    description=_entries_raw_query_docstring,
    response_model=EntriesRawResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_bad_owner_response))
async def get_entries_raw(
        with_query: WithQuery = Depends(query_parameters),
        pagination: Pagination = Depends(pagination_parameters),
        user: User = Depends(get_optional_user)):

    return _answer_entries_raw_request(
        owner=with_query.owner, query=with_query.query, pagination=pagination, user=user)


_entries_raw_download_query_docstring = strip('''
    This operation will perform a search and stream a .zip file with raw input and output
    files of the found entries.

    Each entry on NOMAD represents a set of raw files. These are the input and output
    files (as well as additional auxiliary files) in their original form, i.e. as
    provided by the uploader. More specifically, an entry represents a code-run identified
    by a certain *mainfile*. This is usually the main output file of the code. All other
    files in the same directory are considered the entries *auxiliary* no matter their role
    or if they were actually parsed by NOMAD.

    After performing a search (that uses the same parameters as in all search operations),
    NOMAD will iterate through all results and create a .zip-file with all the entries'
    main and auxiliary files. The files will be organized in the same directory structure
    that they were uploaded in. The respective upload root directories are further prefixed
    with the `upload_id` of the respective uploads. The .zip-file will further contain
    a `manifest.json` with `upload_id`, `calc_id`, and `mainfile` of each entry.
    ''')


@router.post(
    '/raw/download/query',
    tags=[raw_tag],
    summary='Search entries and download their raw files',
    description=_entries_raw_download_query_docstring,
    response_class=StreamingResponse,
    responses=create_responses(_raw_download_response, _bad_owner_response))
async def post_entries_raw_download_query(
        data: EntriesRawDownload, user: User = Depends(get_optional_user)):

    return _answer_entries_raw_download_request(
        owner=data.owner, query=data.query, files=data.files, user=user)


@router.get(
    '/raw/download',
    tags=[raw_tag],
    summary='Search entries and download their raw files',
    description=_entries_raw_download_query_docstring,
    response_class=StreamingResponse,
    responses=create_responses(_raw_download_response, _bad_owner_response))
async def get_entries_raw_download(
        with_query: WithQuery = Depends(query_parameters),
        files: Files = Depends(files_parameters),
        user: User = Depends(get_optional_user)):

    return _answer_entries_raw_download_request(
        owner=with_query.owner, query=with_query.query, files=files, user=user)


def _read_archive(entry_metadata, uploads, required):
    calc_id = entry_metadata['calc_id']
    upload_id = entry_metadata['upload_id']
    upload_files = uploads.get_upload_files(upload_id)

    try:
        with upload_files.read_archive(calc_id) as archive:
            return {
                'calc_id': calc_id,
                'parser_name': entry_metadata['parser_name'],
                'archive': query_archive(archive, {calc_id: required})[calc_id]
            }
    except ArchiveQueryError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


def _answer_entries_archive_request(
        owner: Owner, query: Query, pagination: Pagination, required: ArchiveRequired,
        user: User):

    if owner == Owner.all_:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            The owner=all is not allowed for this operation as it will search for entries
            that you might now be allowed to access.
            '''))

    if required is None:
        required = '*'

    try:
        required_with_references = compute_required_with_referenced(required)
    except KeyError as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=(
            'The required specification contains an unknown quantity or section: %s' % str(e)))

    search_response = perform_search(
        owner=owner, query=query,
        pagination=pagination,
        required=MetadataRequired(include=['calc_id', 'upload_id', 'parser_name']),
        user_id=user.user_id if user is not None else None)

    if required_with_references is not None:
        # We can produce all the required archive data from the partial archives stored
        # in mongodb.
        entry_ids = [entry['calc_id'] for entry in search_response.data]
        partial_archives = cast(dict, read_partial_archives_from_mongo(entry_ids, as_dict=True))

    uploads = _Uploads()
    response_data = {}
    for entry_metadata in search_response.data:
        calc_id, upload_id = entry_metadata['calc_id'], entry_metadata['upload_id']

        archive_data = None
        if required_with_references is not None:
            try:
                partial_archive = partial_archives[calc_id]
                archive_data = filter_archive(required, partial_archive, transform=lambda e: e)
            except KeyError:
                # the partial archive might not exist, e.g. due to processing problems
                pass
            except ArchiveQueryError as e:
                detail = 'The required specification could not be understood: %s' % str(e)
                raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=detail)

        if archive_data is None:
            try:
                archive_data = _read_archive(entry_metadata, uploads, required)['archive']
            except KeyError as e:
                logger.error('missing archive', exc_info=e, calc_id=calc_id)
                continue

        response_data[calc_id] = {
            'calc_id': calc_id,
            'upload_id': upload_id,
            'parser_name': entry_metadata['parser_name'],
            'archive': archive_data}

    uploads.close()

    return EntriesArchiveResponse(
        owner=search_response.owner,
        query=search_response.query,
        pagination=search_response.pagination,
        data=list(response_data.values()))


_entries_archive_docstring = strip('''
    This operation will perform a search with the given `query` and `owner` and return
    the a *page* of `required` archive data. Look at the body schema or parameter documentation
    for more details. The **GET** version of this operation will only allow to provide
    the full archives.
    ''')


@router.post(
    '/archive/query',
    tags=[archive_tag],
    summary='Search entries and access their archives',
    description=_entries_archive_docstring,
    response_model=EntriesArchiveResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_bad_owner_response, _bad_archive_required_response))
async def post_entries_archive_query(
        data: EntriesArchive, user: User = Depends(get_optional_user)):

    return _answer_entries_archive_request(
        owner=data.owner, query=data.query, pagination=data.pagination,
        required=data.required, user=user)


@router.get(
    '/archive',
    tags=[archive_tag],
    summary='Search entries and access their archives',
    description=_entries_archive_docstring,
    response_model=EntriesArchiveResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_bad_owner_response, _bad_archive_required_response))
async def get_entries_archive_query(
        with_query: WithQuery = Depends(query_parameters),
        pagination: Pagination = Depends(pagination_parameters),
        user: User = Depends(get_optional_user)):

    return _answer_entries_archive_request(
        owner=with_query.owner, query=with_query.query, pagination=pagination,
        required=None, user=user)


def _answer_entries_archive_download_request(
        owner: Owner, query: Query, files: Files, user: User):

    if owner == Owner.all_:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            The owner=all is not allowed for this operation as it will search for entries
            that you might now be allowed to access.
            '''))

    files_params = Files() if files is None else files

    response = perform_search(
        owner=owner, query=query,
        pagination=Pagination(size=0),
        required=MetadataRequired(include=[]),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total > config.max_entry_download:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=(
                'The limit of maximum number of entries in a single download (%d) has been '
                'exeeded (%d).' % (config.max_entry_download, response.pagination.total)))

    uploads = _Uploads()
    manifest = []
    search_includes = ['calc_id', 'upload_id', 'parser_name']

    # a generator of File objects to create the streamed zip from
    def file_generator():
        # go through all entries that match the query
        for entry_metadata in _do_exaustive_search(owner, query, include=search_includes, user=user):
            path = os.path.join(entry_metadata['upload_id'], '%s.json' % entry_metadata['calc_id'])
            try:
                archive_data = _read_archive(entry_metadata, uploads, '*')

                f = io.BytesIO(orjson.dumps(
                    archive_data, option=orjson.OPT_INDENT_2 | orjson.OPT_NON_STR_KEYS))

                yield File(path=path, f=f, size=f.getbuffer().nbytes)
            except KeyError as e:
                logger.error('missing archive', calc_id=entry_metadata['calc_id'], exc_info=e)

            entry_metadata['path'] = path
            manifest.append(entry_metadata)

        # add the manifest at the end
        manifest_content = json.dumps(manifest).encode()
        yield File(path='manifest.json', f=io.BytesIO(manifest_content), size=len(manifest_content))

    try:
        # create the streaming response with zip file contents
        content = create_streamed_zipfile(file_generator(), compress=files_params.compress)
        return StreamingResponse(content, media_type='application/zip')
    finally:
        uploads.close()


_entries_archive_download_docstring = strip('''
    This operation will perform a search with the given `query` and `owner` and stream
    a .zip-file with the full archive contents for all matching entries. This is not
    paginated. Look at the body schema or parameter documentation for more details.
    ''')


@router.post(
    '/archive/download/query',
    tags=[archive_tag],
    summary='Search entries and download their archives',
    description=_entries_archive_download_docstring,
    response_class=StreamingResponse,
    responses=create_responses(
        _archive_download_response, _bad_owner_response, _bad_archive_required_response))
async def post_entries_archive_download_query(
        data: EntriesArchiveDownload, user: User = Depends(get_optional_user)):

    return _answer_entries_archive_download_request(
        owner=data.owner, query=data.query, files=data.files, user=user)


@router.get(
    '/archive/download',
    tags=[archive_tag],
    summary='Search entries and download their archives',
    description=_entries_archive_download_docstring,
    response_class=StreamingResponse,
    responses=create_responses(
        _archive_download_response, _bad_owner_response, _bad_archive_required_response))
async def get_entries_archive_download(
        with_query: WithQuery = Depends(query_parameters),
        files: Files = Depends(files_parameters),
        user: User = Depends(get_optional_user)):

    return _answer_entries_archive_download_request(
        owner=with_query.owner, query=with_query.query, files=files, user=user)


@router.get(
    '/{entry_id}', tags=[metadata_tag],
    summary='Get the metadata of an entry by its id',
    response_model=EntryMetadataResponse,
    responses=create_responses(_bad_id_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_entry_metadata(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve metadata from.'),
        required: MetadataRequired = Depends(metadata_required_parameters),
        user: User = Depends(get_optional_user)):
    '''
    Retrives the entry metadata for the given id.
    '''

    query = {'calc_id': entry_id}
    response = perform_search(owner=Owner.all_, query=query, required=required, user_id=user.user_id if user is not None else None)

    if response.pagination.total == 0:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The entry with the given id does not exist or is not visible to you.')

    return {
        'entry_id': entry_id,
        'required': required,
        'data': response.data[0]
    }


@router.get(
    '/{entry_id}/raw',
    tags=[raw_tag],
    summary='Get the raw files metadata for an entry by its id',
    response_model=EntryRawResponse,
    responses=create_responses(_bad_id_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_entry_raw(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve raw data from.'),
        files: Files = Depends(files_parameters),
        user: User = Depends(get_optional_user)):
    '''
    Returns the file metadata for all input and output files (including auxiliary files)
    of the given `entry_id`. The first file will be the *mainfile*.
    '''
    query = dict(calc_id=entry_id)
    response = perform_search(
        owner=Owner.visible, query=query,
        required=MetadataRequired(include=['calc_id', 'upload_id', 'mainfile']),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total == 0:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The entry with the given id does not exist or is not visible to you.')

    uploads = _Uploads()
    try:
        return EntryRawResponse(entry_id=entry_id, data=_create_entry_raw(response.data[0], uploads))
    finally:
        uploads.close()


@router.get(
    '/{entry_id}/raw/download',
    tags=[raw_tag],
    summary='Get the raw data of an entry by its id',
    response_class=StreamingResponse,
    responses=create_responses(_bad_id_response, _raw_download_response))
async def get_entry_raw_download(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve raw data from.'),
        files: Files = Depends(files_parameters),
        user: User = Depends(get_optional_user)):
    '''
    Streams a .zip file with the raw files from the requested entry.
    '''
    query = dict(calc_id=entry_id)
    response = perform_search(
        owner=Owner.visible, query=query,
        required=MetadataRequired(include=['calc_id']),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total == 0:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The entry with the given id does not exist or is not visible to you.')

    return _answer_entries_raw_download_request(owner=Owner.public, query=query, files=files, user=user)


@router.get(
    '/{entry_id}/archive',
    tags=[archive_tag],
    summary='Get the archive for an entry by its id',
    response_model=EntryArchiveResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_bad_id_response))
async def get_entry_archive(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve raw data from.'),
        user: User = Depends(get_optional_user)):
    '''
    Returns the full archive for the given `entry_id`.
    '''
    query = dict(calc_id=entry_id)
    response = perform_search(
        owner=Owner.visible, query=query,
        required=MetadataRequired(include=['calc_id', 'upload_id', 'parser_name']),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total == 0:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The entry with the given id does not exist or is not visible to you.')

    uploads = _Uploads()
    try:
        try:
            archive_data = _read_archive(response.data[0], uploads, required='*')
        except KeyError:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail='The entry with the given id does exist, but it has no archive.')

        return {
            'entry_id': entry_id,
            'data': archive_data['archive']}
    finally:
        uploads.close()
