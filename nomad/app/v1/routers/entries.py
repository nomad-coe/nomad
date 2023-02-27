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
import math
from datetime import datetime

from typing import Optional, Set, Union, Dict, Iterator, Any, List
from fastapi import (
    APIRouter, Depends, Path, status, HTTPException, Request, Query as QueryParameter,
    Body)
from fastapi.responses import StreamingResponse
from fastapi.exceptions import RequestValidationError
from pydantic import BaseModel, Field, validator
import os.path
import io
import json
import orjson
from pydantic.main import create_model
from starlette.responses import Response
from joblib import Parallel, delayed, parallel_backend

from nomad import files, config, utils, metainfo, processing as proc
from nomad import datamodel
from nomad.datamodel import EditableUserMetadata
from nomad.files import StreamedFile, create_zipstream
from nomad.utils import strip
from nomad.archive import RequiredReader, RequiredValidationError, ArchiveQueryError
from nomad.search import AuthenticationRequiredError, SearchError, update_metadata as es_update_metadata
from nomad.search import search, QueryValidationError
from nomad.metainfo.elasticsearch_extension import entry_type

from .auth import create_user_dependency
from ..utils import (
    create_download_stream_zipped, create_download_stream_raw_file, browser_download_headers,
    DownloadItem, create_responses)
from ..models import (
    Aggregation, Pagination, PaginationResponse, MetadataPagination, TermsAggregation,
    WithQuery, WithQueryAndPagination, MetadataRequired, MetadataResponse, Metadata,
    MetadataEditRequest, Files, Query, User, Owner,
    QueryParameters, metadata_required_parameters, files_parameters, metadata_pagination_parameters,
    HTTPExceptionModel)


router = APIRouter()
default_tag = 'entries'
metadata_tag = 'entries/metadata'
raw_tag = 'entries/raw'
archive_tag = 'entries/archive'

logger = utils.get_logger(__name__)

query_parameters = QueryParameters(doc_type=entry_type)

archive_required_documentation = strip('''
The `required` part allows you to specify what parts of the requested archives
should be returned. The NOMAD Archive is a hierarchical data format and
you can *require* certain branches (i.e. *sections*) in the hierarchy.
By specifying certain sections with specific contents or all contents (via
the directive `"*"`), you can determine what sections and what quantities should
be returned. The default is the whole archive, i.e., `"*"`.

For example to specify that you are only interested in the `metadata`
use:

```json
{
    "metadata": "*"
}
```

Or to only get the `energy_total` from each individual entry, use:
```json
{
    "run": {
        "configuration": {
            "energy": "*"
        }
    }
}
```

You can also request certain parts of a list, e.g. the last calculation:
```json
{
    "run": {
        "calculation[-1]": "*"
    }
}
```

These required specifications are also very useful to get workflow results.
This works because we can use references (e.g. workflow to final result calculation)
and the API will resolve these references and return the respective data.
For example just the total energy value and reduced formula from the resulting
calculation:
```json
{
    "workflow": {
        "calculation_result_ref": {
            "energy": "*",
            "system_ref": {
                "value": {
                    "chemical_composition": "*"
                }
            }
        }
    }
}
```

You can also resolve all references in a branch with the `include-resolved`
directive. This will resolve all references in the branch, and also all references
in referenced sections:
```json
{
    "workflow":
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
```json
{
    "resolve-inplace": true,
    "workflow":
        "calculation_result_ref": "include-resolved"
    }
}
```
''')


ArchiveRequired = Union[str, Dict[str, Any]]

_archive_required_field = Body(
    '*',
    embed=True,
    description=archive_required_documentation,
    example={
        'run': {
            'calculation[-1]': {
                'energy': '*'
            },
            'system[-1]': '*'
        },
        'metadata': '*'
    })


class EntriesArchive(WithQueryAndPagination):
    required: Optional[ArchiveRequired] = _archive_required_field


class EntryArchiveRequest(BaseModel):
    required: Optional[ArchiveRequired] = _archive_required_field


class EntriesArchiveDownload(WithQuery):
    files: Optional[Files] = Body(None)


class EntriesRawDir(WithQuery):
    pagination: Optional[MetadataPagination] = Body(None)


class EntriesRaw(WithQuery):
    files: Optional[Files] = Body(
        None,
        example={
            'glob_pattern': 'vasp*.xml*'
        })


class EntryRawDirFile(BaseModel):
    path: str = Field(None)
    size: int = Field(None)


class EntryRawDir(BaseModel):
    entry_id: str = Field(None)
    upload_id: str = Field(None)
    mainfile: str = Field(None)
    mainfile_key: Optional[str] = Field(None)
    files: List[EntryRawDirFile] = Field(None)


class EntriesRawDirResponse(EntriesRawDir):
    pagination: PaginationResponse = Field(None)  # type: ignore
    data: List[EntryRawDir] = Field(None)


class EntryRawDirResponse(BaseModel):
    entry_id: str = Field(...)
    data: EntryRawDir = Field(...)


class EntryArchive(BaseModel):
    entry_id: str = Field(None)
    upload_id: str = Field(None)
    parser_name: str = Field(None)
    archive: Dict[str, Any] = Field(None)


class EntriesArchiveResponse(EntriesArchive):
    pagination: PaginationResponse = Field(None)  # type: ignore
    data: List[EntryArchive] = Field(None)


class EntryArchiveResponse(EntryArchiveRequest):
    entry_id: str = Field(...)
    data: EntryArchive = Field(None)


class EntryMetadataResponse(BaseModel):
    entry_id: str = Field(None)
    required: MetadataRequired = Field(None)
    data: Any = Field(
        None, description=strip('''The entry metadata as dictionary.'''))


class EntryMetadataEditActionField(BaseModel):
    value: str = Field(None, description='The value/values that is set as a string.')
    success: Optional[bool] = Field(None, description='If this can/could be done. Only in API response.')
    message: Optional[str] = Field(None, descriptin='A message that details the action result. Only in API response.')


EntryMetadataEditActions = create_model('EntryMetadataEditActions', **{  # type: ignore
    quantity.name: (
        Optional[EntryMetadataEditActionField]
        if quantity.is_scalar else Optional[List[EntryMetadataEditActionField]], None)
    for quantity in EditableUserMetadata.m_def.definitions
    if isinstance(quantity, metainfo.Quantity)
})


class EntryMetadataEdit(WithQuery):
    verify: Optional[bool] = Field(False, description='If true, no action is performed.')
    actions: EntryMetadataEditActions = Field(  # type: ignore
        None,
        description='Each action specifies a single value (even for multi valued quantities).')

    @validator('owner')
    def validate_query(cls, owner):  # pylint: disable=no-self-argument
        return Owner.user


class EntryMetadataEditResponse(EntryMetadataEdit):
    success: bool = Field(None, description='If the overall edit can/could be done. Only in API response.')
    message: str = Field(None, description='A message that details the overall edit result. Only in API response.')


_bad_owner_response = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Unauthorized. The given owner requires authorization,
        but no or bad authentication credentials are given.''')}

_bad_id_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Entry not found. The given id does not match any entry.''')}

_bad_path_response = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('File or directory not found.')}

_bad_edit_request = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('Edit request could not be executed.')}

_bad_edit_request_authorization = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('Not enough permissions to execute edit request.')}

_bad_edit_request_empty_query = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('No matching entries found.')}

_raw_response = 200, {
    'content': {'application/zip': {}},
    'description': strip('''
        A zip file with the requested raw files. The file is streamed.
        The content length is not known in advance.
    ''')}

_raw_file_response = 200, {
    'content': {'application/octet-stream': {}},
    'description': strip('''
        A byte stream with raw file contents. The content length is not known in advance.
        If the whole file is requested, the mime-type might be more specific, depending
        on the file contents.
    ''')}

_archives_download_response = 200, {
    'content': {'application/zip': {}},
    'description': strip('''
        A zip file with the requested archive files. The file is streamed.
        The content length is not known in advance.
    ''')}

_archive_download_response = 200, {
    'content': {'application/json': {}},
    'description': strip('''
        A json body with the requested archive.
    ''')}


_bad_archive_required_response = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The given required specification could not be understood.''')}


_bad_metadata_edit_response = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The given edit actions cannot be performed by you on the given query.''')}


def perform_search(*args, **kwargs):
    try:
        search_response = search(*args, **kwargs)
        search_response.es_query = None
        return search_response
    except QueryValidationError as e:
        raise RequestValidationError(errors=e.errors)
    except AuthenticationRequiredError as e:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=str(e))
    except SearchError as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Elasticsearch could not process your query: %s' % str(e))


@router.post(
    '/query', tags=[metadata_tag],
    summary='Search entries and retrieve their metadata',
    response_model=MetadataResponse,
    responses=create_responses(_bad_owner_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_entries_metadata_query(
        request: Request,
        data: Metadata,
        user: User = Depends(create_user_dependency())):

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
        aggregations=data.aggregations,
        user_id=user.user_id if user is not None else None)


@router.get(
    '', tags=[metadata_tag],
    summary='Search entries and retrieve their metadata',
    response_model=MetadataResponse,
    responses=create_responses(_bad_owner_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_entries_metadata(
        request: Request,
        with_query: WithQuery = Depends(query_parameters),
        pagination: MetadataPagination = Depends(metadata_pagination_parameters),
        required: MetadataRequired = Depends(metadata_required_parameters),
        user: User = Depends(create_user_dependency())):
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

    res = perform_search(
        owner=with_query.owner, query=with_query.query,
        pagination=pagination, required=required,
        user_id=user.user_id if user is not None else None)
    res.pagination.populate_urls(request)
    return res


def _do_exaustive_search(owner: Owner, query: Query, include: List[str], user: User) -> Iterator[Dict[str, Any]]:
    page_after_value = None
    while True:
        response = perform_search(
            owner=owner, query=query,
            pagination=MetadataPagination(page_size=100, page_after_value=page_after_value, order_by='upload_id'),
            required=MetadataRequired(include=include),
            user_id=user.user_id if user is not None else None)

        page_after_value = response.pagination.next_page_after_value

        for result in response.data:
            yield result

        if page_after_value is None or len(response.data) == 0:
            break


class _Uploads:
    '''
    A helper class that caches subsequent access to upload files the same upload.
    '''

    def __init__(self):
        self._upload_files = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def get_upload_files(self, upload_id: str) -> files.UploadFiles:
        if self._upload_files is not None and self._upload_files.upload_id != upload_id:
            self._upload_files.close()

        if self._upload_files is None or self._upload_files.upload_id != upload_id:
            self._upload_files = files.UploadFiles.get(upload_id)

        return self._upload_files

    def close(self):
        if self._upload_files is not None:
            self._upload_files.close()


def _create_entry_rawdir(entry_metadata: Dict[str, Any], uploads: _Uploads):
    entry_id = entry_metadata['entry_id']
    upload_id = entry_metadata['upload_id']
    mainfile = entry_metadata['mainfile']
    mainfile_key = entry_metadata.get('mainfile_key')

    upload_files = uploads.get_upload_files(upload_id)
    mainfile_dir = os.path.dirname(mainfile)

    files = []
    for path_info in upload_files.raw_directory_list(mainfile_dir, files_only=True):
        files.append(EntryRawDirFile(path=path_info.path, size=path_info.size))

    return EntryRawDir(
        entry_id=entry_id, upload_id=upload_id, mainfile=mainfile, mainfile_key=mainfile_key, files=files)


def _answer_entries_rawdir_request(
        owner: Owner, query: Query, pagination: MetadataPagination, user: User):

    if owner == Owner.all_:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            The owner=all is not allowed for this operation as it will search for entries
            that you might now be allowed to access.
            '''))

    search_response = perform_search(
        owner=owner, query=query,
        pagination=pagination,
        required=MetadataRequired(include=['entry_id', 'upload_id', 'mainfile']),
        user_id=user.user_id if user is not None else None)

    with _Uploads() as uploads:
        response_data = [_create_entry_rawdir(
            entry_metadata, uploads) for entry_metadata in search_response.data]

    return EntriesRawDirResponse(
        owner=search_response.owner,
        query=search_response.query,
        pagination=search_response.pagination,
        data=response_data)


def _answer_entries_raw_request(owner: Owner, query: Query, files: Files, user: User):
    if owner == Owner.all_:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            The owner=all is not allowed for this operation as it will search for entries
            that you might now be allowed to access.
            '''))

    response = perform_search(
        owner=owner, query=query,
        pagination=MetadataPagination(page_size=0),
        required=MetadataRequired(include=[]),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total > config.services.max_entry_download:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='The limit of maximum number of entries in a single download (%d) has been exeeded (%d).' % (
                config.services.max_entry_download, response.pagination.total))

    files_params = Files() if files is None else files
    search_includes = ['entry_id', 'upload_id', 'mainfile']

    try:
        # a generator of File objects to create the streamed zip from
        def download_items_generator():
            # go through all entries that match the query
            for entry_metadata in _do_exaustive_search(owner, query, include=search_includes, user=user):
                upload_id = entry_metadata['upload_id']
                mainfile = entry_metadata['mainfile']
                entry_metadata['mainfile'] = os.path.join(upload_id, mainfile)

                mainfile_dir = os.path.dirname(mainfile)
                yield DownloadItem(
                    upload_id=upload_id,
                    raw_path=mainfile_dir,
                    zip_path=os.path.join(upload_id, mainfile_dir),
                    entry_metadata=entry_metadata)

        # create the streaming response with zip file contents
        content = create_download_stream_zipped(
            download_items=download_items_generator(),
            re_pattern=files_params.re_pattern,
            recursive=False,
            create_manifest_file=True,
            compress=files_params.compress)
        return StreamingResponse(content, headers=browser_download_headers(
            filename='raw_files.zip',
            media_type='application/zip'))
    except Exception as e:
        logger.error('exception while streaming download', exc_info=e)
        raise


_entries_rawdir_query_docstring = strip('''
    Will perform a search and return a *page* of raw file metadata for entries fulfilling
    the query. This allows you to get a complete list of all rawfiles with their full
    path in their respective upload and their sizes. The first returned files for each
    entry, is their respective *mainfile*.

    Each entry on NOMAD has a set of raw files. These are the files in their original form,
    i.e. as provided by the uploader. More specifically, an entry has a *mainfile*, identified as
    parseable. For CMS entries, the mainfile is usually the main output file of the code. All other
    files in the same directory are considered the entries *auxiliary* no matter their role
    or if they were actually parsed by NOMAD.

    This operation supports the usual `owner`, `query`, and `pagination` parameters.
    ''')


@router.post(
    '/rawdir/query',
    tags=[raw_tag],
    summary='Search entries and get their raw files metadata',
    description=_entries_rawdir_query_docstring,
    response_model=EntriesRawDirResponse,
    responses=create_responses(_bad_owner_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_entries_rawdir_query(
        request: Request, data: EntriesRawDir, user: User = Depends(create_user_dependency())):

    return _answer_entries_rawdir_request(
        owner=data.owner, query=data.query, pagination=data.pagination, user=user)


@router.get(
    '/rawdir',
    tags=[raw_tag],
    summary='Search entries and get their raw files metadata',
    description=_entries_rawdir_query_docstring,
    response_model=EntriesRawDirResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_bad_owner_response))
async def get_entries_rawdir(
        request: Request,
        with_query: WithQuery = Depends(query_parameters),
        pagination: MetadataPagination = Depends(metadata_pagination_parameters),
        user: User = Depends(create_user_dependency())):

    res = _answer_entries_rawdir_request(
        owner=with_query.owner, query=with_query.query, pagination=pagination, user=user)
    res.pagination.populate_urls(request)
    return res


_entries_raw_query_docstring = strip('''
    This operation will perform a search and stream a .zip file with the raw files of the
    found entries.

    Each entry on NOMAD has a set of raw files. These are the files in their original form,
    i.e. as provided by the uploader. More specifically, an entry has a *mainfile*, identified as
    parseable. For CMS entries, the mainfile is usually the main output file of the code. All other
    files in the same directory are considered the entries *auxiliary* no matter their role
    or if they were actually parsed by NOMAD.

    After performing a search (that uses the same parameters as in all search operations),
    NOMAD will iterate through all results and create a .zip-file with all the entries'
    main and auxiliary files. The files will be organized in the same directory structure
    that they were uploaded in. The respective upload root directories are further prefixed
    with the `upload_id` of the respective uploads. The .zip-file will further contain
    a `manifest.json` with `upload_id`, `entry_id`, and `mainfile` of each entry.
    ''')


@router.post(
    '/raw/query',
    tags=[raw_tag],
    summary='Search entries and download their raw files',
    description=_entries_raw_query_docstring,
    response_class=StreamingResponse,
    responses=create_responses(_raw_response, _bad_owner_response))
async def post_entries_raw_query(
        data: EntriesRaw, user: User = Depends(create_user_dependency())):

    return _answer_entries_raw_request(
        owner=data.owner, query=data.query, files=data.files, user=user)


@router.get(
    '/raw',
    tags=[raw_tag],
    summary='Search entries and download their raw files',
    description=_entries_raw_query_docstring,
    response_class=StreamingResponse,
    responses=create_responses(_raw_response, _bad_owner_response))
async def get_entries_raw(
        with_query: WithQuery = Depends(query_parameters),
        files: Files = Depends(files_parameters),
        user: User = Depends(create_user_dependency(signature_token_auth_allowed=True))):

    return _answer_entries_raw_request(
        owner=with_query.owner, query=with_query.query, files=files, user=user)


def _read_archive(entry_metadata, uploads, required_reader: RequiredReader):
    entry_id = entry_metadata['entry_id']
    upload_id = entry_metadata['upload_id']
    upload_files = uploads.get_upload_files(upload_id)

    try:
        with upload_files.read_archive(entry_id) as archive:
            return {
                'entry_id': entry_id,
                'parser_name': entry_metadata['parser_name'],
                'archive': required_reader.read(archive, entry_id, upload_id)}
    except ArchiveQueryError as e:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, detail=str(e))


def _validate_required(required: ArchiveRequired, user) -> RequiredReader:
    try:
        return RequiredReader(required, user=user)
    except RequiredValidationError as e:
        raise HTTPException(
            status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=[dict(msg=e.msg, loc=['required'] + e.loc)])


def _read_entry_from_archive(entry: dict, uploads, required_reader: RequiredReader):
    entry_id, upload_id = entry['entry_id'], entry['upload_id']

    # all other exceptions are handled by the caller `_read_entries_from_archive`
    try:
        upload_files = uploads.get_upload_files(upload_id)

        with upload_files.read_archive(entry_id, True) as archive:
            entry['archive'] = required_reader.read(archive, entry_id, upload_id)
            return entry
    except ArchiveQueryError as e:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, detail=str(e))
    except KeyError as e:
        logger.error('missing archive', exc_info=e, entry_id=entry_id)

        return None


def _read_entries_from_archive(entries: Union[list, dict], required: ArchiveRequired, user):
    '''
    Takes pickleable arguments so that it can be offloaded to worker processes.

    It is important to ensure the return values are also pickleable.
    '''
    with _Uploads() as uploads:
        required_reader = _validate_required(required, user)

        if isinstance(entries, dict):
            return _read_entry_from_archive(entries, uploads, required_reader)

        return [_read_entry_from_archive(entry, uploads, required_reader) for entry in entries]


def _answer_entries_archive_request(
        owner: Owner, query: Query, pagination: MetadataPagination, required: ArchiveRequired,
        user: User):
    if owner == Owner.all_:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED, detail=strip(
                '''The owner=all is not allowed for this operation as it will search for entries
                that you might now be allowed to access.'''))

    if required is None:
        required = '*'

    search_response = perform_search(
        owner=owner, query=query,
        pagination=pagination,
        required=MetadataRequired(include=['entry_id', 'upload_id', 'parser_name']),
        user_id=user.user_id if user is not None else None)

    entries: list = [{
        'entry_id': entry['entry_id'], 'upload_id': entry['upload_id'],
        'parser_name': entry['parser_name']} for entry in search_response.data]

    # fewer than config.archive.min_entries_per_process entries per process is not useful
    # more than config.max_process_number processes is too much for the server
    number: int = min(
        int(math.ceil(len(entries) / config.archive.min_entries_per_process)),
        config.archive.max_process_number)

    if number <= 1:
        request_data: list = _read_entries_from_archive(entries, required, user)
    else:
        with parallel_backend('threading', n_jobs=number):
            request_data = Parallel()(delayed(
                _read_entries_from_archive)(i, required, user) for i in entries)

    return EntriesArchiveResponse(
        owner=search_response.owner,
        query=search_response.query,
        pagination=search_response.pagination,
        required=required,
        data=list(filter(None, request_data)))


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
        request: Request, data: EntriesArchive, user: User = Depends(create_user_dependency())):

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
        request: Request,
        with_query: WithQuery = Depends(query_parameters),
        pagination: MetadataPagination = Depends(metadata_pagination_parameters),
        user: User = Depends(create_user_dependency())):

    res = _answer_entries_archive_request(
        owner=with_query.owner, query=with_query.query, pagination=pagination,
        required=None, user=user)
    res.pagination.populate_urls(request)
    return res


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
        pagination=MetadataPagination(page_size=0),
        required=MetadataRequired(include=[]),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total > config.services.max_entry_download:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=(
                'The limit of maximum number of entries in a single download (%d) has been '
                'exeeded (%d).' % (config.services.max_entry_download, response.pagination.total)))

    manifest = []
    search_includes = ['entry_id', 'upload_id', 'parser_name']

    required_reader = RequiredReader('*', user=user)

    # a generator of StreamedFile objects to create the zipstream from
    def streamed_files():
        # go through all entries that match the query
        for entry_metadata in _do_exaustive_search(owner, query, include=search_includes, user=user):
            path = os.path.join(entry_metadata['upload_id'], '%s.json' % entry_metadata['entry_id'])
            try:
                archive_data = _read_archive(entry_metadata, uploads, required_reader)

                f = io.BytesIO(orjson.dumps(
                    archive_data, option=orjson.OPT_INDENT_2 | orjson.OPT_NON_STR_KEYS))

                yield StreamedFile(path=path, f=f, size=f.getbuffer().nbytes)
            except KeyError as e:
                logger.error('missing archive', entry_id=entry_metadata['entry_id'], exc_info=e)

            entry_metadata['path'] = path
            manifest.append(entry_metadata)

        # add the manifest at the end
        manifest_content = json.dumps(manifest, indent=2).encode()
        yield StreamedFile(path='manifest.json', f=io.BytesIO(manifest_content), size=len(manifest_content))

    with _Uploads() as uploads:
        # create the streaming response with zip file contents
        content = create_zipstream(streamed_files(), compress=files_params.compress)
        return StreamingResponse(content, headers=browser_download_headers(
            filename='archives.zip',
            media_type='application/zip'))


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
        _archives_download_response, _bad_owner_response, _bad_archive_required_response))
async def post_entries_archive_download_query(
        data: EntriesArchiveDownload, user: User = Depends(create_user_dependency())):

    return _answer_entries_archive_download_request(
        owner=data.owner, query=data.query, files=data.files, user=user)


@router.get(
    '/archive/download',
    tags=[archive_tag],
    summary='Search entries and download their archives',
    description=_entries_archive_download_docstring,
    response_class=StreamingResponse,
    responses=create_responses(
        _archives_download_response, _bad_owner_response, _bad_archive_required_response))
async def get_entries_archive_download(
        with_query: WithQuery = Depends(query_parameters),
        files: Files = Depends(files_parameters),
        user: User = Depends(create_user_dependency(signature_token_auth_allowed=True))):

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
        user: User = Depends(create_user_dependency())):
    '''
    Retrives the entry metadata for the given id.
    '''

    query = {'entry_id': entry_id}
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
    '/{entry_id}/rawdir',
    tags=[raw_tag],
    summary='Get the raw files metadata for an entry by its id',
    response_model=EntryRawDirResponse,
    responses=create_responses(_bad_id_response),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_entry_rawdir(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve raw data from.'),
        user: User = Depends(create_user_dependency())):
    '''
    Returns the file metadata for all input and output files (including auxiliary files)
    of the given `entry_id`. The first file will be the *mainfile*.
    '''
    query = dict(entry_id=entry_id)
    response = perform_search(
        owner=Owner.visible, query=query,
        required=MetadataRequired(include=['entry_id', 'upload_id', 'mainfile']),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total == 0:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The entry with the given id does not exist or is not visible to you.')

    with _Uploads() as uploads:
        return EntryRawDirResponse(
            entry_id=entry_id, data=_create_entry_rawdir(response.data[0], uploads))


@router.get(
    '/{entry_id}/raw',
    tags=[raw_tag],
    summary='Get the raw data of an entry by its id',
    response_class=StreamingResponse,
    responses=create_responses(_bad_id_response, _raw_response))
async def get_entry_raw(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve raw data from.'),
        files: Files = Depends(files_parameters),
        user: User = Depends(create_user_dependency(signature_token_auth_allowed=True))):
    '''
    Streams a .zip file with the raw files from the requested entry.
    '''
    query = dict(entry_id=entry_id)
    response = perform_search(
        owner=Owner.visible, query=query,
        required=MetadataRequired(include=['entry_id']),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total == 0:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The entry with the given id does not exist or is not visible to you.')

    return _answer_entries_raw_request(owner=Owner.visible, query=query, files=files, user=user)


@router.get(
    '/{entry_id}/raw/{path}',
    tags=[raw_tag],
    summary='Get the raw data of an entry by its id',
    response_class=StreamingResponse,
    responses=create_responses(_bad_id_response, _bad_path_response, _raw_file_response))
async def get_entry_raw_file(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve raw data from.'),
        path: str = Path(..., description='A relative path to a file based on the directory of the entry\'s mainfile.'),
        offset: Optional[int] = QueryParameter(
            0, ge=0, description=strip('''
                Integer offset that marks the start of the contents to retrieve. Default
                is the start of the file.''')),
        length: Optional[int] = QueryParameter(
            -1, ge=0, description=strip('''
                The amounts of contents in bytes to stream. By default, the remainder of
                the file is streamed.''')),
        decompress: Optional[bool] = QueryParameter(
            False, description=strip('''
                Attempt to decompress the contents, if the file is .gz or .xz.''')),
        user: User = Depends(create_user_dependency(signature_token_auth_allowed=True))):
    '''
    Streams the contents of an individual file from the requested entry.
    '''
    query = dict(entry_id=entry_id)
    response = perform_search(
        owner=Owner.visible, query=query,
        required=MetadataRequired(include=['entry_id', 'upload_id', 'mainfile']),
        user_id=user.user_id if user is not None else None)

    if response.pagination.total == 0:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The entry with the given id does not exist or is not visible to you.')

    entry_metadata = response.data[0]
    upload_id, mainfile = entry_metadata['upload_id'], entry_metadata['mainfile']
    # The user is allowed to access all files, because the entry is in the "visible" scope
    upload_files = files.UploadFiles.get(upload_id)

    entry_path = os.path.dirname(mainfile)
    path = os.path.join(entry_path, path)

    if not upload_files.raw_path_exists(path):
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='The requested file does not exist.')
    # We only provide a specific mime-type, if the whole file is requested. Otherwise,
    # it is unlikely that the provided contents will match the overall file mime-type.
    mime_type = 'application/octet-stream'
    if offset == 0 and length < 0:
        mime_type = upload_files.raw_file_mime_type(path)

    raw_file_content = create_download_stream_raw_file(
        upload_files, path, offset, length, decompress)
    return StreamingResponse(raw_file_content, media_type=mime_type)


def answer_entry_archive_request(
        query: Dict[str, Any], required: ArchiveRequired, user: User, entry_metadata=None):
    required_reader = _validate_required(required, user)

    if not entry_metadata:
        response = perform_search(
            owner=Owner.visible, query=query,
            required=MetadataRequired(include=['entry_id', 'upload_id', 'parser_name']),
            user_id=user.user_id if user is not None else None)

        if response.pagination.total == 0:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail='The entry does not exist or is not visible to you.')

        entry_metadata = response.data[0]

    entry_id = entry_metadata['entry_id']

    with _Uploads() as uploads:
        try:
            archive_data = _read_archive(entry_metadata, uploads, required_reader)['archive']
        except KeyError:
            raise HTTPException(
                status.HTTP_404_NOT_FOUND,
                detail='The entry does exist, but it has no archive.')

        return {
            'entry_id': entry_id,
            'required': required,
            'data': {
                'entry_id': entry_id,
                'upload_id': entry_metadata['upload_id'],
                'parser_name': entry_metadata['parser_name'],
                'archive': archive_data}}


@router.get(
    '/{entry_id}/archive',
    tags=[archive_tag],
    summary='Get the archive for an entry by its id',
    response_model=EntryArchiveResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_bad_id_response))
async def get_entry_archive(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve archive data from.'),
        user: User = Depends(create_user_dependency())):
    '''
    Returns the full archive for the given `entry_id`.
    '''
    return answer_entry_archive_request(dict(entry_id=entry_id), required='*', user=user)


@router.get(
    '/{entry_id}/archive/download',
    tags=[archive_tag],
    summary='Get the archive for an entry by its id as plain archive json',
    responses=create_responses(_bad_id_response, _archive_download_response))
async def get_entry_archive_download(
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve archive data from.'),
        ignore_mime_type: bool = QueryParameter(
            False, description=strip('''
                Sets the mime type specified in the response headers to `application/octet-stream`
                instead of the actual mime type (i.e. `application/json`).''')),
        user: User = Depends(create_user_dependency(signature_token_auth_allowed=True))):
    '''
    Returns the full archive for the given `entry_id`.
    '''
    response = answer_entry_archive_request(dict(entry_id=entry_id), required='*', user=user)
    archive = response['data']['archive']
    return StreamingResponse(
        io.BytesIO(json.dumps(archive, indent=2).encode()),
        headers=browser_download_headers(
            filename=f'{entry_id}.json',
            media_type='application/octet-stream' if ignore_mime_type else 'application/json'))


@router.post(
    '/{entry_id}/archive/query',
    tags=[archive_tag],
    summary='Get the archive for an entry by its id',
    response_model=EntryArchiveResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_bad_id_response, _bad_archive_required_response))
async def post_entry_archive_query(
        data: EntryArchiveRequest, user: User = Depends(create_user_dependency()),
        entry_id: str = Path(..., description='The unique entry id of the entry to retrieve archive data from.')):

    '''
    Returns a partial archive for the given `entry_id` based on the `required` specified
    in the body.
    '''
    return answer_entry_archive_request(dict(entry_id=entry_id), required=data.required, user=user)


def edit(query: Query, user: User, mongo_update: Dict[str, Any] = None, re_index=True) -> List[str]:
    # get all entries that have to change
    entry_ids: List[str] = []
    upload_ids: Set[str] = set()
    with utils.timer(logger, 'edit query executed'):
        all_entries = _do_exaustive_search(
            owner=Owner.user, query=query, include=['entry_id', 'upload_id'], user=user)

        for entry_dict in all_entries:
            entry_ids.append(entry_dict['entry_id'])
            upload_ids.add(entry_dict['upload_id'])

    # perform the update on the mongo db
    with utils.timer(logger, 'edit mongo update executed', size=len(entry_ids)):
        if mongo_update is not None:
            n_updated = proc.Entry.objects(entry_id__in=entry_ids).update(multi=True, **mongo_update)
            if n_updated != len(entry_ids):
                logger.error('edit repo did not update all entries', payload=mongo_update)

    # re-index the affected entries in elastic search
    with utils.timer(logger, 'edit elastic update executed', size=len(entry_ids)):
        if re_index:
            updated_metadata: List[datamodel.EntryMetadata] = []
            for entry in proc.Entry.objects(entry_id__in=entry_ids):
                entry_metadata = entry.mongo_metadata(entry.upload)
                # Ensure that updated fields are marked as "set", even if they are cleared
                entry_metadata.m_update_from_dict(mongo_update)
                # Add to list
                updated_metadata.append(entry_metadata)

            failed = es_update_metadata(updated_metadata, update_materials=True, refresh=True)

            if failed > 0:
                logger.error(
                    'edit repo with failed elastic updates',
                    payload=mongo_update, nfailed=failed)

    return list(upload_ids)


def get_quantity_values(quantity, **kwargs):
    '''
    Performs the search defined by `kwargs`, aggregated by quantity, and returns the encountered
    values of this quantity.
    '''
    response = perform_search(
        **kwargs,
        aggregations=dict(agg=Aggregation(terms=TermsAggregation(quantity=quantity))),
        pagination=Pagination(page_size=0))
    terms = response.aggregations['agg'].terms  # pylint: disable=no-member
    return [bucket.value for bucket in terms.data]


_editable_quantities = {
    quantity.name: quantity for quantity in EditableUserMetadata.m_def.definitions}


@router.post(
    '/edit_v0',
    tags=[metadata_tag],
    summary='Edit the user metadata of a set of entries',
    response_model=EntryMetadataEditResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_bad_metadata_edit_response))
async def post_entry_metadata_edit(
        response: Response,
        data: EntryMetadataEdit,
        user: User = Depends(create_user_dependency())):

    '''
    Performs or validates edit actions on a set of entries that match a given query.
    '''

    # checking the edit actions and preparing a mongo update on the fly
    query = data.query
    data = EntryMetadataEditResponse(**data.dict())
    data.query = query  # to dict from dict does not work with the op aliases in queries
    actions = data.actions
    verify = data.verify
    data.success = True
    mongo_update = {}
    main_author_ids = None
    has_error = False
    removed_datasets = None

    with utils.timer(logger, 'edit verified'):
        for action_quantity_name in actions.dict():  # type: ignore
            quantity_actions = getattr(actions, action_quantity_name, None)
            if quantity_actions is None:
                continue

            quantity = _editable_quantities.get(action_quantity_name)
            if quantity is None:
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail='Unknown quantity %s' % action_quantity_name)

            # TODO this does not work. Because the quantities are not in EditableUserMetadata
            # they are also not in the model and ignored by fastapi. This probably
            # also did not work in the old API.
            if action_quantity_name in ['main_author', 'upload_create_time']:
                if not user.is_admin():
                    raise HTTPException(
                        status_code=status.HTTP_400_BAD_REQUEST,
                        detail='Only the admin user can set %s' % quantity.name)

            if isinstance(quantity_actions, list) == quantity.is_scalar:
                raise HTTPException(
                    status_code=status.HTTP_400_BAD_REQUEST,
                    detail='Wrong shape for quantity %s' % action_quantity_name)

            if not isinstance(quantity_actions, list):
                quantity_actions = [quantity_actions]

            verify_reference = None
            if isinstance(quantity.type, metainfo.Reference):
                verify_reference = quantity.type.target_section_def.section_cls
            mongo_key = quantity.name
            has_error = False
            for action in quantity_actions:
                action.success = True
                action.message = None
                action_value = action.value
                action_value = action_value if action_value is None else action_value.strip()

                if action_quantity_name == 'with_embargo':
                    raise HTTPException(
                        status_code=status.HTTP_400_BAD_REQUEST,
                        detail='Updating the embargo flag on entry level is no longer allowed.')

                if action_value is None:
                    mongo_value = None

                elif action_value == '':
                    mongo_value = None

                elif verify_reference in [datamodel.User, datamodel.Author]:
                    try:
                        mongo_value = datamodel.User.get(user_id=action_value).user_id
                    except KeyError:
                        action.success = False
                        has_error = True
                        action.message = 'User does not exist'
                        continue

                    if main_author_ids is None:
                        main_author_ids = get_quantity_values(
                            quantity='main_author.user_id', owner=Owner.user, query=data.query, user_id=user.user_id)
                    if action_value in main_author_ids:
                        action.success = False
                        has_error = True
                        action.message = 'This user is already the main author of an entry in the query'
                        continue

                elif verify_reference == datamodel.Dataset:
                    try:
                        mongo_value = datamodel.Dataset.m_def.a_mongo.get(
                            user_id=user.user_id, dataset_name=action_value).dataset_id
                    except KeyError:
                        action.message = 'Dataset does not exist and will be created'
                        mongo_value = None
                        if not verify:
                            dataset = datamodel.Dataset(
                                dataset_id=utils.create_uuid(), user_id=user.user_id,
                                dataset_name=action_value, dataset_create_time=datetime.utcnow())
                            dataset.a_mongo.create()
                            mongo_value = dataset.dataset_id

                else:
                    mongo_value = action_value

                if len(quantity.shape) == 0:
                    mongo_update[mongo_key] = mongo_value
                else:
                    mongo_values = mongo_update.setdefault(mongo_key, [])
                    if mongo_value is not None:
                        if mongo_value in mongo_values:
                            action.success = False
                            has_error = True
                            action.message = 'Duplicate values are not allowed'
                            continue
                        mongo_values.append(mongo_value)

            if len(quantity_actions) == 0 and len(quantity.shape) > 0:
                mongo_update[mongo_key] = []

            if action_quantity_name == 'datasets':
                # check if datasets edit is allowed and if datasets have to be removed
                old_datasets = get_quantity_values(
                    quantity='datasets.dataset_id', owner=Owner.user, query=data.query, user_id=user.user_id)

                removed_datasets = []
                for dataset_id in old_datasets:
                    if dataset_id not in mongo_update.get(mongo_key, []):
                        removed_datasets.append(dataset_id)

                doi_ds = datamodel.Dataset.m_def.a_mongo.objects(
                    dataset_id__in=removed_datasets, doi__ne=None).first()
                if doi_ds is not None and not user.is_admin:
                    data.success = False
                    data.message = (data.message if data.message else '') + (
                        'Edit would remove entries from a dataset with DOI (%s) ' % doi_ds.dataset_name)
                    has_error = True

    # stop here, if client just wants to verify its actions
    if verify:
        return data

    # stop if the action were not ok
    if has_error:
        response.status_code = status.HTTP_400_BAD_REQUEST
        return data

    # perform the change
    mongo_update['last_edit_time'] = datetime.utcnow()
    edit(data.query, user, mongo_update, True)

    # remove potentially empty old datasets
    if removed_datasets is not None:
        for dataset in removed_datasets:
            if proc.Entry.objects(datasets=dataset).first() is None:
                datamodel.Dataset.m_def.a_mongo.objects(dataset_id=dataset).delete()

    return data


@router.post(
    '/edit',
    tags=[metadata_tag],
    summary='Edit the user metadata of a set of entries',
    response_model=MetadataEditRequest,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(
        _bad_edit_request, _bad_edit_request_authorization, _bad_edit_request_empty_query))
async def post_entries_edit(
        request: Request,
        data: MetadataEditRequest,
        user: User = Depends(create_user_dependency(required=True))):
    '''
    Updates the metadata of the specified entries.

    **Note:**
      - Only admins can edit some of the fields.
      - Only entry level attributes (like `comment`, `references` etc.) can be set using
        this endpoint; upload level attributes (like `upload_name`, `coauthors`, embargo
        settings, etc) need to be set through the endpoint **uploads/upload_id/edit**.
      - If the upload is published, the only operation permitted using this endpoint is to
        edit the entries in datasets that where created by the current user.
    '''
    edit_request_json = await request.json()
    try:
        verified_json = proc.MetadataEditRequestHandler.edit_metadata(edit_request_json, None, user)
        return verified_json
    except RequestValidationError:
        raise  # A problem which we have handled explicitly. Fastapi does json conversion.
    except Exception as e:
        # The upload is processing or some kind of unexpected error has occured
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
