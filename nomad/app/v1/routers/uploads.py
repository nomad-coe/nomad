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
import os
import io
import shutil
from datetime import datetime
from typing import Tuple, List, Set, Dict, Any, Optional, Union
from pydantic import BaseModel, Field, validator
from mongoengine.queryset.visitor import Q
from urllib.parse import unquote
from fastapi import (
    APIRouter, Request, File, UploadFile, status, Depends, Body, Path, Query as FastApiQuery,
    HTTPException)
from fastapi.responses import StreamingResponse
from fastapi.exceptions import RequestValidationError

from nomad import utils, config, files
from nomad.files import StagingUploadFiles, is_safe_relative_path, is_safe_basename
from nomad.bundles import BundleExporter, BundleImporter
from nomad.processing import Upload, Entry, ProcessAlreadyRunning, ProcessStatus, MetadataEditRequestHandler
from nomad.utils import strip
from nomad.search import search, search_iterator, refresh as search_refresh, QueryValidationError

from .auth import create_user_dependency, generate_upload_token
from ..models import (
    MetadataPagination, User, Direction, Pagination, PaginationResponse, HTTPExceptionModel,
    Files, files_parameters, Owner, WithQuery, MetadataRequired, MetadataEditRequest,
    restrict_query_to_upload)
from .entries import EntryArchiveResponse, answer_entry_archive_request
from ..utils import (
    parameter_dependency_from_model, create_responses, DownloadItem, browser_download_headers,
    create_download_stream_zipped, create_download_stream_raw_file, create_stream_from_string)

router = APIRouter()
default_tag = 'uploads'
metadata_tag = 'uploads/metadata'
raw_tag = 'uploads/raw'
archive_tag = 'uploads/archive'
action_tag = 'uploads/action'
bundle_tag = 'uploads/bundle'

logger = utils.get_logger(__name__)


class ProcData(BaseModel):
    process_running: bool = Field(
        description='If a process is running')
    current_process: Optional[str] = Field(
        description='Name of the current or last completed process')
    process_status: str = Field(
        description='The status of the current or last completed process')
    last_status_message: Optional[str] = Field(
        description='A short, human readable message from the current process, with '
                    'information about what the current process is doing, or information '
                    'about the completion (successful or not) of the last process, if no '
                    'process is currently running.')
    errors: List[str] = Field(
        descriptions='A list of error messages that occurred during the last processing')
    warnings: List[str] = Field(
        description='A list of warning messages that occurred during the last processing')
    complete_time: Optional[datetime] = Field(
        description='Date and time of the completion of the last process')

    class Config:
        orm_mode = True


class UploadProcData(ProcData):
    upload_id: str = Field(
        None,
        description='The unique id for the upload.')
    upload_name: Optional[str] = Field(
        description='The name of the upload. This can be provided during upload '
                    'using the `upload_name` query parameter.')
    upload_create_time: datetime = Field(
        None,
        description='Date and time of the creation of the upload.')
    main_author: str = Field(
        None,
        description=strip('The main author of the upload.'))
    coauthors: List[str] = Field(
        None,
        description=strip('A list of upload coauthors.'))
    reviewers: List[str] = Field(
        None,
        description=strip('A user provided list of reviewers.'))
    viewers: List[str] = Field(
        None,
        description=strip('All viewers (main author, upload coauthors, and reviewers)'))
    writers: List[str] = Field(
        None,
        description=strip('All writers (main author, upload coauthors)'))
    published: bool = Field(
        False,
        description='If this upload is already published.')
    published_to: List[str] = Field(
        None,
        description='A list of other NOMAD deployments that this upload was uploaded to already.')
    publish_time: Optional[datetime] = Field(
        'Date and time of publication, if the upload has been published.')
    with_embargo: bool = Field(
        description='If the upload has an embargo set (embargo_length not equal to zero).')
    embargo_length: int = Field(
        description='The length of the requested embargo, in months. 0 if no embargo is requested.')
    license: str = Field(
        description='The license under which this upload is distributed.')
    entries: int = Field(
        0,
        description='The number of identified entries in this upload.')
    upload_files_server_path: Optional[str] = Field(
        None, description='The path to the uploads files on the server.')


class EntryProcData(ProcData):
    entry_id: str = Field()
    entry_create_time: datetime = Field()
    mainfile: str = Field()
    mainfile_key: Optional[str] = Field()
    upload_id: str = Field()
    parser_name: str = Field()
    entry_metadata: Optional[dict] = Field()


class UploadProcDataPagination(Pagination):
    @validator('order_by')
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        if order_by is None:
            return 'upload_create_time'  # Default value
        assert order_by in ('upload_create_time', 'publish_time', 'upload_name', 'last_status_message'), 'order_by must be a valid attribute'
        return order_by

    @validator('page_after_value')
    def validate_page_after_value(cls, page_after_value, values):  # pylint: disable=no-self-argument
        # Validation handled elsewhere
        return page_after_value


upload_proc_data_pagination_parameters = parameter_dependency_from_model(
    'upload_proc_data_pagination_parameters', UploadProcDataPagination)


class EntryProcDataPagination(Pagination):
    @validator('order_by')
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        if order_by is None:
            return 'mainfile'  # Default value
        assert order_by in ('mainfile', 'parser_name', 'process_status', 'current_process'), 'order_by must be a valid attribute'
        return order_by

    @validator('page_after_value')
    def validate_page_after_value(cls, page_after_value, values):  # pylint: disable=no-self-argument
        # Validation handled elsewhere
        return page_after_value


entry_proc_data_pagination_parameters = parameter_dependency_from_model(
    'entry_proc_data_pagination_parameters', EntryProcDataPagination)


class UploadProcDataResponse(BaseModel):
    upload_id: str = Field(None, description=strip('''
        Unique id of the upload.'''))
    data: UploadProcData = Field(
        None, description=strip('''
        The upload data as a dictionary.'''))


class UploadProcDataQuery(BaseModel):
    upload_id: Optional[List[str]] = Field(
        description='Search for uploads matching the given id. Multiple values can be specified.')
    upload_name: Optional[List[str]] = Field(
        description='Search for uploads matching the given upload_name. Multiple values can be specified.')
    is_processing: Optional[bool] = Field(
        description=strip('''
            If True, only include currently processing uploads.
            If False, do not include currently processing uploads.
            If unset, include everything.'''))
    is_published: Optional[bool] = Field(
        description=strip('''
            If True: only include published uploads.
            If False: only include unpublished uploads.
            If unset: include everything.'''))


upload_proc_data_query_parameters = parameter_dependency_from_model(
    'upload_proc_data_query_parameters', UploadProcDataQuery)


class UploadProcDataQueryResponse(BaseModel):
    query: UploadProcDataQuery = Field()
    pagination: PaginationResponse = Field()
    data: List[UploadProcData] = Field(
        None, description=strip('''
        The upload data as a list. Each item is a dictionary with the data for each
        upload.'''))


class EntryProcDataResponse(BaseModel):
    entry_id: str = Field()
    data: EntryProcData = Field()


class EntryProcDataQueryResponse(BaseModel):
    pagination: PaginationResponse = Field()
    processing_successful: int = Field(
        None, description=strip('''
        Number of entries that has been processed successfully.
        '''))
    processing_failed: int = Field(
        None, description=strip('''
        Number of entries that failed to process.
        '''))
    upload: UploadProcData = Field(
        None, description=strip('''
        The upload processing data of the upload.
        '''))
    data: List[EntryProcData] = Field(
        None, description=strip('''
        The entries data as a list. Each item is a dictionary with the data for one entry.
        '''))


class RawDirPagination(Pagination):
    @validator('order_by')
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        assert not order_by, 'Cannot specify `order_by` for rawdir calls'
        return None

    @validator('page_after_value')
    def validate_page_after_value(cls, page_after_value, values):  # pylint: disable=no-self-argument
        # Validation handled elsewhere
        return page_after_value


rawdir_pagination_parameters = parameter_dependency_from_model(
    'rawdir_pagination_parameters', RawDirPagination, exclude=['order', 'order_by'])


class RawDirFileMetadata(BaseModel):
    ''' Metadata about a file '''
    name: str = Field()
    size: Optional[int] = Field()
    entry_id: Optional[str] = Field(description=strip('''
        If this is a mainfile: the ID of the corresponding entry.'''))
    parser_name: Optional[str] = Field(description=strip('''
        If this is a mainfile: the name of the matched parser.'''))


class RawDirElementMetadata(RawDirFileMetadata):
    ''' Metadata about an directory *element*, i.e. a file or a directory '''
    is_file: bool = Field()


class RawDirDirectoryMetadata(BaseModel):
    ''' Metadata about a directory '''
    name: str = Field()
    size: Optional[int] = Field()
    content: List[RawDirElementMetadata] = Field(
        example=[
            {'name': 'a_directory', 'is_file': False, 'size': 456},
            {'name': 'a_file.json', 'is_file': True, 'size': 123, 'entry_id': 'XYZ', 'parser_name': 'parsers/vasp'}])


class RawDirResponse(BaseModel):
    path: str = Field(example='The/requested/path')
    access: str = Field()
    file_metadata: Optional[RawDirFileMetadata] = Field()
    directory_metadata: Optional[RawDirDirectoryMetadata] = Field()
    pagination: Optional[PaginationResponse] = Field()


class ProcessingData(BaseModel):
    upload_id: str = Field()
    path: str = Field()
    entry_id: Optional[str] = Field()
    parser_name: Optional[str] = Field()
    entry: Optional[EntryProcData] = Field()
    archive: Optional[Dict[str, Any]] = Field()


class PutRawFileResponse(BaseModel):
    upload_id: str = Field(None, description=strip('''
        Unique id of the upload.'''))
    data: UploadProcData = Field(None, description=strip('''
        The upload data as a dictionary.'''))
    processing: Optional[ProcessingData] = Field(None, description=strip('''
        Information about the processing, including the entry (if one was generated) and
        [optionally] the archive data of this entry.'''))


class DeleteEntryFilesRequest(WithQuery):
    ''' Defines a request to delete entry files. '''
    owner: Optional[Owner] = Body('all')
    include_parent_folders: Optional[bool] = Field(
        False,
        description=strip('''
            If the delete operation should include not only the mainfiles of the selected entries,
            but also their folders.'''))


class UploadCommandExamplesResponse(BaseModel):
    upload_url: str = Field()
    upload_command: str = Field()
    upload_command_with_name: str = Field()
    upload_progress_command: str = Field()
    upload_command_form: str = Field()
    upload_tar_command: str = Field()


_not_authorized = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Unauthorized. Authorization is required, but no or bad authentication credentials provided.''')}

_not_authorized_to_upload = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Unauthorized. No credentials provided, or you do not have permissions to the
        specified upload.''')}

_not_authorized_to_entry = status.HTTP_401_UNAUTHORIZED, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Unauthorized. No credentials provided, or you do not have permissions to the
        specified upload or entry.''')}

_bad_request = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Bad request. The request could not be processed because of some error/invalid argument.''')}

_bad_pagination = status.HTTP_400_BAD_REQUEST, {
    'model': HTTPExceptionModel,
    'description': strip('''
        Bad request. Invalid pagination arguments supplied.''')}

_upload_not_found = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The specified upload could not be found.''')}

_entry_not_found = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The specified upload or entry could not be found.''')}

_upload_or_path_not_found = status.HTTP_404_NOT_FOUND, {
    'model': HTTPExceptionModel,
    'description': strip('''
        The specified upload, or a resource with the specified path within the upload,
        could not be found.''')}

_post_upload_response = 200, {
    'model': UploadProcDataResponse,
    'content': {
        'application/json': {},
        'text/plain': {'example': 'Thanks for uploading your data to nomad.'}
    },
    'description': strip('''
        A json structure with upload data, or a plain text information string.
        It will be a json structure if the request headers specifies `Accept = application/json`.''')}

_put_raw_file_response = 200, {
    'model': PutRawFileResponse,
    'content': {
        'application/json': {},
        'text/plain': {'example': 'Thanks for uploading your data to nomad.'}
    },
    'description': strip('''
        A json structure with upload data and possibly information from the processing,
        or a plain text information string.
        It will be a json structure if the request headers specifies `Accept = application/json`
        or if `wait_for_processing` is set.''')}


_raw_path_response = 200, {
    'content': {
        'application/octet-stream': {'example': 'file data'},
        'application/zip': {'example': '<zipped file or directory content>'}},
    'description': strip('''
        If `path` denotes a file: a stream with the file content, zipped if `compress = true`.
        If `path` denotes a directory, and `compress = true`, the directory content, zipped.''')}

_upload_bundle_response = 200, {
    'content': {
        'application/zip': {'example': '<zipped bundle data>'}}}


_thank_you_message = f'''
Thanks for uploading your data to nomad.
Go back to {config.gui_url()} and press
reload to see the progress on your upload
and publish your data.'''


@router.get(
    '/command-examples', tags=[default_tag],
    summary='Get example commands for shell based uploads.',
    response_model=UploadCommandExamplesResponse,
    responses=create_responses(_not_authorized),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_command_examples(user: User = Depends(create_user_dependency(required=True))):
    ''' Get url and example command for shell based uploads. '''
    token = generate_upload_token(user)
    api_url = config.api_url(ssl=config.services.https_upload, api='api/v1')
    upload_url = f'{api_url}/uploads?token={token}'
    upload_url_with_name = upload_url + '&upload_name=<name>'
    # Upload via streaming data tends to work much easier, e.g. no mime type issues, etc.
    # It is also easier for the user to unterstand IMHO.
    upload_command = f"curl -X POST '{upload_url}' -T <local_file>"
    rv = UploadCommandExamplesResponse(
        upload_url=upload_url,
        upload_command=upload_command,
        upload_command_form=f"curl -X POST '{upload_url}' -F file=@<local_file>",
        upload_command_with_name=f"curl -X POST '{upload_url_with_name}' -T <local_file>",
        upload_progress_command=upload_command + ' | xargs echo',
        upload_tar_command=f"tar -cf - <local_folder> | curl -# '{upload_url}' -X POST -T - | xargs echo"
    )
    return rv


@router.get(
    '', tags=[metadata_tag],
    summary='List uploads of authenticated user.',
    response_model=UploadProcDataQueryResponse,
    responses=create_responses(_not_authorized, _bad_pagination),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_uploads(
        request: Request,
        query: UploadProcDataQuery = Depends(upload_proc_data_query_parameters),
        pagination: UploadProcDataPagination = Depends(upload_proc_data_pagination_parameters),
        user: User = Depends(create_user_dependency(required=True))):
    '''
    Retrieves metadata about all uploads that match the given query criteria.
    '''
    # Build query
    mongo_query = Q()
    user_id = str(user.user_id)
    mongo_query &= Q(main_author=user_id) | Q(reviewers=user_id) | Q(coauthors=user_id)

    if query.upload_id:
        mongo_query &= Q(upload_id__in=query.upload_id)

    if query.upload_name:
        mongo_query &= Q(upload_name__in=query.upload_name)

    if query.is_processing is True:
        mongo_query &= Q(process_status__in=ProcessStatus.STATUSES_PROCESSING)
    elif query.is_processing is False:
        mongo_query &= Q(process_status__in=ProcessStatus.STATUSES_NOT_PROCESSING)

    if query.is_published is True:
        mongo_query &= Q(publish_time__ne=None)
    elif query.is_published is False:
        mongo_query &= Q(publish_time=None)

    # Fetch data from DB
    mongodb_query = Upload.objects.filter(mongo_query)
    # Create response
    start = pagination.get_simple_index()
    end = start + pagination.page_size

    order_by = pagination.order_by
    order_by_with_sign = order_by if pagination.order == Direction.asc else '-' + order_by
    if order_by == 'upload_create_time':
        order_by_args = [order_by_with_sign, 'upload_id']  # Use upload_id as tie breaker
    else:
        order_by_args = [order_by_with_sign, 'upload_create_time', 'upload_id']

    mongodb_query = mongodb_query.order_by(*order_by_args)

    data = [upload_to_pydantic(upload) for upload in mongodb_query[start:end]]

    pagination_response = PaginationResponse(total=mongodb_query.count(), **pagination.dict())
    pagination_response.populate_simple_index_and_urls(request)

    return UploadProcDataQueryResponse(
        query=query,
        pagination=pagination_response,
        data=data)


@router.get(
    '/{upload_id}', tags=[metadata_tag],
    summary='Get a specific upload',
    response_model=UploadProcDataResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_upload(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload to retrieve.'),
        user: User = Depends(create_user_dependency())):
    '''
    Fetches a specific upload by its upload_id.
    '''
    # Get upload (or throw exception if nonexistent/no access)
    upload = get_upload_with_read_access(upload_id, user, include_others=True)

    return UploadProcDataResponse(
        upload_id=upload_id,
        data=upload_to_pydantic(upload))


@router.get(
    '/{upload_id}/entries', tags=[metadata_tag],
    summary='Get the entries of the specific upload as a list',
    response_model=EntryProcDataQueryResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload, _bad_pagination),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_upload_entries(
        request: Request,
        upload_id: str = Path(
            ...,
            description='The unique id of the upload to retrieve entries for.'),
        pagination: EntryProcDataPagination = Depends(entry_proc_data_pagination_parameters),
        user: User = Depends(create_user_dependency())):
    '''
    Fetches the entries of a specific upload. Pagination is used to browse through the
    results.
    '''
    upload = get_upload_with_read_access(upload_id, user, include_others=True)

    order_by = pagination.order_by
    order_by_with_sign = order_by if pagination.order == Direction.asc else '-' + order_by

    start = pagination.get_simple_index()
    end = start + pagination.page_size

    # load upload's entries. Use entry_id as tie breaker for ordering.
    entries = list(upload.entries_sublist(start, end, order_by=(order_by_with_sign, 'entry_id')))
    failed_entries_count = upload.failed_entries_count

    # load entries's metadata from search
    metadata_entries_query = WithQuery(
        query={
            'entry_id:any': list(entry.entry_id for entry in entries)
        }).query
    metadata_entries = search(
        pagination=MetadataPagination(page_size=len(entries)),
        owner='admin' if user and user.is_admin else 'visible',
        user_id=user.user_id if user else None,
        query=metadata_entries_query)
    metadata_entries_map = {
        metadata_entry['entry_id']: metadata_entry
        for metadata_entry in metadata_entries.data}

    # convert data to pydantic
    data = []
    for entry in entries:
        pydantic_entry = entry_to_pydantic(entry)
        pydantic_entry.entry_metadata = metadata_entries_map.get(entry.entry_id)
        data.append(pydantic_entry)

    pagination_response = PaginationResponse(total=upload.total_entries_count, **pagination.dict())
    pagination_response.populate_simple_index_and_urls(request)

    return EntryProcDataQueryResponse(
        pagination=pagination_response,
        processing_successful=upload.processed_entries_count - failed_entries_count,
        processing_failed=failed_entries_count,
        upload=upload_to_pydantic(upload),
        data=data)


@router.get(
    '/{upload_id}/entries/{entry_id}', tags=[metadata_tag],
    summary='Get a specific entry for a specific upload',
    response_model=EntryProcDataResponse,
    responses=create_responses(_entry_not_found, _not_authorized_to_entry),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_upload_entry(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload.'),
        entry_id: str = Path(
            ...,
            description='The unique id of the entry, belonging to the specified upload.'),
        user: User = Depends(create_user_dependency(required=True))):
    '''
    Fetches a specific entry for a specific upload.
    '''
    upload = get_upload_with_read_access(upload_id, user)
    entry = upload.get_entry(entry_id)
    if not entry:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=strip('''
            An entry by that id could not be found in the specified upload.'''))

    data = entry_to_pydantic(entry, add_es_metadata=True, user=user)

    return EntryProcDataResponse(entry_id=entry_id, data=data)


@router.get(
    '/{upload_id}/rawdir/{path:path}', tags=[raw_tag],
    summary='Get the metadata for the raw file or folder located at the specified path in the specified upload.',
    response_model=RawDirResponse,
    responses=create_responses(_upload_or_path_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_upload_rawdir_path(
        request: Request,
        upload_id: str = Path(
            ...,
            description='The unique id of the upload.'),
        path: str = Path(
            ...,
            description='The path within the upload raw files.'),
        pagination: RawDirPagination = Depends(rawdir_pagination_parameters),
        include_entry_info: bool = FastApiQuery(
            False,
            description=strip('''
                If the fields `entry_id` and `parser_name` should be populated for all
                encountered mainfiles.''')),
        user: User = Depends(create_user_dependency(required=False))):
    '''
    For the upload specified by `upload_id`, gets the raw file or directory metadata
    located at the given `path`. The response will either contain a `file_metadata` or
    `directory_metadata` key. For files, basic data about the file is returned, such as its
    name and size. For directories, the response includes a list of elements
    (files and folders) in the directory. For directories, the result is paginated.
    '''
    # Get upload
    upload = get_upload_with_read_access(upload_id, user, include_others=True)
    upload_files = None
    try:
        # Get upload files
        upload_files = upload.upload_files
        if not upload_files.raw_path_exists(path):
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=strip('''
                Not found. Invalid path?'''))

        response = RawDirResponse(
            path=path.rstrip('/'),
            access='unpublished' if not upload.published else (
                'embargoed' if upload.embargo_length else 'public'))

        if upload_files.raw_path_is_file(path):
            # Path denotes a file
            response.file_metadata = RawDirFileMetadata(
                name=os.path.basename(path),
                size=upload_files.raw_file_size(path))
            if include_entry_info:
                entry: Entry = Entry.objects(upload_id=upload_id, mainfile=path, mainfile_key=None).first()
                if entry:
                    response.file_metadata.entry_id = entry.entry_id
                    response.file_metadata.parser_name = entry.parser_name
        else:
            # Path denotes a directory
            start = pagination.get_simple_index()
            end = start + pagination.page_size
            directory_list = upload_files.raw_directory_list(path)
            upload_files.close()
            content = []
            path_to_element: Dict[str, RawDirElementMetadata] = {}
            total = 0
            total_size = 0
            for i, path_info in enumerate(directory_list):
                total += 1
                total_size += path_info.size
                if start <= i < end:
                    element = RawDirElementMetadata(
                        name=os.path.basename(path_info.path),
                        is_file=path_info.is_file,
                        size=path_info.size)
                    content.append(element)
                    if include_entry_info:
                        path_to_element[path_info.path] = element

            if include_entry_info and content:
                for entry in Entry.objects(
                        upload_id=upload_id,
                        mainfile__in=path_to_element.keys(),
                        mainfile_key=None):
                    element = path_to_element[entry.mainfile]
                    element.entry_id = entry.entry_id
                    element.parser_name = entry.parser_name

            response.directory_metadata = RawDirDirectoryMetadata(
                name=os.path.basename(path),
                size=total_size,
                content=content)

            pagination_response = PaginationResponse(total=total, **pagination.dict())
            pagination_response.populate_simple_index_and_urls(request)
            response.pagination = pagination_response

        return response
    except Exception:
        if upload_files:
            upload_files.close()
        raise


@router.get(
    '/{upload_id}/raw/{path:path}', tags=[raw_tag],
    summary='Download the raw file or folder located at the specified path in the specified upload.',
    response_class=StreamingResponse,
    responses=create_responses(
        _raw_path_response, _upload_or_path_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_upload_raw_path(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload.'),
        path: str = Path(
            ...,
            description='The path within the upload raw files.'),
        files_params: Files = Depends(files_parameters),
        offset: Optional[int] = FastApiQuery(
            0,
            description=strip('''
                When dowloading individual files with `compress = false`, this can be
                used to seek to a specified position within the file in question. Default
                is 0, i.e. the start of the file.''')),
        length: Optional[int] = FastApiQuery(
            -1,
            description=strip('''
                When dowloading individual files with `compress = false`, this can be
                used to specify the number of bytes to read. By default, the value is -1,
                which means that the remainder of the file is streamed.''')),
        decompress: bool = FastApiQuery(
            False,
            description=strip('''
                Set if compressed files should be decompressed before streaming the
                content (that is: if there are compressed files *within* the raw files).
                Note, only some compression formats are supported.''')),
        ignore_mime_type: bool = FastApiQuery(
            False,
            description=strip('''
                Sets the mime type specified in the response headers to `application/octet-stream`
                instead of the actual mime type.''')),
        user: User = Depends(create_user_dependency(required=False, signature_token_auth_allowed=True))):
    '''
    For the upload specified by `upload_id`, gets the raw file or directory content located
    at the given `path`. The data is zipped if `compress = true`.

    It is possible to download both individual files and directories, but directories can
    only be downloaded if `compress = true`. When downloading a directory, it is also
    possible to specify `re_pattern`, `glob_pattern` or `include_files` to filter the files
    based on the file names.

    When downloading a file, you can specify `decompress` to attempt to decompress the data
    if the file is compressed before streaming it. You can also specify `offset` and `length`
    to download only a segment of the file (*Note:* `offset` and `length` does not work if
    `compress` is set to true).
    '''
    if files_params.compress and (offset != 0 or length != -1):
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
            Cannot specify `offset` or `length` when `compress` is true'''))
    # Get upload
    upload = get_upload_with_read_access(upload_id, user, include_others=True)
    # Get upload files
    upload_files = upload.upload_files
    try:
        if not upload_files.raw_path_exists(path):
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=strip('''
                Not found. Invalid path?'''))
        if upload_files.raw_path_is_file(path):
            # File
            if files_params.compress:
                media_type = 'application/zip'
                download_item = DownloadItem(
                    upload_id=upload_id,
                    raw_path=path,
                    zip_path=os.path.basename(path))
                content = create_download_stream_zipped(
                    download_item, upload_files, compress=True)
            else:
                if offset < 0:
                    raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                        Invalid offset provided.'''))
                if length <= 0 and length != -1:
                    raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                        Invalid length provided. Should be greater than 0, or -1 if the remainder
                        of the file should be read.'''))
                if ignore_mime_type or not (offset == 0 and length == -1):
                    media_type = 'application/octet-stream'
                else:
                    media_type = upload_files.raw_file_mime_type(path)
                content = create_download_stream_raw_file(
                    upload_files, path, offset, length, decompress)
            return StreamingResponse(content, headers=browser_download_headers(
                filename=os.path.basename(path) + ('.zip' if files_params.compress else ''),
                media_type=media_type))
        else:
            # Directory
            if not files_params.compress:
                raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                    Path is a directory, `compress` must be set to true'''))
            # Stream directory content, compressed.
            download_item = DownloadItem(
                upload_id=upload_id, raw_path=path, zip_path='')
            content = create_download_stream_zipped(
                download_item, upload_files,
                re_pattern=files_params.re_pattern, recursive=True,
                create_manifest_file=False, compress=True)
            return StreamingResponse(content, headers=browser_download_headers(
                (upload.upload_id if not path else os.path.basename(path.rstrip('/'))) + '.zip',
                media_type='application/zip'))
    except Exception as e:
        if not isinstance(e, HTTPException):
            logger.error('exception while streaming download', exc_info=e)
        upload_files.close()
        raise


@router.put(
    '/{upload_id}/raw/{path:path}', tags=[raw_tag],
    summary='Upload a raw file to the specified path (directory) in the specified upload.',
    response_class=StreamingResponse,
    responses=create_responses(
        _put_raw_file_response, _upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def put_upload_raw_path(
        request: Request,
        upload_id: str = Path(
            ...,
            description='The unique id of the upload.'),
        path: str = Path(
            ...,
            description='The path within the upload raw files.'),
        file: List[UploadFile] = File(None),
        local_path: str = FastApiQuery(
            None,
            description=strip('''
            Internal/Admin use only.''')),
        file_name: str = FastApiQuery(
            None,
            description=strip('''
            Specifies the name of the file, when using method 2.''')),
        overwrite_if_exists: bool = FastApiQuery(
            True,
            description=strip('''
            If set to True (default), overwrites the file if it already exists.''')),
        copy_or_move: str = FastApiQuery(
            None,
            description=strip('''
             If moving or copying a file within the same upload,
             specify which operation to do: move or copy''')),
        copy_or_move_source_path: str = FastApiQuery(
            None,
            description=strip('''
            If moving or copying a file within the same upload, specify the path to the
            source file.''')),
        wait_for_processing: bool = FastApiQuery(
            False,
            description=strip('''
            Waits for the processing to complete and return information about the outcome
            in the response (**USE WITH CARE**).''')),
        include_archive: bool = FastApiQuery(
            False,
            description=strip('''
            If the archive data should be included in the response when using
            `wait_for_processing` (**USE WITH CARE**).''')),
        entry_hash: str = FastApiQuery(
            None,
            description=strip('''
            The hash code of the not modified entry.''')),
        user: User = Depends(create_user_dependency(required=True, upload_token_auth_allowed=True))):
    '''
    Upload one or more files to the directory specified by `path` in the upload specified by `upload_id`.

    When uploading a zip or tar archive, it will first be extracted, and the content will be
    *merged* with the existing content, i.e. new files are added, and if there is a collision
    (an old file with the same path and name as one of the new files), the old file will
    be overwritten, but the rest of the old files will remain untouched. If the file is not
    a zip or tar archive, the file will just be uploaded as it is, overwriting the existing
    file if there is one.

    The `path` should denote a directory. The empty string gives the "root" directory.

    If a single file is uploaded (and it is not a zip or tar archive), it is possible to specify
    `wait_for_processing`. This means that the file (and only this file) will be matched and
    processed, and information about the outcome will be returned with the response. **NOTE**:
    this should be used with caution! When this option is set, the call will block until
    processing is complete, which may take some time. Also note, that just processing the
    new/modified file may not be enough in some cases (since adding/modifying a file somewhere
    in the directory structure may affect other entries). Also note that
    processing.entry.entry_metadata will not be populated in the response.

    There are two basic ways to upload files: in the multipart-formdata or streaming the
    file data in the http body. Both are supported. Note, however, that the second method
    only allows the upload of a single file, and that it does not transfer a filename. If a
    transfer is made using method 2, you can specify the query argument `file_name` to name it.
    This *needs* to be specified when using method 2, unless you are uploading a zip/tar file
    (for zip/tar files the names don't matter since they are extracted). See the POST `uploads`
    endpoint for examples of curl commands for uploading files.

    Also, this path can be used to copy/move a file from one directory to another. Three
    query parameters are required for a successful operation: 1) `copy_or_move` param to specify
    if the file needs to be moved (if set to move then the original file will be removed), 2)
    `file_name` param that contains the new name for the file moved/copied file and 3) `copy_or_move_source_path`
    param that contains the path of the original/existing local file to be copied or moved.
    '''
    if include_archive and not wait_for_processing:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='`include_archive` requires `wait_for_processing`.')

    upload = _get_upload_with_write_access(upload_id, user, include_published=False)

    if not is_safe_relative_path(path):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Bad path provided.')

    if copy_or_move is not None or copy_or_move_source_path is not None:
        if copy_or_move not in ['copy', 'move']:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='The copy_or_move query parameter should be one of \'copy\' or \'move\' options.')

        if copy_or_move is None or copy_or_move_source_path is None or file_name is None:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='''For a successful copy/move operation, all three query parameters:
                 file_name, copy_or_move and copy_or_move_source_path are required.''')

        if not is_safe_basename(file_name):
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='Bad file name provided')

        if not is_safe_relative_path(copy_or_move_source_path):
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='Bad source path provided.')

    upload_paths, method = await _get_files_if_provided(
        upload_id, request, file, local_path, file_name, user)

    if not upload_paths and not (copy_or_move and copy_or_move_source_path and file_name):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Either an upload file or the query parameters for moving/copying a file should be provided.')

    if entry_hash:
        upload_path = upload_paths[0]
        full_path = os.path.join(path, os.path.basename(upload_path))
        entry_id = utils.generate_entry_id(upload_id, full_path)
        entry = upload.get_entry(entry_id)
        if entry and entry_hash != entry.entry_hash or not entry:
            raise HTTPException(
                status_code=status.HTTP_409_CONFLICT,
                detail='The provided hash did not match the current file.')

    upload_files = StagingUploadFiles(upload_id)

    for upload_path in upload_paths:
        decompress = files.auto_decompress(upload_path)
        if decompress == 'error':
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='Cannot extract file. Bad file format or file extension?')
        if not decompress and not overwrite_if_exists:
            full_path = os.path.join(path, os.path.basename(upload_path))
            if upload_files.raw_path_exists(full_path):
                raise HTTPException(
                    status_code=status.HTTP_409_CONFLICT,
                    detail='The provided path already exists and overwrite_if_exists is set to False.')

    if not wait_for_processing:
        # Process on worker (normal case)
        if copy_or_move:  # the case for move/copy an existing file
            path_to_target_file = os.path.join(path, file_name)
            if upload_files.raw_path_exists(path_to_target_file):
                raise HTTPException(
                    status_code=status.HTTP_409_CONFLICT,
                    detail='The provided path already exists.')
            if not upload_files.raw_path_exists(path):
                raise HTTPException(
                    status_code=status.HTTP_404_NOT_FOUND,
                    detail='No file or folder with that path found.')
            file_operations: Any = [
                dict(
                    op=copy_or_move.upper(),
                    path_to_existing_file=copy_or_move_source_path,
                    path_to_target_file=path_to_target_file)]
        else:
            file_operations = [
                dict(op='ADD', path=upload_path, target_dir=path, temporary=(method != 0))
                for upload_path in upload_paths]

        # Initiate processing
        try:
            upload.process_upload(file_operations=file_operations, only_updated_files=True)
        except ProcessAlreadyRunning:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='The upload is currently blocked by another process.')
        # Create response
        if request.headers.get('Accept') == 'application/json':
            response = PutRawFileResponse(
                upload_id=upload_id,
                data=upload_to_pydantic(upload))
            response_text = response.json()
            media_type = 'application/json'
        else:
            response_text = _thank_you_message
            media_type = 'text/plain'
    else:
        # Process locally
        if copy_or_move:  # case for move/copy an existing file
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='Cannot move/copy the file with wait_for_processing set to true.')

        if len(upload_paths) != 1 or decompress:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='`wait_for_processing` can only be used with single files, and not with compressed files.')

        upload_path = upload_paths[0]
        full_path = os.path.join(path, os.path.basename(upload_path))
        try:
            reprocess_settings = config.Reprocess(
                index_individual_entries=True, reprocess_existing_entries=True)
            entry = upload.put_file_and_process_local(
                upload_path, path, reprocess_settings=reprocess_settings)

            search_refresh()

            archive = None
            if entry and entry.process_status == ProcessStatus.SUCCESS and include_archive:
                # NOTE: We can't rely on ES to get the metadata for the entry, since it may
                # not have hade enough time to update its index etc. For now, we will just
                # ignore this, as we do not need it.
                entry_metadata = dict(
                    upload_id=upload_id,
                    entry_id=entry.entry_id,
                    parser_name=entry.parser_name)
                archive = answer_entry_archive_request(
                    dict(upload_id=upload_id, mainfile=full_path),
                    required='*', user=user,
                    entry_metadata=entry_metadata)['data']['archive']

            response = PutRawFileResponse(
                upload_id=upload_id,
                data=upload_to_pydantic(upload),
                processing=ProcessingData(
                    upload_id=upload_id,
                    path=full_path,
                    entry_id=entry.entry_id if entry else None,
                    parser_name=entry.parser_name if entry else None,
                    entry=entry_to_pydantic(entry) if entry else None,
                    archive=archive))
            response_text = response.json()
            media_type = 'application/json'

        except HTTPException:
            raise
        except ProcessAlreadyRunning:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='The upload is currently being processed, operation not allowed.')
        except Exception as e:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f'Unexpected exception occurred: {e}')
        finally:
            try:
                shutil.rmtree(os.path.dirname(upload_path))
            except Exception:
                pass

    return StreamingResponse(create_stream_from_string(response_text), media_type=media_type)


@router.delete(
    '/{upload_id}/raw/{path:path}', tags=[raw_tag],
    summary='Delete the raw file or folder located at the specified path in the specified upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def delete_upload_raw_path(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload.'),
        path: str = Path(
            ...,
            description='The path within the upload raw files.'),
        user: User = Depends(create_user_dependency(required=True, upload_token_auth_allowed=True))):
    '''
    Delete file or folder located at the specified path in the specified upload. The upload
    must not be published. This also automatically triggers a reprocessing of the upload.
    Choosing the empty string as `path` deletes all files.
    '''
    upload = _get_upload_with_write_access(upload_id, user, include_published=False)

    if not is_safe_relative_path(path):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Bad path provided.')

    upload_files = StagingUploadFiles(upload_id)

    if not upload_files.raw_path_exists(path):
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail='No file or folder with that path found.')

    try:
        upload.process_upload(file_operations=[dict(op='DELETE', path=path)], only_updated_files=True)
    except ProcessAlreadyRunning:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='The upload is currently blocked by another process.')

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/raw-create-dir/{path:path}', tags=[raw_tag],
    summary='Create a new empty directory with the specified path in the specified upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_upload_raw_create_dir_path(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload.'),
        path: str = Path(
            ...,
            description='The path within the upload raw files.'),
        user: User = Depends(create_user_dependency(required=True, upload_token_auth_allowed=True))):
    '''
    Create a new empty directory in the specified upload. The `path` should be the full path
    to the new directory (i.e. ending with the name of the new directory). The api call returns
    immediately (no processing is neccessary).
    '''
    upload = _get_upload_with_write_access(upload_id, user, include_published=False)

    if not path or not is_safe_relative_path(path):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Bad path provided.')
    if upload.staging_upload_files.raw_path_exists(path):
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f'Path `{path}` already exists.')
    try:
        upload.staging_upload_files.raw_create_directory(path)
        # No real processing is needed when just adding a folder, but we should signal that
        # the upload has changed.
        upload.complete_time = datetime.utcnow()
        upload.save()
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f'Failed to create directory: {e}')

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.get(
    '/{upload_id}/archive/mainfile/{mainfile:path}', tags=[archive_tag],
    summary='Get the full archive for the given upload and mainfile path.',
    response_model=EntryArchiveResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_upload_or_path_not_found, _not_authorized_to_upload))
async def get_upload_entry_archive_mainfile(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload.'),
        mainfile: str = Path(
            ...,
            description='The mainfile path within the upload\'s raw files.'),
        mainfile_key: Optional[str] = FastApiQuery(
            None,
            description='The mainfile_key, for accessing child entries.'),
        user: User = Depends(create_user_dependency(required=False))):
    '''
    For the upload specified by `upload_id`, gets the full archive of a single entry that
    is identified by the given `mainfile`.
    '''
    get_upload_with_read_access(upload_id, user, include_others=True)
    query = dict(upload_id=upload_id, mainfile=mainfile)
    if mainfile_key:
        query.update(mainfile_key=mainfile_key)
    return answer_entry_archive_request(query, required='*', user=user)


@router.get(
    '/{upload_id}/archive/{entry_id}', tags=[archive_tag],
    summary='Get the full archive for the given upload and entry.',
    response_model=EntryArchiveResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_upload_or_path_not_found, _not_authorized_to_upload))
async def get_upload_entry_archive(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload.'),
        entry_id: str = Path(
            ...,
            description='The unique entry id.'),
        user: User = Depends(create_user_dependency(required=False))):
    '''
    For the upload specified by `upload_id`, gets the full archive of a single entry that
    is identified by the given `entry_id`.
    '''
    get_upload_with_read_access(upload_id, user, include_others=True)
    return answer_entry_archive_request(
        dict(upload_id=upload_id, entry_id=entry_id),
        required='*', user=user)


@router.post(
    '', tags=[default_tag],
    summary='Submit a new upload',
    response_class=StreamingResponse,
    responses=create_responses(_post_upload_response, _not_authorized, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_upload(
        request: Request,
        file: List[UploadFile] = File(None),
        local_path: str = FastApiQuery(
            None,
            description=strip('''
            Internal/Admin use only.''')),
        file_name: str = FastApiQuery(
            None,
            description=strip('''
            Specifies the name of the file, when using method 2.''')),
        upload_name: str = FastApiQuery(
            None,
            description=strip('''
            A human readable name for the upload.''')),
        embargo_length: int = FastApiQuery(
            0,
            description=strip('''
            The requested embargo length, in months, if any (0-36).''')),
        publish_directly: bool = FastApiQuery(
            None,
            description=strip('''
            If the upload should be published directly. False by default.''')),
        user: User = Depends(create_user_dependency(required=True, upload_token_auth_allowed=True))):
    '''
    Creates a new, empty upload and, optionally, uploads one or more files to it. If zip or
    tar files are uploaded, they will first be extracted, then added.

    It is recommended to give the upload itself a descriptive `upload_name`. If not specified,
    and a single file is provided, `upload_name` will be set to the name of this file. The
    `upload_name` can also be edited afterwards (as long as the upload is not published).

    There are two basic ways to upload files: in the multipart-formdata or streaming the
    file data in the http body. Both are supported. Note, however, that the second method
    only allows the upload of a single file, and that it does not transfer a filename. If a
    transfer is made using method 2, you can specify the query argument `file_name` to name it.
    This *needs* to be specified when using method 2, unless you are uploading a zip/tar file
    (for zip/tar files the names don't matter since they are extracted).

    Example curl commands for creating an upload and uploading a file:

    Method 1: multipart-formdata

        curl -X 'POST' "url" -F file=@local_file

    Method 2: streaming data

        curl -X 'POST' "url" -T local_file

    Authentication is required. This can either be done using the regular bearer token,
    or using the simplified upload token. To use the simplified upload token, just
    specify it as a query parameter in the url, i.e.

        curl -X 'POST' "baseurl?token=ABC.XYZ" ...

    Note, there is a limit on how many unpublished uploads a user can have. If exceeded,
    error code 400 will be returned.
    '''
    if not user.is_admin:
        # Check upload limit
        if _query_mongodb(main_author=str(user.user_id), publish_time=None).count() >= config.services.upload_limit:  # type: ignore
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                Limit of unpublished uploads exceeded for user.'''))

    if not 0 <= embargo_length <= 36:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='`embargo_length` must be between 0 and 36 months.')

    upload_id = utils.create_uuid()

    upload_paths, method = await _get_files_if_provided(
        upload_id, request, file, local_path, file_name, user)

    if not upload_name:
        # Try to default upload_name
        if method == 2:
            upload_name = file_name or None
        elif len(upload_paths) == 1:
            upload_name = os.path.basename(upload_paths[0])

    upload: Upload = Upload.create(
        upload_id=upload_id,
        main_author=user,
        upload_name=upload_name,
        upload_create_time=datetime.utcnow(),
        embargo_length=embargo_length,
        publish_directly=publish_directly)

    # Create staging files
    files.StagingUploadFiles(upload_id=upload_id, create=True)

    logger.info('upload created', upload_id=upload_id)

    if upload_paths:
        upload.process_upload(
            file_operations=[
                dict(op='ADD', path=upload_path, target_dir='', temporary=(method != 0))
                for upload_path in upload_paths])

    if request.headers.get('Accept') == 'application/json':
        upload_proc_data_response = UploadProcDataResponse(
            upload_id=upload_id,
            data=upload_to_pydantic(upload))
        response_text = upload_proc_data_response.json()
        media_type = 'application/json'
    else:
        response_text = _thank_you_message
        media_type = 'text/plain'

    return StreamingResponse(create_stream_from_string(response_text), media_type=media_type)


@router.post(
    '/{upload_id}/edit', tags=[metadata_tag],
    summary='Updates the metadata of the specified upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_upload_edit(
        request: Request,
        data: MetadataEditRequest,
        upload_id: str = Path(..., description='The unique id of the upload.'),
        user: User = Depends(create_user_dependency(required=True))):
    '''
    Updates the metadata of the specified upload and entries. An optional `query` can be
    specified to select only some of the entries of the upload (the query results are
    automatically restricted to the specified upload).

    **Note:**
      - Only admins can edit some of the fields.
      - The embargo of a published upload is lifted by setting the `embargo_length` attribute
        to 0.
      - If the upload is published, the only operations permitted using this endpoint is to
        lift the embargo, i.e. set `embargo_length` to 0, and to edit the entries in datasets
        that where created by the current user.
      - If a query is specified, it is not possible to edit upload level metadata (like
        `upload_name`, `coauthors`, etc.), as the purpose of queries is to select only a
        subset of the upload entries to edit, but changing upload level metadata would affect
        **all** entries of the upload.
    '''
    edit_request_json = await request.json()
    try:
        MetadataEditRequestHandler.edit_metadata(edit_request_json, upload_id, user)
        return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(Upload.get(upload_id)))
    except RequestValidationError:
        raise  # A problem which we have handled explicitly. Fastapi does json conversion.
    except Exception as e:
        # The upload is processing or some kind of unexpected error has occured
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.delete(
    '/{upload_id}', tags=[default_tag],
    summary='Delete an upload',
    response_model=UploadProcDataResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def delete_upload(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload to delete.'),
        user: User = Depends(create_user_dependency(required=True))):
    '''
    Delete an existing upload.

    Only uploads that are sill in staging, not already deleted, not still uploaded, and
    not currently processed, can be deleted.
    '''
    upload = _get_upload_with_write_access(
        upload_id, user, include_published=True, published_requires_admin=True, include_failed_imports=True)

    try:
        upload.delete_upload()
    except ProcessAlreadyRunning:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
            The upload is still being processed.'''))
    except Exception as e:
        logger.error('could not delete processing upload', exc_info=e)
        raise

    return UploadProcDataResponse(
        upload_id=upload_id,
        data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/publish', tags=[action_tag],
    summary='Publish an upload',
    response_model=UploadProcDataResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_upload_action_publish(
        upload_id: str = Path(
            ...,
            description=strip('''
                The unique id of the upload to publish.''')),
        embargo_length: int = FastApiQuery(
            None,
            description=strip('''
                If provided, updates the embargo length of the upload. The value should
                be between 0 and 36 months. 0 means no embargo.''')),
        to_central_nomad: bool = FastApiQuery(
            False,
            description=strip('''
                Will send the upload to the central NOMAD repository and publish it. This
                option is only available on an OASIS. The upload must already be published
                on the OASIS.''')),
        user: User = Depends(create_user_dependency(required=True))):
    '''
    Publishes an upload. The upload cannot be modified after this point (except for special
    cases, like when lifting the embargo prematurely, and by admins). After the upload is
    published and the embargo period (if any) is expired, the generated archive entries
    will be publicly visible.
    '''
    upload = _get_upload_with_write_access(
        upload_id, user, include_published=True, published_requires_admin=False)

    if upload.published and not user.is_admin and not to_central_nomad:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Upload already published.')

    _check_upload_not_processing(upload)

    if upload.process_status == ProcessStatus.FAILURE:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Cannot publish an upload that failed processing.')
    if upload.processed_entries_count == 0:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Cannot publish an upload without any resulting entries.')
    if embargo_length is not None and not 0 <= embargo_length <= 36:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Invalid embargo_length. Must be between 0 and 36 months.')

    if to_central_nomad:
        # Publish from an OASIS to the central repository
        if not config.oasis.is_oasis:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='Must be on an OASIS to publish to the central NOMAD repository.')
        if not upload.published:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail='The upload must be published on the OASIS first.')
        # Everything looks ok, try to publish it to the central NOMAD!
        upload.publish_externally(embargo_length=embargo_length)
    else:
        # Publish to this repository
        if upload.published:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail='The upload is already published.')
        try:
            upload.publish_upload(embargo_length=embargo_length)
        except ProcessAlreadyRunning:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='The upload is still/already processed.')

    return UploadProcDataResponse(
        upload_id=upload_id,
        data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/process', tags=[action_tag],
    summary='Manually triggers processing of an upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_upload_action_process(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload to process.'),
        user: User = Depends(create_user_dependency(required=True))):
    '''
    Processes an upload, i.e. parses the files and updates the NOMAD archive. Only admins
    can process an already published upload.
    '''
    upload = _get_upload_with_write_access(
        upload_id, user, include_published=True, published_requires_admin=True)

    _check_upload_not_processing(upload)

    upload.process_upload()
    return UploadProcDataResponse(
        upload_id=upload_id,
        data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/delete-entry-files', tags=[action_tag],
    summary='Deletes the files of the entries specified by a query.',
    response_model=UploadProcDataResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_upload_action_delete_entry_files(
        data: DeleteEntryFilesRequest,
        upload_id: str = Path(
            ...,
            description='The unique id of the upload within which to delete entry files.'),
        user: User = Depends(create_user_dependency(required=True))):
    ''' Deletes the files of the entries specified by the provided query. '''

    upload = _get_upload_with_write_access(upload_id, user, include_published=False)

    # Evaluate query
    if not data.query:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
            A query must be specified.'''))
    restricted_query = restrict_query_to_upload(data.query, upload_id)
    es_entries = search_iterator(
        user_id=user.user_id,
        owner=data.owner,
        query=restricted_query,
        required=MetadataRequired(include=['mainfile']))

    # Determine paths to delete
    try:
        paths_to_delete: Set[str] = set()
        for es_entry in es_entries:
            mainfile = es_entry['mainfile']
            path_to_delete = os.path.dirname(mainfile) if data.include_parent_folders else mainfile
            paths_to_delete.add(path_to_delete)
    except QueryValidationError as e:
        raise RequestValidationError(errors=e.errors)

    # Execute operation
    if paths_to_delete:
        try:
            upload.process_upload(
                file_operations=[
                    dict(op='DELETE', path=path_to_delete)
                    for path_to_delete in sorted(paths_to_delete)],
                only_updated_files=True)
        except ProcessAlreadyRunning:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='The upload is currently blocked by another process.')

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/lift-embargo', tags=[action_tag],
    summary='Lifts the embargo of an upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_upload_action_lift_embargo(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload to lift the embargo for.'),
        user: User = Depends(create_user_dependency(required=True))):
    ''' Lifts the embargo of an upload. '''
    upload = _get_upload_with_write_access(
        upload_id, user, include_published=True, published_requires_admin=False)
    _check_upload_not_processing(upload)
    if not upload.published:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
            Upload is not published, no embargo to lift.'''))
    if not upload.with_embargo:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
            Upload has no embargo.'''))
    # Lift the embargo using MetadataEditRequestHandler.edit_metadata
    try:
        MetadataEditRequestHandler.edit_metadata({'metadata': {'embargo_length': 0}}, upload_id, user)
        upload.reload()
        return UploadProcDataResponse(
            upload_id=upload_id,
            data=upload_to_pydantic(upload))
    except Exception as e:
        # Should only happen if the upload just started processing or something unexpected happens
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.get(
    '/{upload_id}/bundle', tags=[bundle_tag],
    summary='Gets an *upload bundle* for the specified upload.',
    response_class=StreamingResponse,
    responses=create_responses(
        _upload_bundle_response, _upload_not_found, _not_authorized_to_upload, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_upload_bundle(
        upload_id: str = Path(
        ...,
        description='The unique id of the upload.'),
        include_raw_files: Optional[bool] = FastApiQuery(
            True,
            description=strip('''
                If raw files should be included in the bundle (true by default).''')),
        include_archive_files: Optional[bool] = FastApiQuery(
            True,
            description=strip('''
                If archive files (i.e. parsed entries data) should be included in the bundle
                (true by default).''')),
        include_datasets: Optional[bool] = FastApiQuery(
            True,
            description=strip('''
                If datasets references to this upload should be included in the bundle
                (true by default).''')),
        user: User = Depends(create_user_dependency(required=False))):
    '''
    Get an *upload bundle* for the specified upload. An upload bundle is a file bundle which
    can be used to export and import uploads between different NOMAD deployments.
    '''
    upload = get_upload_with_read_access(upload_id, user, include_others=True)
    _check_upload_not_processing(upload)

    export_settings = config.bundle_export.default_settings.customize(
        dict(
            include_raw_files=include_raw_files,
            include_archive_files=include_archive_files,
            include_datasets=include_datasets))

    try:
        stream = BundleExporter(
            upload,
            export_as_stream=True,
            export_path=None,
            zipped=True,
            overwrite=False,
            export_settings=export_settings).export_bundle()
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
            Could not export due to error: ''' + str(e)))

    return StreamingResponse(iter(stream), media_type='application/zip')


@router.post(
    '/bundle', tags=[bundle_tag],
    summary='Posts an *upload bundle* to this NOMAD deployment.',
    response_model=UploadProcDataResponse,
    responses=create_responses(_not_authorized, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_upload_bundle(
        request: Request,
        file: List[UploadFile] = File(None),
        local_path: str = FastApiQuery(
            None,
            description=strip('''
            Internal/Admin use only.''')),
        embargo_length: Optional[int] = FastApiQuery(
            None,
            description=strip('''
                Specifies the embargo length in months to set on the upload. If omitted,
                the value specified in the bundle will be used. A value of 0 means no
                embargo.''')),
        include_raw_files: Optional[bool] = FastApiQuery(
            None,
            description=strip('''
                If raw files should be imported from the bundle
                *(only admins can change this setting)*.''')),
        include_archive_files: Optional[bool] = FastApiQuery(
            None,
            description=strip('''
                If archive files (i.e. parsed entries data) should be imported from the bundle
                *(only admins can change this setting)*.''')),
        include_datasets: Optional[bool] = FastApiQuery(
            None,
            description=strip('''
                If dataset references to this upload should be imported from the bundle
                *(only admins can change this setting)*.''')),
        include_bundle_info: Optional[bool] = FastApiQuery(
            None,
            description=strip('''
                If the bundle_info.json file should be kept
                *(only admins can change this setting)*.''')),
        keep_original_timestamps: Optional[bool] = FastApiQuery(
            None,
            description=strip('''
                If all original timestamps, including `upload_create_time`, `entry_create_time`
                and `publish_time`, should be kept
                *(only admins can change this setting)*.''')),
        set_from_oasis: Optional[bool] = FastApiQuery(
            None,
            description=strip('''
                If the `from_oasis` flag and `oasis_deployment_url` should be set
                *(only admins can change this setting)*.''')),
        trigger_processing: Optional[bool] = FastApiQuery(
            None,
            description=strip('''
                If processing should be triggered after the bundle has been imported
                *(only admins can change this setting)*.''')),
        user: User = Depends(create_user_dependency(required=True, upload_token_auth_allowed=True))):
    '''
    Posts an *upload bundle* to this NOMAD deployment. An upload bundle is a file bundle which
    can be used to export and import uploads between different NOMAD installations. The
    endpoint expects an upload bundle attached as a zipfile.

    **NOTE:** This endpoint is restricted to admin users and oasis admins. Further, all
    settings except `embargo_length` requires an admin user to change (these settings
    have default values specified by the system configuration).

    There are two basic ways to upload files: in the multipart-formdata or streaming the
    file data in the http body. Both are supported. See the POST `uploads` endpoint for
    examples of curl commands for uploading files.
    '''
    import_settings = config.bundle_import.default_settings.customize(
        dict(
            include_raw_files=include_raw_files,
            include_archive_files=include_archive_files,
            include_datasets=include_datasets,
            include_bundle_info=include_bundle_info,
            keep_original_timestamps=keep_original_timestamps,
            set_from_oasis=set_from_oasis,
            trigger_processing=trigger_processing))

    bundle_importer: BundleImporter = None
    bundle_path: str = None
    method = None

    try:
        bundle_importer = BundleImporter(user, import_settings)
        bundle_importer.check_api_permissions()

        bundle_paths, method = await _get_files_if_provided(
            tmp_dir_prefix='bundle', request=request, file=file, local_path=local_path, file_name=None, user=user)

        if not bundle_paths:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST, detail='No bundle file provided')
        if len(bundle_paths) > 1:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST, detail='Can only provide one bundle file at a time')
        bundle_path = bundle_paths[0]

        bundle_importer.open(bundle_path)
        upload = bundle_importer.create_upload_skeleton()
        bundle_importer.close()
        # Run the import as a @process
        upload.import_bundle(
            bundle_path=bundle_path,
            import_settings=import_settings.dict(),
            embargo_length=embargo_length)
        return UploadProcDataResponse(
            upload_id=upload.upload_id,
            data=upload_to_pydantic(upload))

    except Exception as e:
        if bundle_importer:
            bundle_importer.close()
            if bundle_path and method != 0:
                bundle_importer.delete_bundle()
        if isinstance(e, HTTPException):
            raise
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail='Failed to import bundle: ' + str(e))


async def _get_files_if_provided(
        tmp_dir_prefix: str, request: Request, file: List[UploadFile], local_path: str, file_name: str,
        user: User) -> Tuple[List[str], Union[None, int]]:
    '''
    If the user provides one or more files with the api call, load and save them to a temporary
    folder (or, if method 0 is used, just "forward" the file path). The method thus needs to identify
    which file transfer method was used (0 - 2), and save the data to disk (if method is 1 or 2).

    Returns a list of os paths to the resulting files and method (0-2), or ([], None) if no file
    data was provided with the api call.
    '''
    # Determine the source data stream
    sources: List[Tuple[Any, str]] = []  # List of tuples (source, filename)
    if local_path:
        # Method 0: Local file - only for local path under /tests/data/proc/
        if not local_path.startswith('examples') and not user.is_admin:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
                You are not authorized to access this path.'''))
        if not os.path.exists(local_path) or not os.path.isfile(local_path):
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                The specified local_path cannot be found or is not a file.'''))
        method = 0
    elif file:
        # Method 1: Data provided as formdata
        method = 1
        sources = [
            (_asyncronous_file_reader(multipart_file), unquote(multipart_file.filename))
            for multipart_file in file]
    else:
        # Method 2: Data has to be sent streamed in the body
        method = 2
        sources = [(request.stream(), file_name or 'NO NAME')]

    no_file_name_info_provided = not file_name

    for _, file_name in sources:
        if not files.is_safe_basename(file_name):
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                Bad file name provided.'''))

    # Forward the file path (if method == 0) or save the file(s)
    if method == 0:
        upload_paths = [local_path]  # Use the provided path
        uploaded_bytes = os.path.getsize(local_path)
    else:
        tmp_dir = files.create_tmp_dir(tmp_dir_prefix)
        upload_paths = []
        uploaded_bytes = 0
        for source_stream, file_name in sources:
            upload_path = os.path.join(tmp_dir, file_name)
            try:
                with open(upload_path, 'wb') as f:
                    uploaded_bytes = 0
                    log_interval = 1e9
                    log_unit = 'GB'
                    next_log_at = log_interval
                    async for chunk in source_stream:
                        if not chunk:
                            # End of data stream
                            break
                        uploaded_bytes += len(chunk)
                        f.write(chunk)
                        if uploaded_bytes > next_log_at:
                            logger.info('Large upload in progress - uploaded: '
                                        f'{uploaded_bytes // log_interval} {log_unit}')
                            next_log_at += log_interval
                    logger.info(f'Uploaded {uploaded_bytes} bytes')
            except Exception as e:
                if not (isinstance(e, RuntimeError) and 'Stream consumed' in str(e)):
                    if os.path.exists(tmp_dir):
                        shutil.rmtree(tmp_dir)
                    logger.warn('IO error receiving upload data', exc_info=e)
                    raise HTTPException(
                        status_code=status.HTTP_400_BAD_REQUEST,
                        detail='Some IO went wrong, upload probably aborted/disrupted.')
            upload_paths.append(upload_path)

        if not uploaded_bytes and method == 2:
            # No data was provided
            shutil.rmtree(tmp_dir)
            return [], None

    logger.info(f'received {len(upload_paths)} uploaded file(s)')
    if method == 2 and no_file_name_info_provided:
        # Only ok if uploaded file is a zip or a tar archive.
        ext = '.zip' if files.zipfile.is_zipfile(upload_path) else '.tar' if files.tarfile.is_tarfile(upload_path) else None
        if not ext:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                No file name provided, and the file does not look like a zip or tar file.'''))
        # Add the correct extension
        shutil.move(upload_path, upload_path + ext)
        upload_paths = [upload_path + ext]

    return upload_paths, method


async def _asyncronous_file_reader(f):
    ''' Asynchronous generator to read file-like objects. '''
    while True:
        try:
            data: bytes = await f.read(io.DEFAULT_BUFFER_SIZE)
        except Exception:
            await f.close()
            raise
        if not data:
            await f.close()
            return
        yield data


def _query_mongodb(**kwargs):
    return Upload.objects(**kwargs)


def get_upload_with_read_access(upload_id: str, user: User, include_others: bool = False) -> Upload:
    '''
    Determines if the specified user has read access to the specified upload. If so, the
    corresponding Upload object is returned. If the upload does not exist, or the user has
    no read access to it, a HTTPException is raised.

    Arguments:
        upload_id: The id of the requested upload.
        user: The authenticated user, if any.
        include_others: If uploads owned by others should be included. Access to the uploads
            of other users is only granted if the upload is published and not under embargo.
    '''
    mongodb_query = _query_mongodb(upload_id=upload_id)
    upload = mongodb_query.first()
    if upload is None:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=strip('''
            The specified upload_id was not found.'''))
    if user and (user.is_admin or (str(user.user_id) in upload.viewers)):
        # Ok, the user a viewer, or we have an admin user
        return upload
    elif include_others:
        if not upload.published:
            if user:
                raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
                    You do not have access to the specified upload.'''))
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
                You need to log in to access the specified upload.'''))
        if upload.published and upload.with_embargo:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
                You do not have access to the specified upload - published with embargo.'''))
        return upload
    else:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            You do not have access to the specified upload.'''))


def _get_upload_with_write_access(
        upload_id: str, user: User, include_published: bool = False,
        published_requires_admin: bool = True, include_failed_imports: bool = False) -> Upload:
    '''
    Determines if the specified user has write access to the specified upload. If so, the
    corresponding Upload object is returned. If the upload does not exist, or the user has
    no write access to it, a HTTPException is raised.
    '''
    if not user:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            User authentication required to access uploads.'''))
    mongodb_query = _query_mongodb(upload_id=upload_id)
    if not mongodb_query.count():
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=strip('''
            The specified upload_id was not found.'''))
    upload = mongodb_query.first()
    if not user.is_admin and upload.main_author != str(user.user_id):
        if not upload.coauthors or str(user.user_id) not in upload.coauthors:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
                You do not have write access to the specified upload.'''))
    if upload.published:
        is_failed_import = (
            upload.current_process.startswith('import_bundle') and upload.process_status == ProcessStatus.FAILURE)
        if not include_published:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
                Upload is already published, operation not possible.'''))
        if published_requires_admin and not user.is_admin and not (is_failed_import and include_failed_imports):
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
                Upload is already published, only admins can perform this operation.'''))
    return upload


def upload_to_pydantic(upload: Upload) -> UploadProcData:
    ''' Converts the mongo db object to an UploadProcData object. '''
    pydantic_upload = UploadProcData.from_orm(upload)
    pydantic_upload.entries = upload.total_entries_count
    try:
        pydantic_upload.upload_files_server_path = upload.upload_files.external_os_path
    except KeyError:
        # In case the files are missing for one reason or another
        pass

    return pydantic_upload


def entry_to_pydantic(entry: Entry, add_es_metadata: bool = False, user=None) -> EntryProcData:
    '''
    Converts the mongo db object to an EntryProcData object, and optionally also adds metadata
    from ES
    '''
    rv = EntryProcData.from_orm(entry)
    if add_es_metadata:
        # load entries's metadata from search
        metadata_entries = search(
            pagination=MetadataPagination(page_size=1),
            owner='admin' if user.is_admin else 'visible',
            user_id=user.user_id,
            query=dict(entry_id=entry.entry_id))
        if len(metadata_entries.data) == 1:
            rv.entry_metadata = metadata_entries.data[0]
    return rv


def _check_upload_not_processing(upload: Upload):
    '''
    Checks if the upload is processing, and raises a HTTPException (err code 400) if so.
    '''
    if upload.process_running:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='The upload is currently being processed, operation not allowed.')
