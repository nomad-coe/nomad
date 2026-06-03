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
import io
import os
import shutil
import tarfile
import zipfile
from asyncio import sleep
from datetime import datetime, timezone
from enum import Enum
from typing import Annotated, Any, cast
from urllib.parse import unquote, urlparse

import requests
from fastapi import (
    APIRouter,
    Body,
    Depends,
    File,
    HTTPException,
    Path,
    Request,
    UploadFile,
    status,
)
from fastapi import Query as FastApiQuery
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, StreamingResponse
from mongoengine.errors import MongoEngineException
from mongoengine.queryset.visitor import Q
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator
from pydantic_core import PydanticCustomError

from nomad import files, utils
from nomad.app.v1.models.models import TransferBundleRequest
from nomad.auth.scopes import Scope
from nomad.auth.tokens import generate_upload_token
from nomad.bundles import BundleExporter, BundleImporter
from nomad.common import get_compression_format, is_safe_basename, is_safe_relative_path
from nomad.config import config
from nomad.config.models.config import Reprocess
from nomad.config.models.plugins import ExampleUploadEntryPoint
from nomad.datacite import DataCiteException
from nomad.datacite.service import create_doi_for_upload, publish_doi
from nomad.files import FSUtility, PublicUploadFiles, StagingUploadFiles
from nomad.models.common import UTCDateTime
from nomad.mongo.doi import EmbeddedDOI
from nomad.mongo.groups import MongoUserGroup
from nomad.mongo.search import MongoQueryError, create_mongo_query
from nomad.processing import (
    Entry,
    MetadataEditRequestHandler,
    ProcessAlreadyRunning,
    ProcessStatus,
    Upload,
)
from nomad.search import QueryValidationError, search, search_iterator
from nomad.search import refresh as search_refresh
from nomad.utils import strip

from ..models import (
    Direction,
    Files,
    HTTPExceptionModel,
    MetadataEditRequest,
    MetadataPagination,
    MetadataRequired,
    Owner,
    Pagination,
    PaginationResponse,
    User,
    WithQuery,
    files_parameters,
    restrict_query_to_upload,
)
from ..utils import (
    DownloadItem,
    browser_download_headers,
    create_download_stream_raw_file,
    create_download_stream_zipped,
    create_responses,
    create_stream_from_string,
    parameter_dependency_from_model,
)
from .auth import get_current_user
from .entries import EntryArchiveResponse, answer_entry_archive_request

router = APIRouter()


class APITag(str, Enum):
    DEFAULT = 'uploads'
    METADATA = 'uploads/metadata'
    RAW = 'uploads/raw'
    ARCHIVE = 'uploads/archive'
    ACTION = 'uploads/action'
    BUNDLE = 'uploads/bundle'


logger = utils.get_logger(__name__)


async def async_wrapper(content):
    for x in content:
        yield x


class UploadRole(str, Enum):
    main_author = 'main_author'
    reviewer = 'reviewer'
    coauthor = 'coauthor'


class DOI(BaseModel):
    model_config = ConfigDict(from_attributes=True)

    id: str = Field(description='The DOI name, e.g. 10.2345/nomad.6789-wxyz')


class ProcData(BaseModel):
    process_running: bool = Field(description='If a process is running')
    current_process: str | None = Field(
        None, description='Name of the current or last completed process'
    )
    process_status: str = Field(
        ProcessStatus.READY,
        description='The status of the current or last completed process',
    )
    last_status_message: str | None = Field(
        None,
        description='A short, human readable message from the current process, with '
        'information about what the current process is doing, or information '
        'about the completion (successful or not) of the last process, if no '
        'process is currently running.',
    )
    errors: list[str] = Field(
        description='A list of error messages that occurred during the last processing'
    )
    warnings: list[str] = Field(
        description='A list of warning messages that occurred during the last processing'
    )
    complete_time: UTCDateTime | None = Field(
        None, description='Date and time of the completion of the last process'
    )
    model_config = ConfigDict(from_attributes=True)


class UploadProcData(ProcData):
    upload_id: str = Field(description='The unique id for the upload.')
    upload_name: str | None = Field(
        None,
        description='The name of the upload. This can be provided during upload '
        'using the `upload_name` query parameter.',
    )
    upload_create_time: UTCDateTime | None = Field(
        None, description='Date and time of the creation of the upload.'
    )
    description: str | None = Field(
        None, description='Further information about the upload.'
    )
    main_author: str | None = Field(
        None, description=strip('The main author of the upload.')
    )
    coauthors: list[str] | None = Field(
        None, description=strip('A list of upload coauthors.')
    )
    coauthor_groups: list[str] | None = Field(
        None, description=strip('A list of upload coauthor groups.')
    )
    reviewers: list[str] | None = Field(
        None, description=strip('A list of upload reviewers.')
    )
    reviewer_groups: list[str] | None = Field(
        None, description=strip('A list of upload reviewer groups.')
    )
    writers: list[str] | None = Field(
        None, description=strip('All writer users (main author, upload coauthors).')
    )
    writer_groups: list[str] | None = Field(
        None, description=strip('All writer groups (coauthor groups).')
    )
    viewers: list[str] | None = Field(
        None,
        description=strip(
            'All viewer users (main author, upload coauthors, and reviewers)'
        ),
    )
    viewer_groups: list[str] | None = Field(
        None,
        description=strip('All viewer groups (coauthor groups, reviewer groups).'),
    )
    published: bool = Field(False, description='If this upload is already published.')
    published_to: list[str] | None = Field(
        None,
        description='A list of other NOMAD deployments that this upload was uploaded to already.',
    )
    publish_time: UTCDateTime | None = Field(
        None,
        description='Date and time of publication, if the upload has been published.',
    )
    with_embargo: bool = Field(
        description='If the upload has an embargo set (embargo_length not equal to zero).'
    )
    embargo_length: int = Field(
        description='The length of the requested embargo, in months. 0 if no embargo is requested.'
    )
    license: str = Field(
        description='The license under which this upload is distributed.'
    )
    doi: DOI | None = Field(None, description='The DOI assigned to this upload.')
    entries: int = Field(
        0, description='The number of identified entries in this upload.'
    )
    upload_files_server_path: str | None = Field(
        None, description='The path to the uploads files on the server.'
    )


class EntryProcData(ProcData):
    entry_id: str = Field()
    entry_create_time: UTCDateTime = Field()
    mainfile: str = Field()
    mainfile_key: str | None = Field(None)
    upload_id: str = Field()
    parser_name: str = Field()
    entry_metadata: dict | None = Field(None)


class UploadProcDataPagination(Pagination):
    @model_validator(mode='before')
    @classmethod
    def check_order_by(cls, data):
        order_by_choices = (
            'upload_create_time',
            'publish_time',
            'upload_name',
            'last_status_message',
            'process_status',
        )
        if isinstance(data, dict):
            order_by = data.get('order_by')
            if order_by is None:
                order_by = 'upload_create_time'  # Default value
                data['order_by'] = order_by
            if order_by not in order_by_choices:
                raise PydanticCustomError(
                    'invalid_order_by',
                    f"order_by is '{order_by}', must be one of {order_by_choices}",
                )
        return data

    @field_validator('page_after_value')
    @classmethod
    def validate_page_after_value(cls, page_after_value, values):
        # Validation handled elsewhere
        return page_after_value

    @field_validator('order_by')
    @classmethod
    def validate_order_by(cls, order_by, values):
        # Validation handled elsewhere
        return order_by

    def order_result(self, result):
        if self.order_by is None:
            return result

        prefix: str = '-' if self.order == Direction.desc else '+'
        order_list: list = [f'{prefix}{self.order_by}']
        if self.order_by == 'upload_create_time':
            order_list.append('upload_id')
        else:
            order_list.extend(['upload_create_time', 'upload_id'])

        return result.order_by(*order_list)


upload_proc_data_pagination_parameters = parameter_dependency_from_model(
    'upload_proc_data_pagination_parameters',
    UploadProcDataPagination,
)


class EntryProcDataPagination(Pagination):
    @field_validator('order_by')
    @classmethod
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        order_by_choices = (
            'mainfile',
            'parser_name',
            'process_status',
            'current_process',
            'entry_create_time',
        )

        if order_by == 'mainfile_path':
            return 'mainfile'
        if order_by is None:
            return 'mainfile'  # Default value
        if order_by not in order_by_choices:
            raise PydanticCustomError(
                'invalid_order_by',
                f"order_by is '{order_by}', must be one of {order_by_choices}",
            )
        return order_by

    @field_validator('page_after_value')
    @classmethod
    def validate_page_after_value(cls, page_after_value, values):
        # Validation handled elsewhere
        return page_after_value

    def order_result(self, result):
        if self.order_by is None:
            return result

        prefix: str = '-' if self.order == Direction.desc else '+'
        order_list: list = [f'{prefix}{self.order_by}', 'entry_id']

        return result.order_by(*order_list)


entry_proc_data_pagination_parameters = parameter_dependency_from_model(
    'entry_proc_data_pagination_parameters',
    EntryProcDataPagination,
)


class UploadProcDataResponse(BaseModel):
    upload_id: str | None = Field(
        None,
        description=strip(
            """
        Unique id of the upload."""
        ),
    )
    data: UploadProcData | None = Field(
        None,
        description=strip(
            """
        The upload data as a dictionary."""
        ),
    )


class UploadProcDataQuery(BaseModel):
    upload_id: list[str] | None = Field(
        None,
        description='Search for uploads matching the given id. Multiple values can be specified.',
    )
    doi: list[str] | None = Field(
        None,
        description='Search for uploads matching the given doi. Multiple values can be specified.',
    )
    upload_name: list[str] | None = Field(
        None,
        description=strip(
            """
            Search for uploads by upload_name.

            Implicit exact form:
            - `upload_name=value`

            Explicit form:
            - `upload_name={"value":"name","type":"exact"}`
            - `upload_name={"value":"term","type":"fuzzy"}`

            Multiple values can be specified.
            """
        ),
    )
    is_processing: bool | None = Field(
        None,
        description=strip(
            """
            If True, only include currently processing uploads.
            If False, do not include currently processing uploads.
            If unset, include everything."""
        ),
    )
    is_published: bool | None = Field(
        None,
        description=strip(
            """
            If True: only include published uploads.
            If False: only include unpublished uploads.
            If unset: include everything."""
        ),
    )
    process_status: str | None = Field(
        None, description=strip('Search by the process status.')
    )
    is_owned: bool | None = Field(
        None,
        description=strip(
            """
            If True: only include owned uploads.
            If False: only include shared uploads.
            If unset: include everything."""
        ),
    )

    @field_validator('process_status')
    @classmethod
    def upper_process_status(cls, process_status: str):  # pylint: disable=no-self-argument
        return process_status.upper() if process_status else None


upload_proc_data_query_parameters = parameter_dependency_from_model(
    'upload_proc_data_query_parameters',
    UploadProcDataQuery,
)


class UploadProcDataQueryResponse(BaseModel):
    query: UploadProcDataQuery = Field()
    pagination: PaginationResponse = Field()
    data: list[UploadProcData] | None = Field(
        None,
        description=strip(
            """
        The upload data as a list. Each item is a dictionary with the data for each
        upload."""
        ),
    )


class EntryProcDataResponse(BaseModel):
    entry_id: str = Field()
    data: EntryProcData = Field()


class EntryProcDataQueryResponse(BaseModel):
    pagination: PaginationResponse = Field()
    processing_successful: int | None = Field(
        None,
        description=strip(
            """
        Number of entries that has been processed successfully.
        """
        ),
    )
    processing_failed: int | None = Field(
        None,
        description=strip(
            """
        Number of entries that failed to process.
        """
        ),
    )
    upload: UploadProcData | None = Field(
        None,
        description=strip(
            """
        The upload processing data of the upload.
        """
        ),
    )
    data: list[EntryProcData] | None = Field(
        None,
        description=strip(
            """
        The entries data as a list. Each item is a dictionary with the data for one entry.
        """
        ),
    )


class RawDirPagination(Pagination):
    @field_validator('order_by')
    @classmethod
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        assert not order_by, 'Cannot specify `order_by` for rawdir calls'
        if order_by:
            raise PydanticCustomError(
                'invalid_order_by', 'Cannot specify `order_by` for rawdir calls'
            )

    @field_validator('page_after_value')
    @classmethod
    def validate_page_after_value(cls, page_after_value, values):
        # Validation handled elsewhere
        return page_after_value


rawdir_pagination_parameters = parameter_dependency_from_model(
    'rawdir_pagination_parameters',
    RawDirPagination,
    exclude=['order', 'order_by'],
)


class RawDirFileMetadata(BaseModel):
    """Metadata about a file"""

    name: str = Field()
    size: int | None = Field(None)
    entry_id: str | None = Field(
        None,
        description=strip(
            """
        If this is a mainfile: the ID of the corresponding entry."""
        ),
    )
    parser_name: str | None = Field(
        None,
        description=strip(
            """
        If this is a mainfile: the name of the matched parser."""
        ),
    )


class RawDirElementMetadata(RawDirFileMetadata):
    """Metadata about an directory *element*, i.e. a file or a directory"""

    is_file: bool = Field()


class RawDirDirectoryMetadata(BaseModel):
    """Metadata about a directory"""

    name: str = Field()
    size: int | None = Field(None)
    content: list[RawDirElementMetadata] = Field(
        examples=[
            [
                {'name': 'a_directory', 'is_file': False, 'size': 456},
                {
                    'name': 'a_file.json',
                    'is_file': True,
                    'size': 123,
                    'entry_id': 'XYZ',
                    'parser_name': 'parsers/vasp',
                },
            ]
        ]
    )


class RawDirResponse(BaseModel):
    path: str = Field(examples=['The/requested/path'])
    access: str = Field()
    file_metadata: RawDirFileMetadata | None = Field(None)
    directory_metadata: RawDirDirectoryMetadata | None = Field(None)
    pagination: PaginationResponse | None = Field(None)


class ProcessingData(BaseModel):
    upload_id: str = Field()
    path: str = Field()
    entry_id: str | None = Field(None)
    parser_name: str | None = Field(None)
    entry: EntryProcData | None = Field(None)
    archive: dict[str, Any] | None = Field(None)


class PutRawFileResponse(BaseModel):
    upload_id: str | None = Field(
        None,
        description=strip(
            """
        Unique id of the upload."""
        ),
    )
    data: UploadProcData | None = Field(
        None,
        description=strip(
            """
        The upload data as a dictionary."""
        ),
    )
    processing: ProcessingData | None = Field(
        None,
        description=strip(
            """
        Information about the processing, including the entry (if one was generated) and
        [optionally] the archive data of this entry."""
        ),
    )


class DeleteEntryFilesRequest(WithQuery):
    """Defines a request to delete entry files."""

    owner: Owner | None = Body('all')
    include_parent_folders: bool | None = Field(
        False,
        description=strip(
            """
            If the delete operation should include not only the mainfiles of the selected entries,
            but also their folders."""
        ),
    )


class UploadCommandExamplesResponse(BaseModel):
    upload_url: str = Field()
    upload_command: str = Field()
    upload_command_with_name: str = Field()
    upload_progress_command: str = Field()
    upload_command_form: str = Field()
    upload_tar_command: str = Field()


_not_authorized = (
    status.HTTP_401_UNAUTHORIZED,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        Unauthorized. Authorization is required, but no or bad authentication credentials provided."""
        ),
    },
)

_not_authorized_to_upload = (
    status.HTTP_401_UNAUTHORIZED,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        Unauthorized. No credentials provided, or you do not have permissions to the
        specified upload."""
        ),
    },
)

_not_authorized_to_entry = (
    status.HTTP_401_UNAUTHORIZED,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        Unauthorized. No credentials provided, or you do not have permissions to the
        specified upload or entry."""
        ),
    },
)

_bad_request = (
    status.HTTP_400_BAD_REQUEST,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        Bad request. The request could not be processed because of some error/invalid argument."""
        ),
    },
)

_bad_pagination = (
    status.HTTP_400_BAD_REQUEST,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        Bad request. Invalid pagination arguments supplied."""
        ),
    },
)

_upload_not_found = (
    status.HTTP_404_NOT_FOUND,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        The specified upload could not be found."""
        ),
    },
)

_entry_not_found = (
    status.HTTP_404_NOT_FOUND,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        The specified upload or entry could not be found."""
        ),
    },
)

_upload_or_path_not_found = (
    status.HTTP_404_NOT_FOUND,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        The specified upload, or a resource with the specified path within the upload,
        could not be found."""
        ),
    },
)

_post_upload_response = (
    200,
    {
        'model': UploadProcDataResponse,
        'content': {
            'application/json': {},
            'text/plain': {'example': 'Thanks for uploading your data to nomad.'},
        },
        'description': strip(
            """
        A json structure with upload data, or a plain text information string.
        It will be a json structure if the request headers specifies `Accept = application/json`."""
        ),
    },
)

_put_raw_file_response = (
    200,
    {
        'model': PutRawFileResponse,
        'content': {
            'application/json': {},
            'text/plain': {'example': 'Thanks for uploading your data to nomad.'},
        },
        'description': strip(
            """
        A json structure with upload data and possibly information from the processing,
        or a plain text information string.
        It will be a json structure if the request headers specifies `Accept = application/json`
        or if `wait_for_processing` is set."""
        ),
    },
)


_raw_path_response = (
    200,
    {
        'content': {
            'application/octet-stream': {'example': 'file data'},
            'application/zip': {'example': '<zipped file or directory content>'},
        },
        'description': strip(
            """
        If `path` denotes a file: a stream with the file content, zipped if `compress = true`.
        If `path` denotes a directory, and `compress = true`, the directory content, zipped."""
        ),
    },
)

_upload_bundle_response = (
    200,
    {'content': {'application/zip': {'example': '<zipped bundle data>'}}},
)

_existing_upload_with_findable_state = (
    status.HTTP_400_BAD_REQUEST,
    {
        'model': HTTPExceptionModel,
        'description': 'The upload has failed to be submitted previously.',
    },
)

_upload_already_has_doi = (
    status.HTTP_400_BAD_REQUEST,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        The upload already has a DOI and cannot be changed anymore.
    """
        ),
    },
)

_upload_is_empty = (
    status.HTTP_400_BAD_REQUEST,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        The Upload is empty. No DOI can be assigned at this moment. Add some published
        contents to the upload first.
    """
        ),
    },
)

_upload_is_unpublished = (
    status.HTTP_400_BAD_REQUEST,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        The upload is unpublished. No DOI can be assigned at the moment.
        Publish the upload first.
    """
        ),
    },
)

_datacite_not_enabled = (
    status.HTTP_403_FORBIDDEN,
    {
        'model': HTTPExceptionModel,
        'description': 'The DataCite DOI service is not enabled on this deployment.',
    },
)

_datacite_draft_failed = (
    status.HTTP_500_INTERNAL_SERVER_ERROR,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        An error occurred while creating the DOI draft at DataCite. Please contact the administrator.
    """
        ),
    },
)

_db_set_upload_doi_failed = (
    status.HTTP_500_INTERNAL_SERVER_ERROR,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        An error occurred while saving the upload doi to the database. Please contact the administrator.
    """
        ),
    },
)

_datacite_publish_failed = (
    status.HTTP_500_INTERNAL_SERVER_ERROR,
    {
        'model': HTTPExceptionModel,
        'description': strip(
            """
        An error occurred while publishing the DOI at DataCite. Please contact the administrator.
    """
        ),
    },
)

_thank_you_message = f"""
Thanks for uploading your data to nomad.
Go back to {config.gui_url()} and press
reload to see the progress on your upload
and publish your data."""


def _create_exception(status_code, response_dict):
    return HTTPException(status_code, detail=response_dict.get('description'))


@router.get(
    '/command-examples',
    tags=[APITag.DEFAULT],
    summary='Get example commands for shell based uploads.',
    response_model=UploadCommandExamplesResponse,
    responses=create_responses(_not_authorized),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def get_command_examples(
    user: Annotated[
        User,
        Depends(get_current_user([Scope.TOKENS_CREATE], allow_anonymous=False)),
    ],
):
    """Get URL and example command for shell based uploads."""
    token = generate_upload_token(user)
    api_url = config.api_url(ssl=config.services.https_upload, api='api/v1')
    upload_url = f'{api_url}/uploads'
    header_flag = f"-H 'Upload-Token: {token}'"

    upload_command = f"curl -X POST {header_flag} '{upload_url}' -T <local_file>"
    upload_command_form = (
        f"curl -X POST {header_flag} '{upload_url}' -F file=@<local_file>"
    )
    upload_command_with_name = (
        f"curl -X POST {header_flag} '{upload_url}?upload_name=<name>' -T <local_file>"
    )
    upload_progress_command = upload_command + ' | xargs echo'
    upload_tar_command = f"tar -cf - <local_folder> | curl -# {header_flag} '{upload_url}' -X POST -T - | xargs echo"

    return UploadCommandExamplesResponse(
        upload_url=upload_url,
        upload_command=upload_command,
        upload_command_form=upload_command_form,
        upload_command_with_name=upload_command_with_name,
        upload_progress_command=upload_progress_command,
        upload_tar_command=upload_tar_command,
    )


@router.get(
    '',
    tags=[APITag.METADATA],
    summary='List uploads of authenticated user.',
    response_model=UploadProcDataQueryResponse,
    responses=create_responses(_not_authorized, _bad_pagination, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def get_uploads(
    request: Request,
    query: Annotated[UploadProcDataQuery, Depends(upload_proc_data_query_parameters)],
    pagination: Annotated[
        UploadProcDataPagination, Depends(upload_proc_data_pagination_parameters)
    ],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_READ], allow_anonymous=False)),
    ],
    roles: Annotated[
        list[UploadRole] | None,
        FastApiQuery(
            description='Only return uploads where the user has one of the given roles.'
        ),
    ] = None,
    include_all: Annotated[
        bool,
        FastApiQuery(description='Include uploads that are shared with all users.'),
    ] = False,
):
    """
    Retrieves metadata about all uploads that match the given query criteria.
    """
    # Build query
    role_query = get_role_query(roles, user, include_all=include_all)
    try:
        mongo_query = create_mongo_query(
            query,
            base_query=role_query,
            auth_user_id=str(user.user_id),
        )
    except MongoQueryError as e:
        raise HTTPException(status.HTTP_400_BAD_REQUEST, detail=str(e))

    # Create response
    start = pagination.get_simple_index()
    end = start + pagination.page_size

    # Fetch data from DB
    mongodb_query = pagination.order_result(Upload.objects.filter(mongo_query))  # type: ignore
    uploads = list(mongodb_query[start:end])
    upload_ids = [upload.upload_id for upload in uploads]

    # Batch fetch entry counts
    counts = {
        item['_id']: item['count']
        for item in Entry.objects(upload_id__in=upload_ids).aggregate(
            [{'$group': {'_id': '$upload_id', 'count': {'$sum': 1}}}]
        )
    }

    data = [upload_to_pydantic(upload, include_total_count=False) for upload in uploads]
    for pydantic_upload in data:
        pydantic_upload.entries = counts.get(pydantic_upload.upload_id, 0)

    pagination_response = PaginationResponse(
        total=mongodb_query.count(), **pagination.dict()
    )
    pagination_response.populate_simple_index_and_urls(request)

    return UploadProcDataQueryResponse(
        query=query, pagination=pagination_response, data=data
    )


@router.get(
    '/{upload_id}',
    tags=[APITag.METADATA],
    summary='Get a specific upload',
    response_model=UploadProcDataResponse,
    responses=create_responses(_upload_not_found, _not_authorized_to_upload),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def get_upload(
    upload_id: Annotated[
        str, Path(description='The unique id of the upload to retrieve.')
    ],
    user: Annotated[User, Depends(get_current_user([Scope.UPLOADS_READ]))],
):
    """
    Fetches a specific upload by its upload_id.
    """
    # Get upload (or throw exception if nonexistent/no access)
    upload = get_upload_with_read_access(upload_id, user, include_others=True)

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.get(
    '/{upload_id}/entries',
    tags=[APITag.METADATA],
    summary='Get the entries of the specific upload as a list',
    response_model=EntryProcDataQueryResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_pagination
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def get_upload_entries(
    request: Request,
    upload_id: Annotated[
        str, Path(description='The unique id of the upload to retrieve entries for.')
    ],
    pagination: Annotated[
        EntryProcDataPagination, Depends(entry_proc_data_pagination_parameters)
    ],
    user: Annotated[User, Depends(get_current_user([Scope.UPLOADS_READ]))],
):
    """
    Fetches the entries of a specific upload. Pagination is used to browse through the
    results.
    """
    upload = get_upload_with_read_access(upload_id, user, include_others=True)

    order_by = pagination.order_by
    assert order_by is not None
    order_by_with_sign = (
        order_by if pagination.order == Direction.asc else '-' + order_by
    )

    start = pagination.get_simple_index()
    end = start + pagination.page_size

    # load upload's entries. Use entry_id as tie breaker for ordering.
    entries = list(
        upload.entries_sublist(start, end, order_by=(order_by_with_sign, 'entry_id'))
    )
    failed_entries_count = upload.failed_entries_count

    # load entries's metadata from search
    metadata_entries_query = WithQuery(
        query={'entry_id:any': list(entry.entry_id for entry in entries)}
    ).query
    metadata_entries = search(
        pagination=MetadataPagination(page_size=len(entries)),
        owner='admin' if user is not None and user.is_admin else 'visible',
        user_id=user.user_id if user is not None else None,
        query=metadata_entries_query,
    )
    metadata_entries_map = {
        metadata_entry['entry_id']: metadata_entry
        for metadata_entry in metadata_entries.data
    }

    # convert data to pydantic
    data = []
    for entry in entries:
        pydantic_entry = entry_to_pydantic(entry)
        pydantic_entry.entry_metadata = metadata_entries_map.get(entry.entry_id)
        data.append(pydantic_entry)

    pagination_response = PaginationResponse(
        total=upload.total_entries_count, **pagination.dict()
    )
    pagination_response.populate_simple_index_and_urls(request)

    return EntryProcDataQueryResponse(
        pagination=pagination_response,
        processing_successful=upload.processed_entries_count - failed_entries_count,
        processing_failed=failed_entries_count,
        upload=upload_to_pydantic(upload),
        data=data,
    )


@router.get(
    '/{upload_id}/entries/{entry_id}',
    tags=[APITag.METADATA],
    summary='Get a specific entry for a specific upload',
    response_model=EntryProcDataResponse,
    responses=create_responses(_entry_not_found, _not_authorized_to_entry),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def get_upload_entry(
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    entry_id: Annotated[
        str,
        Path(
            description='The unique id of the entry, belonging to the specified upload.'
        ),
    ],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_READ], allow_anonymous=False)),
    ],
):
    """
    Fetches a specific entry for a specific upload.
    """
    upload = get_upload_with_read_access(upload_id, user)
    entry = upload.get_entry(entry_id)
    if not entry:
        raise HTTPException(
            status.HTTP_404_NOT_FOUND,
            detail=strip(
                """
            An entry by that id could not be found in the specified upload."""
            ),
        )

    data = entry_to_pydantic(entry, add_es_metadata=True, user=user)

    return EntryProcDataResponse(entry_id=entry_id, data=data)


@router.get(
    '/{upload_id}/rawdir/{path:path}',
    tags=[APITag.RAW],
    summary='Get the metadata for the raw file or folder located at the specified path in the specified upload.',
    response_model=RawDirResponse,
    responses=create_responses(
        _upload_or_path_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def get_upload_rawdir_path(
    request: Request,
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    path: Annotated[str, Path(description='The path within the upload raw files.')],
    pagination: Annotated[RawDirPagination, Depends(rawdir_pagination_parameters)],
    user: Annotated[User, Depends(get_current_user([Scope.UPLOADS_READ]))],
    include_entry_info: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """
                If the fields `entry_id` and `parser_name` should be populated for all
                encountered mainfiles."""
            )
        ),
    ] = False,
):
    """
    For the upload specified by `upload_id`, gets the raw file or directory metadata
    located at the given `path`. The response will either contain a `file_metadata` or
    `directory_metadata` key. For files, basic data about the file is returned, such as its
    name and size. For directories, the response includes a list of elements
    (files and folders) in the directory. For directories, the result is paginated.
    """
    # Get upload
    upload = get_upload_with_read_access(upload_id, user, include_others=True)
    upload_files = None
    try:
        # Get upload files
        upload_files = upload.upload_files
        if not upload_files.raw_exists(path):
            raise HTTPException(
                status.HTTP_404_NOT_FOUND,
                detail=strip(
                    """
                Not found. Invalid path?"""
                ),
            )

        response = RawDirResponse(
            path=path.rstrip('/'),
            access='unpublished'
            if not upload.published
            else ('embargoed' if upload.embargo_length else 'public'),
        )

        if upload_files.raw_isfile(path):
            # Path denotes a file
            response.file_metadata = RawDirFileMetadata(
                name=os.path.basename(path), size=upload_files.raw_file_size(path)
            )
            if include_entry_info:
                entry: Entry = Entry.objects(  # type: ignore
                    upload_id=upload_id, mainfile=path, mainfile_key=None
                ).first()
                if entry:
                    response.file_metadata.entry_id = entry.entry_id
                    response.file_metadata.parser_name = entry.parser_name
        else:
            # Path denotes a directory
            start = pagination.get_simple_index()
            end = start + pagination.page_size
            directory_list = upload_files.raw_listdir(path)
            upload_files.close()
            content = []
            path_to_element: dict[str, RawDirElementMetadata] = {}
            total = 0
            total_size = 0
            for i, path_info in enumerate(directory_list):
                total += 1
                total_size += path_info.size
                if start <= i < end:
                    element = RawDirElementMetadata(
                        name=os.path.basename(path_info.path),
                        is_file=path_info.is_file,
                        size=path_info.size,
                    )
                    content.append(element)
                    if include_entry_info:
                        path_to_element[path_info.path] = element

            if include_entry_info and content:
                for entry in Entry.objects(  # type: ignore  # type: ignore
                    upload_id=upload_id,
                    mainfile__in=path_to_element.keys(),
                    mainfile_key=None,
                ):
                    element = path_to_element[entry.mainfile]
                    element.entry_id = entry.entry_id
                    element.parser_name = entry.parser_name

            response.directory_metadata = RawDirDirectoryMetadata(
                name=os.path.basename(path), size=total_size, content=content
            )

            pagination_response = PaginationResponse(total=total, **pagination.dict())
            pagination_response.populate_simple_index_and_urls(request)
            response.pagination = pagination_response

        return response
    except Exception:
        if upload_files:
            upload_files.close()
        raise


@router.get(
    '/{upload_id}/raw',
    tags=[APITag.RAW],
    summary='Downloads the published upload .zip file with all the raw files of the upload.',
    response_class=StreamingResponse,
    responses=create_responses(
        _raw_path_response,
        _upload_or_path_not_found,
        _not_authorized_to_upload,
        _bad_request,
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
async def get_upload_raw(
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    user: Annotated[User, Depends(get_current_user([Scope.UPLOADS_READ]))],
):
    """
    NOMAD manages the raw files of published uploads as a .zip file. This endpoint
    allows to download it. While the outcome is similar to using `/uploads/<upload_id>/raw/`
    which creates a .zip file on the fly, this endpoint is more efficient
    because it simply streams an already existing .zip file. On the other hand, this
    endpoint is only available for published uploads and does not allow to selectively
    filter the files.
    """

    # Get upload
    upload = get_upload_with_read_access(upload_id, user, include_others=True)
    # Get upload files
    upload_files = upload.upload_files
    if not isinstance(upload_files, PublicUploadFiles):
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=strip(
                """
            Cannot download raw files .zip from non published uploads. Use '/{upload_id}/raw/' instead
            to recursively create and download a .zip file with all files."""
            ),
        )

    if FSUtility.is_local(file_path := upload_files.raw_zip_file_object().os_path):
        return FileResponse(file_path, media_type='application/zip')

    async def file_stream():
        with FSUtility.open(file_path) as file_obj:
            while chunk := file_obj.read(2**20):
                yield chunk
                await sleep(0)

    return StreamingResponse(file_stream(), media_type='application/zip')


@router.get(
    '/{upload_id}/raw/{path:path}',
    tags=[APITag.RAW],
    summary='Download the raw file or folder located at the specified path in the specified upload.',
    response_class=StreamingResponse,
    responses=create_responses(
        _raw_path_response,
        _upload_or_path_not_found,
        _not_authorized_to_upload,
        _bad_request,
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def get_upload_raw_path(
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    path: Annotated[str, Path(description='The path within the upload raw files.')],
    files_params: Annotated[Files, Depends(files_parameters)],
    user: Annotated[User, Depends(get_current_user([Scope.UPLOADS_READ]))],
    offset: Annotated[
        int | None,
        FastApiQuery(
            description=strip(
                """
                When dowloading individual files with `compress = false`, this can be
                used to seek to a specified position within the file in question. Default
                is 0, i.e. the start of the file."""
            )
        ),
    ] = 0,
    length: Annotated[
        int | None,
        FastApiQuery(
            description=strip(
                """
                When dowloading individual files with `compress = false`, this can be
                used to specify the number of bytes to read. By default, the value is -1,
                which means that the remainder of the file is streamed."""
            )
        ),
    ] = -1,
    decompress: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """
                Set if compressed files should be decompressed before streaming the
                content (that is: if there are compressed files *within* the raw files).
                Note, only some compression formats are supported."""
            )
        ),
    ] = False,
    ignore_mime_type: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """
                Sets the mime type specified in the response headers to `application/octet-stream`
                instead of the actual mime type."""
            )
        ),
    ] = False,
):
    """
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
    """
    if files_params.compress and (offset != 0 or length != -1):
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=strip(
                """
            Cannot specify `offset` or `length` when `compress` is true"""
            ),
        )
    # Get upload
    upload = get_upload_with_read_access(upload_id, user, include_others=True)
    # Get upload files
    upload_files = upload.upload_files
    try:
        if not upload_files.raw_exists(path):
            raise HTTPException(
                status.HTTP_404_NOT_FOUND,
                detail=strip(
                    """
                Not found. Invalid path?"""
                ),
            )
        if upload_files.raw_isfile(path):
            # File
            if files_params.compress:
                media_type = 'application/zip'
                content = create_download_stream_zipped(
                    DownloadItem(
                        upload_id=upload_id,
                        raw_path=path,
                        zip_path=os.path.basename(path),
                    ),
                    upload_files,
                    compress=True,
                )
            else:
                if offset is not None and offset < 0:
                    raise HTTPException(
                        status.HTTP_400_BAD_REQUEST,
                        detail=strip(
                            """
                        Invalid offset provided."""
                        ),
                    )
                if length is not None and length <= 0 and length != -1:
                    raise HTTPException(
                        status.HTTP_400_BAD_REQUEST,
                        detail=strip(
                            """
                        Invalid length provided. Should be greater than 0, or -1 if the remainder
                        of the file should be read."""
                        ),
                    )
                if ignore_mime_type or not (offset == 0 and length == -1):
                    media_type = 'application/octet-stream'
                else:
                    media_type = upload_files.raw_file_mime_type(path)
                content = create_download_stream_raw_file(
                    upload_files, path, offset, length, decompress
                )
            return StreamingResponse(
                content,
                headers=browser_download_headers(
                    filename=os.path.basename(path)
                    + ('.zip' if files_params.compress else ''),
                    media_type=media_type,
                ),
            )
        else:
            # Directory
            if not files_params.compress:
                raise HTTPException(
                    status.HTTP_400_BAD_REQUEST,
                    detail=strip(
                        """
                    Path is a directory, `compress` must be set to true"""
                    ),
                )
            # Stream directory content, compressed.
            return StreamingResponse(
                create_download_stream_zipped(
                    DownloadItem(upload_id=upload_id, raw_path=path, zip_path=''),
                    upload_files,
                    re_pattern=files_params.re_pattern,
                    recursive=True,
                    create_manifest_file=False,
                    compress=True,
                ),
                headers=browser_download_headers(
                    (
                        upload.upload_id
                        if not path
                        else os.path.basename(path.rstrip('/'))
                    )
                    + '.zip',
                    media_type='application/zip',
                ),
            )
    except Exception as e:
        if not isinstance(e, HTTPException):
            logger.error('exception while streaming download', exc_info=e)
        upload_files.close()
        raise


@router.put(
    '/{upload_id}/raw/{path:path}',
    tags=[APITag.RAW],
    summary='Upload a raw file to the specified path (directory) in the specified upload.',
    response_class=StreamingResponse,
    responses=create_responses(
        _put_raw_file_response,
        _upload_not_found,
        _not_authorized_to_upload,
        _bad_request,
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
async def put_upload_raw_path(
    request: Request,
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    path: Annotated[str, Path(description='The path within the upload raw files.')],
    user: Annotated[
        User,
        Depends(
            get_current_user(
                [Scope.UPLOADS_WRITE], allow_anonymous=False, allow_upload_token=True
            )
        ),
    ],
    file: Annotated[list[UploadFile] | None, File()] = None,
    local_path: Annotated[
        str | None, FastApiQuery(description=strip("""Internal/Admin use only."""))
    ] = None,
    file_name: Annotated[
        str | None,
        FastApiQuery(
            description=strip(
                """Specifies the name of the file, when using method 2."""
            )
        ),
    ] = None,
    overwrite_if_exists: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """If set to True (default), overwrites the file if it already exists."""
            )
        ),
    ] = True,
    copy_or_move: Annotated[
        str | None,
        FastApiQuery(
            description=strip(
                """If moving or copying a file within the same upload, specify which operation to do: move or copy"""
            )
        ),
    ] = None,
    copy_or_move_source_path: Annotated[
        str | None,
        FastApiQuery(
            description=strip(
                """If moving or copying a file within the same upload, specify the path to the source file."""
            )
        ),
    ] = None,
    wait_for_processing: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """Waits for the processing to complete and return information about the outcome in the response (**USE WITH CARE**)."""
            )
        ),
    ] = False,
    include_archive: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """If the archive data should be included in the response when using `wait_for_processing` (**USE WITH CARE**)."""
            )
        ),
    ] = False,
    entry_hash: Annotated[
        str | None,
        FastApiQuery(description=strip("""The hash code of the not modified entry.""")),
    ] = None,
    auto_decompress: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """
            Automatically decompress uploaded files upon receiving (ZIP or TAR). True by default."""
            )
        ),
    ] = True,
    trigger_processing: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """
            If set to true (default), reprocesses the upload after deleting the file/folder."""
            ),
        ),
    ] = True,
):
    """
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
    """
    if include_archive and not wait_for_processing:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='`include_archive` requires `wait_for_processing`.',
        )
    if wait_for_processing and not trigger_processing:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='`trigger_processing` must be true when `wait_for_processing` is set.`',
        )

    upload = _get_upload_with_write_access(upload_id, user, include_published=False)

    if local_path and not os.path.isfile(local_path):
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='Uploading folders with local_path is not yet supported.',
        )

    if not is_safe_relative_path(path):
        raise HTTPException(status.HTTP_400_BAD_REQUEST, detail='Bad path provided.')

    if copy_or_move is not None or copy_or_move_source_path is not None:
        if copy_or_move not in ['copy', 'move']:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail="The copy_or_move query parameter should be one of 'copy' or 'move' options.",
            )

        if (
            copy_or_move is None
            or copy_or_move_source_path is None
            or file_name is None
        ):
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail="""For a successful copy/move operation, all three query parameters: file_name, copy_or_move and copy_or_move_source_path are required.""",
            )

        if not is_safe_basename(file_name):
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST, detail='Bad file name provided'
            )

        if not is_safe_relative_path(copy_or_move_source_path):
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='Bad source path provided.',
            )

    upload_paths, _, method = await _get_files_if_provided(
        upload_id, request, file, local_path, file_name, user
    )

    if not upload_paths and not (
        copy_or_move and copy_or_move_source_path and file_name
    ):
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='Either an upload file or the query parameters for moving/copying a file should be provided.',
        )

    if entry_hash:
        upload_path = upload_paths[0]
        full_path = os.path.join(path, os.path.basename(upload_path))
        entry_id = utils.generate_entry_id(upload_id, full_path)
        entry = upload.get_entry(entry_id)
        if entry and entry_hash != entry.entry_hash or not entry:
            raise HTTPException(
                status.HTTP_409_CONFLICT,
                detail='The provided hash did not match the current file.',
            )

    upload_files = StagingUploadFiles(upload_id)

    compression_format = None
    for upload_path in upload_paths:
        compression_format = get_compression_format(upload_path)
        if compression_format == 'error':
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='Cannot extract file. Bad file format or file extension?',
            )
        if not compression_format and not overwrite_if_exists:
            full_path = os.path.join(path, os.path.basename(upload_path))
            if upload_files.raw_exists(full_path):
                raise HTTPException(
                    status.HTTP_409_CONFLICT,
                    detail='The provided path already exists and overwrite_if_exists is set to False.',
                )

    if not wait_for_processing:
        # Process on worker (normal case)
        if copy_or_move:  # the case for move/copy an existing file
            path_to_target_file = os.path.join(path, file_name)
            if upload_files.raw_exists(path_to_target_file):
                raise HTTPException(
                    status.HTTP_409_CONFLICT,
                    detail='The provided path already exists.',
                )
            if not upload_files.raw_exists(path):
                raise HTTPException(
                    status.HTTP_404_NOT_FOUND,
                    detail='No file or folder with that path found.',
                )
            if not upload_files.raw_exists(copy_or_move_source_path):
                raise HTTPException(
                    status.HTTP_409_CONFLICT,
                    detail=f'No file or folder with that source path: {copy_or_move_source_path}',
                )
            file_operations: Any = [
                dict(
                    op=copy_or_move.upper(),
                    path_to_existing_file=copy_or_move_source_path,
                    path_to_target_file=path_to_target_file,
                )
            ]
        else:
            file_operations = [
                dict(
                    op='ADD',
                    path=upload_path,
                    target_dir=path,
                    temporary=(method != 0),
                    auto_decompress=auto_decompress,
                )
                for upload_path in upload_paths
            ]

        # Initiate processing
        try:
            upload.process_upload(
                file_operations=file_operations,
                only_updated_files=True,
                trigger_processing=trigger_processing,
            )
        except ProcessAlreadyRunning:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='The upload is currently blocked by another process.',
            )

        # Create response
        if request.headers.get('Accept') == 'application/json':
            response = PutRawFileResponse(
                upload_id=upload_id, data=upload_to_pydantic(upload)
            )
            response_text = response.model_dump_json()
            media_type = 'application/json'
        else:
            response_text = _thank_you_message
            media_type = 'text/plain'

        return StreamingResponse(
            create_stream_from_string(response_text), media_type=media_type
        )

    # Process locally
    if copy_or_move:  # case for move/copy an existing file
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='Cannot move/copy the file with wait_for_processing set to true.',
        )

    if len(upload_paths) != 1 or compression_format:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='`wait_for_processing` can only be used with single files, and not with compressed files.',
        )

    upload_path = upload_paths[0]
    full_path = os.path.join(path, os.path.basename(upload_path))
    try:
        entry = upload.put_file_and_process_local(
            upload_path,
            path,
            reprocess_settings=Reprocess(
                index_individual_entries=True, reprocess_existing_entries=True
            ),
        )

        search_refresh()

        archive = None
        if entry and entry.process_status == ProcessStatus.SUCCESS and include_archive:
            # NOTE: We can't rely on ES to get the metadata for the entry, since it may
            # not have had enough time to update its index etc. For now, we will just
            # ignore this, as we do not need it.
            archive = answer_entry_archive_request(
                dict(upload_id=upload_id, mainfile=full_path),
                required='*',
                user=user,
                entry_metadata=dict(
                    upload_id=upload_id,
                    entry_id=entry.entry_id,
                    parser_name=entry.parser_name,
                ),
            )['data']['archive']

        response = PutRawFileResponse(
            upload_id=upload_id,
            data=upload_to_pydantic(upload),
            processing=ProcessingData(
                upload_id=upload_id,
                path=full_path,
                entry_id=entry.entry_id if entry else None,
                parser_name=entry.parser_name if entry else None,
                entry=entry_to_pydantic(entry) if entry else None,
                archive=archive,
            ),
        )

        return StreamingResponse(
            create_stream_from_string(response.model_dump_json()),
            media_type='application/json',
        )
    except HTTPException:
        raise
    except ProcessAlreadyRunning:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='The upload is currently being processed, operation not allowed.',
        )
    except Exception as e:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=f'Unexpected exception occurred: {e}',
        )
    finally:
        try:
            shutil.rmtree(os.path.dirname(upload_path))
        except Exception:  # noqa
            pass


@router.delete(
    '/{upload_id}/raw/{path:path}',
    tags=[APITag.RAW],
    summary='Delete the raw file or folder located at the specified path in the specified upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def delete_upload_raw_path(
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    path: Annotated[str, Path(description='The path within the upload raw files.')],
    user: Annotated[
        User,
        Depends(
            get_current_user(
                [Scope.UPLOADS_WRITE], allow_anonymous=False, allow_upload_token=True
            )
        ),
    ],
    trigger_processing: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """
            If set to true (default), reprocesses the upload after deleting the file/folder."""
            ),
        ),
    ] = True,
):
    """
    Delete file or folder located at the specified path in the specified upload. The upload
    must not be published. This also automatically triggers a reprocessing of the upload.
    Choosing the empty string as `path` deletes all files.
    """
    upload = _get_upload_with_write_access(upload_id, user, include_published=False)

    if not is_safe_relative_path(path):
        raise HTTPException(status.HTTP_400_BAD_REQUEST, detail='Bad path provided.')

    upload_files = StagingUploadFiles(upload_id)

    if not upload_files.raw_exists(path):
        raise HTTPException(
            status.HTTP_404_NOT_FOUND,
            detail='No file or folder with that path found.',
        )

    try:
        upload.process_upload(
            file_operations=[dict(op='DELETE', path=path)],
            only_updated_files=True,
            trigger_processing=trigger_processing,
        )
    except ProcessAlreadyRunning:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='The upload is currently blocked by another process.',
        )

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/raw-create-dir/{path:path}',
    tags=[APITag.RAW],
    summary='Create a new empty directory with the specified path in the specified upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def post_upload_raw_create_dir_path(
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    path: Annotated[str, Path(description='The path within the upload raw files.')],
    user: Annotated[
        User,
        Depends(
            get_current_user(
                [Scope.UPLOADS_WRITE], allow_anonymous=False, allow_upload_token=True
            )
        ),
    ],
):
    """
    Create a new empty directory in the specified upload. The `path` should be the full path
    to the new directory (i.e. ending with the name of the new directory). The api call returns
    immediately (no processing is necessary).
    """
    upload = _get_upload_with_write_access(upload_id, user, include_published=False)

    if not path or not is_safe_relative_path(path):
        raise HTTPException(status.HTTP_400_BAD_REQUEST, detail='Bad path provided.')
    if upload.staging_upload_files.raw_exists(path):
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=f'Path `{path}` already exists.',
        )
    try:
        upload.staging_upload_files.raw_create_directory(path)
        # No real processing is needed when just adding a folder, but we should signal that
        # the upload has changed.
        upload.complete_time = datetime.now(timezone.utc)
        upload.save()
    except Exception as e:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=f'Failed to create directory: {e}',
        )

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.get(
    '/{upload_id}/archive/mainfile/{mainfile:path}',
    tags=[APITag.ARCHIVE],
    summary='Get the full archive for the given upload and mainfile path.',
    response_model=EntryArchiveResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_upload_or_path_not_found, _not_authorized_to_upload),
)
def get_upload_entry_archive_mainfile(
    user: Annotated[User, Depends(get_current_user([Scope.UPLOADS_READ]))],
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    mainfile: Annotated[
        str, Path(description="The mainfile path within the upload's raw files.")
    ],
    mainfile_key: Annotated[
        str | None,
        FastApiQuery(description='The mainfile_key, for accessing child entries.'),
    ] = None,
):
    """
    For the upload specified by `upload_id`, gets the full archive of a single entry that
    is identified by the given `mainfile`.
    """
    get_upload_with_read_access(upload_id, user, include_others=True)
    query = dict(upload_id=upload_id, mainfile=mainfile)
    if mainfile_key:
        query.update(mainfile_key=mainfile_key)
    return answer_entry_archive_request(query, required='*', user=user)


@router.get(
    '/{upload_id}/archive/{entry_id}',
    tags=[APITag.ARCHIVE],
    summary='Get the full archive for the given upload and entry.',
    response_model=EntryArchiveResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
    responses=create_responses(_upload_or_path_not_found, _not_authorized_to_upload),
)
def get_upload_entry_archive(
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    entry_id: Annotated[str, Path(description='The unique entry id.')],
    user: Annotated[User, Depends(get_current_user([Scope.UPLOADS_READ]))],
):
    """
    For the upload specified by `upload_id`, gets the full archive of a single entry that
    is identified by the given `entry_id`.
    """
    get_upload_with_read_access(upload_id, user, include_others=True)
    return answer_entry_archive_request(
        dict(upload_id=upload_id, entry_id=entry_id), required='*', user=user
    )


@router.post(
    '',
    tags=[APITag.DEFAULT],
    summary='Submit a new upload',
    response_class=StreamingResponse,
    responses=create_responses(_post_upload_response, _not_authorized, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
async def post_upload(
    request: Request,
    user: Annotated[
        User,
        Depends(
            get_current_user(
                [Scope.UPLOADS_WRITE], allow_anonymous=False, allow_upload_token=True
            )
        ),
    ],
    file: Annotated[list[UploadFile] | None, File()] = None,
    local_path: Annotated[
        str | None,
        FastApiQuery(
            description=strip(
                """
            Internal/Admin use only."""
            )
        ),
    ] = None,
    example_upload_id: Annotated[
        str | None,
        FastApiQuery(
            description=strip(
                """
            If provided, instantiates a new upload from the given example upload
            entry point id. You may use this parameter in combination with other
            file sources.
            """
            )
        ),
    ] = None,
    file_name: Annotated[
        str | None,
        FastApiQuery(
            description=strip(
                """
            Specifies the name of the file, when using method 2."""
            )
        ),
    ] = None,
    upload_name: Annotated[
        str | None,
        FastApiQuery(
            description=strip(
                """
            A human readable name for the upload."""
            )
        ),
    ] = None,
    embargo_length: Annotated[
        int,
        FastApiQuery(
            description=strip(
                """
            The requested embargo length, in months, if any (0-36)."""
            )
        ),
    ] = 0,
    publish_directly: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
            If the upload should be published directly. False by default."""
            )
        ),
    ] = None,
    auto_decompress: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """
            Automatically decompress uploaded files upon receiving (ZIP or TAR). True by default."""
            )
        ),
    ] = True,
):
    """
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
    """
    if not user.is_admin:
        # Check upload limit
        if (
            _query_mongodb(main_author=str(user.user_id), publish_time=None).count()
            >= config.services.upload_limit
        ):  # type: ignore
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail=strip(
                    """
                Limit of unpublished uploads exceeded for user."""
                ),
            )

    if not 0 <= embargo_length <= 36:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='`embargo_length` must be between 0 and 36 months.',
        )

    upload_id = utils.create_uuid()

    upload_paths, upload_folders, method = await _get_files_if_provided(
        upload_id, request, file, local_path, file_name, user
    )

    if not upload_name:
        # Try to default upload_name
        if example_upload_id:
            try:
                entry_point = cast(
                    ExampleUploadEntryPoint,
                    config.get_plugin_entry_point(example_upload_id),
                )
            except Exception:
                raise HTTPException(
                    status.HTTP_400_BAD_REQUEST,
                    detail=f'Could not find example upload with id "{example_upload_id}"',
                )
            upload_name = entry_point.title
        elif method == 2:
            upload_name = file_name or None
        elif len(upload_paths) == 1:
            upload_name = os.path.basename(upload_paths[0])

    upload: Upload = Upload.create(
        upload_id=upload_id,
        main_author=user,
        upload_name=upload_name,
        upload_create_time=datetime.now(timezone.utc),
        embargo_length=embargo_length,
        publish_directly=publish_directly,
    )

    # Create staging files
    files.StagingUploadFiles(upload_id=upload_id, create=True)

    logger.info('upload created', upload_id=upload_id)

    file_operations = [
        dict(
            op='ADD',
            path=upload_path,
            target_dir=upload_folders[i_path],
            temporary=(method != 0),
            auto_decompress=auto_decompress,
        )
        for i_path, upload_path in enumerate(upload_paths)
    ]

    # If creating an example upload, the contents are loaded only during the
    # first processing: they should not be loaded anymore in later reprocessing.
    if example_upload_id is not None:
        upload.process_example_upload(example_upload_id, file_operations)
    elif upload_paths:
        upload.process_upload(file_operations)

    if request.headers.get('Accept') == 'application/json':
        upload_proc_data_response = UploadProcDataResponse(
            upload_id=upload_id, data=upload_to_pydantic(upload)
        )
        response_text = upload_proc_data_response.model_dump_json()
        media_type = 'application/json'
    else:
        response_text = _thank_you_message
        media_type = 'text/plain'

    return StreamingResponse(
        create_stream_from_string(response_text), media_type=media_type
    )


@router.post(
    '/{upload_id}/edit',
    tags=[APITag.METADATA],
    summary='Updates the metadata of the specified upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
async def post_upload_edit(
    request: Request,
    data: MetadataEditRequest,
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_WRITE], allow_anonymous=False)),
    ],
    wait_for_processing: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """Waits for the processing to complete and return information about the outcome in the response (**USE WITH CARE**)."""
            )
        ),
    ] = True,
):
    """
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
    """
    edit_request_json = await request.json()
    try:
        MetadataEditRequestHandler.edit_metadata(
            edit_request_json, upload_id, user, wait_for_processing=wait_for_processing
        )
        return UploadProcDataResponse(
            upload_id=upload_id, data=upload_to_pydantic(Upload.get(upload_id))
        )
    except RequestValidationError:
        raise  # A problem which we have handled explicitly. Fastapi does json conversion.
    except Exception as e:
        # The upload is processing or some kind of unexpected error has occurred
        raise HTTPException(status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.delete(
    '/{upload_id}',
    tags=[APITag.DEFAULT],
    summary='Delete an upload',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def delete_upload(
    upload_id: Annotated[
        str, Path(description='The unique id of the upload to delete.')
    ],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_WRITE], allow_anonymous=False)),
    ],
    wait_for_processing: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """Waits for the processing to complete and return information about the outcome in the response (**USE WITH CARE**)."""
            )
        ),
    ] = True,
):
    """
    Delete an existing upload.

    Only uploads that are sill in staging, not already deleted, not still uploaded, and
    not currently processed, can be deleted.
    """
    upload = _get_upload_with_write_access(
        upload_id,
        user,
        include_published=True,
        published_requires_admin=True,
        include_failed_imports=True,
        only_main_author=True,
    )
    try:
        upload.delete_upload(wait_for_processing=wait_for_processing)
    except ProcessAlreadyRunning:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=strip(
                """
            The upload is still being processed."""
            ),
        )
    except Exception as e:
        logger.error('could not delete processing upload', exc_info=e)
        raise

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/publish',
    tags=[APITag.ACTION],
    summary='Publish an upload',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def post_upload_action_publish(
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_PUBLISH], allow_anonymous=False)),
    ],
    upload_id: Annotated[
        str,
        Path(
            description=strip(
                """
                The unique id of the upload to publish."""
            )
        ),
    ],
    embargo_length: Annotated[
        int | None,
        FastApiQuery(
            description=strip(
                """
                If provided, updates the embargo length of the upload. The value should
                be between 0 and 36 months. 0 means no embargo."""
            )
        ),
    ] = None,
    to_central_nomad: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """
            DEPRECATED
            To publish to an external oasis or to the central nomad you can use the new entpoint /uploads/{upload_id}/action/transfer.
                """
            ),
            deprecated=True,
        ),
    ] = False,
    wait_for_processing: Annotated[
        bool,
        FastApiQuery(
            description=strip(
                """Waits for the processing to complete and return information about the outcome in the response (**USE WITH CARE**)."""
            )
        ),
    ] = True,
):
    """
    Publishes an upload. The upload cannot be modified after this point (except for special
    cases, like when lifting the embargo prematurely, and by admins). After the upload is
    published and the embargo period (if any) is expired, the generated archive entries
    will be publicly visible.
    """
    upload = _get_upload_with_write_access(
        upload_id, user, include_published=True, published_requires_admin=False
    )

    if upload.published and not user.is_admin and not to_central_nomad:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST, detail='Upload already published.'
        )

    _check_upload_not_processing(upload)

    if upload.process_status == ProcessStatus.FAILURE:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='Cannot publish an upload that failed processing.',
        )
    if upload.processed_entries_count == 0:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='Cannot publish an upload without any resulting entries.',
        )
    if embargo_length is not None and not 0 <= embargo_length <= 36:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='Invalid embargo_length. Must be between 0 and 36 months.',
        )

    if to_central_nomad:
        # Publish from an OASIS to the central repository
        if not config.oasis.is_oasis:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='Must be on an OASIS to publish to the central NOMAD repository.',
            )
        if not upload.published:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='The upload must be published on the OASIS first.',
            )
        if not user.is_admin:
            raise HTTPException(
                status.HTTP_403_FORBIDDEN,
                detail='Only admin of OASIS can publish to the central NOMAD.',
            )
        # Everything looks ok, try to publish it to the central NOMAD!
        upload.publish_externally(embargo_length=embargo_length)
    else:
        # Publish to this repository
        if upload.published:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='The upload is already published.',
            )
        try:
            upload.publish_upload(
                embargo_length=embargo_length, wait_for_processing=wait_for_processing
            )
        except ProcessAlreadyRunning:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='The upload is still/already processed.',
            )

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/process',
    tags=[APITag.ACTION],
    summary='Manually triggers processing of an upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def post_upload_action_process(
    upload_id: Annotated[
        str, Path(description='The unique id of the upload to process.')
    ],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_PROCESS], allow_anonymous=False)),
    ],
):
    """
    Processes an upload, i.e. parses the files and updates the NOMAD archive. Only admins
    can process an already published upload.
    """
    upload = _get_upload_with_write_access(
        upload_id, user, include_published=True, published_requires_admin=True
    )

    _check_upload_not_processing(upload)

    upload.process_upload()
    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/delete-entry-files',
    tags=[APITag.ACTION],
    summary='Deletes the files of the entries specified by a query.',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def post_upload_action_delete_entry_files(
    data: DeleteEntryFilesRequest,
    upload_id: Annotated[
        str,
        Path(
            description='The unique id of the upload within which to delete entry files.'
        ),
    ],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_WRITE], allow_anonymous=False)),
    ],
):
    """Deletes the files of the entries specified by the provided query."""

    upload = _get_upload_with_write_access(upload_id, user, include_published=False)

    # Evaluate query
    if not data.query:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=strip(
                """
            A query must be specified."""
            ),
        )
    restricted_query = restrict_query_to_upload(data.query, upload_id)
    es_entries = search_iterator(
        user_id=user.user_id,
        owner=data.owner,
        query=restricted_query,
        required=MetadataRequired(include=['mainfile']),
    )

    # Determine paths to delete
    try:
        paths_to_delete: set[str] = set()
        for es_entry in es_entries:
            mainfile = es_entry['mainfile']
            path_to_delete = (
                os.path.dirname(mainfile) if data.include_parent_folders else mainfile
            )
            paths_to_delete.add(path_to_delete)
    except QueryValidationError as e:
        raise RequestValidationError(errors=e.errors)

    # Execute operation
    if paths_to_delete:
        try:
            upload.process_upload(
                file_operations=[
                    dict(op='DELETE', path=path_to_delete)
                    for path_to_delete in sorted(paths_to_delete)
                ],
                only_updated_files=True,
            )
        except ProcessAlreadyRunning:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='The upload is currently blocked by another process.',
            )

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/lift-embargo',
    tags=[APITag.ACTION],
    summary='Lifts the embargo of an upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def post_upload_action_lift_embargo(
    upload_id: Annotated[
        str, Path(description='The unique id of the upload to lift the embargo for.')
    ],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_PUBLISH], allow_anonymous=False)),
    ],
):
    """Lifts the embargo of an upload."""
    upload = _get_upload_with_write_access(
        upload_id, user, include_published=True, published_requires_admin=False
    )
    _check_upload_not_processing(upload)
    if not upload.published:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=strip(
                """
            Upload is not published, no embargo to lift."""
            ),
        )
    if not upload.with_embargo:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=strip(
                """
            Upload has no embargo."""
            ),
        )
    # Lift the embargo using MetadataEditRequestHandler.edit_metadata
    try:
        MetadataEditRequestHandler.edit_metadata(
            {'metadata': {'embargo_length': 0}}, upload_id, user
        )
        upload.reload()
        return UploadProcDataResponse(
            upload_id=upload_id, data=upload_to_pydantic(upload)
        )
    except Exception as e:
        # Should only happen if the upload just started processing or something unexpected happens
        raise HTTPException(status.HTTP_400_BAD_REQUEST, detail=str(e))


@router.get(
    '/{upload_id}/bundle',
    tags=[APITag.BUNDLE],
    summary='Gets an *upload bundle* for the specified upload.',
    response_class=StreamingResponse,
    responses=create_responses(
        _upload_bundle_response,
        _upload_not_found,
        _not_authorized_to_upload,
        _bad_request,
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def get_upload_bundle(
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_BUNDLE_READ])),
    ],
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    include_raw_files: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If raw files should be included in the bundle (true by default)."""
            )
        ),
    ] = True,
    include_archive_files: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If archive files (i.e. parsed entries data) should be included in the bundle
                (true by default)."""
            )
        ),
    ] = True,
    include_datasets: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If datasets references to this upload should be included in the bundle
                (true by default)."""
            )
        ),
    ] = True,
):
    """
    Get an *upload bundle* for the specified upload. An upload bundle is a file bundle which
    can be used to export and import uploads between different NOMAD deployments.
    """
    upload = get_upload_with_read_access(upload_id, user, include_others=True)
    _check_upload_not_processing(upload)

    export_settings = config.bundle_export.default_settings.customize(
        dict(
            include_raw_files=include_raw_files,
            include_archive_files=include_archive_files,
            include_datasets=include_datasets,
        )
    )

    try:
        stream = BundleExporter(
            upload,
            export_as_stream=True,
            export_path=None,
            zipped=True,
            overwrite=False,
            export_settings=export_settings,
        ).export_bundle()
    except Exception as e:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=strip(f'Could not export due to error: {e}'),
        )

    return StreamingResponse(async_wrapper(stream), media_type='application/zip')


@router.post(
    '/bundle',
    tags=[APITag.BUNDLE],
    summary='Posts an *upload bundle* to this NOMAD deployment.',
    response_model=UploadProcDataResponse,
    responses=create_responses(_not_authorized, _bad_request),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
async def post_upload_bundle(
    request: Request,
    user: Annotated[
        User,
        Depends(
            get_current_user(
                [Scope.UPLOADS_BUNDLE_WRITE],
                allow_anonymous=False,
                allow_upload_token=True,
            )
        ),
    ],
    file: Annotated[list[UploadFile] | None, File()] = None,
    local_path: Annotated[
        str | None,
        FastApiQuery(
            description=strip(
                """
            Internal/Admin use only."""
            )
        ),
    ] = None,
    embargo_length: Annotated[
        int | None,
        FastApiQuery(
            description=strip(
                """
                Specifies the embargo length in months to set on the upload. If omitted,
                the value specified in the bundle will be used. A value of 0 means no
                embargo."""
            )
        ),
    ] = None,
    include_raw_files: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If raw files should be imported from the bundle
                *(only admins can change this setting)*."""
            )
        ),
    ] = None,
    include_archive_files: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If archive files (i.e. parsed entries data) should be imported from the bundle
                *(only admins can change this setting)*."""
            )
        ),
    ] = None,
    include_datasets: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If dataset references to this upload should be imported from the bundle
                *(only admins can change this setting)*."""
            )
        ),
    ] = None,
    include_bundle_info: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If the bundle_info.json file should be kept
                *(only admins can change this setting)*."""
            )
        ),
    ] = None,
    keep_original_timestamps: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If all original timestamps, including `upload_create_time`, `entry_create_time`
                and `publish_time`, should be kept
                *(only admins can change this setting)*."""
            )
        ),
    ] = None,
    set_from_oasis: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If the `from_oasis` flag and `oasis_deployment_url` should be set
                *(only admins can change this setting)*."""
            )
        ),
    ] = None,
    trigger_processing: Annotated[
        bool | None,
        FastApiQuery(
            description=strip(
                """
                If processing should be triggered after the bundle has been imported
                *(only admins can change this setting)*."""
            )
        ),
    ] = None,
):
    """
    Posts an *upload bundle* to this NOMAD deployment. An upload bundle is a file bundle which
    can be used to export and import uploads between different NOMAD installations. The
    endpoint expects an upload bundle attached as a zipfile.

    **NOTE:** This endpoint is restricted to admin users and oasis admins. Further, all
    settings except `embargo_length` requires an admin user to change (these settings
    have default values specified by the system configuration).

    There are two basic ways to upload files: in the multipart-formdata or streaming the
    file data in the http body. Both are supported. See the POST `uploads` endpoint for
    examples of curl commands for uploading files.
    """
    import_settings = config.bundle_import.default_settings.customize(
        dict(
            include_raw_files=include_raw_files,
            include_archive_files=include_archive_files,
            include_datasets=include_datasets,
            include_bundle_info=include_bundle_info,
            keep_original_timestamps=keep_original_timestamps,
            set_from_oasis=set_from_oasis,
            trigger_processing=trigger_processing,
        )
    )

    bundle_importer: BundleImporter | None = None
    bundle_path: str | None = None
    method = None

    if local_path:
        if not os.path.isfile(local_path):
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='You can only target a single bundle file using local_path.',
            )

    try:
        bundle_importer = BundleImporter(user, import_settings)
        bundle_importer.check_api_permissions()

        bundle_paths, _, method = await _get_files_if_provided(
            tmp_dir_prefix='bundle',
            request=request,
            file=file,
            local_path=local_path,
            file_name=None,
            user=user,
        )

        if not bundle_paths:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='No bundle file provided',
            )
        if len(bundle_paths) > 1:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='Can only provide one bundle file at a time',
            )
        bundle_path = bundle_paths[0]

        bundle_importer.open(bundle_path)
        upload = bundle_importer.create_upload_skeleton()
        bundle_importer.close()
        # Import the bundle using the unified method
        upload.import_bundle(
            bundle_path=bundle_path,
            import_settings=import_settings.model_dump()
            if import_settings is not None
            else {},
            embargo_length=embargo_length,
        )

        return UploadProcDataResponse(
            upload_id=upload.upload_id, data=upload_to_pydantic(upload)
        )
    except Exception as e:
        if bundle_importer:
            bundle_importer.close()
            if bundle_path and method != 0:
                bundle_importer.delete_bundle()
        if isinstance(e, HTTPException):
            raise
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=f'Failed to import bundle: {str(e)}',
        )


@router.post(
    '/{upload_id}/action/transfer',
    tags=[APITag.ACTION],
    summary='Transfer upload to another NOMAD deployment.',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def transfer_upload_bundle(
    transfer_options: TransferBundleRequest,
    upload_id: Annotated[
        str,
        Path(
            description=strip(
                """
                The unique id of the upload to transfer."""
            )
        ),
    ],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_BUNDLE_READ], allow_anonymous=False)),
    ],
):
    """
    Start a transfer of an upload to another NOMAD deployment.
    By default the transfer will target the central nomad deployment if no `target_deployment_url` is provided.
    `auth_token` is required to authenticate the transfer process in the target deploment (custom OASIS or central NOMAD).
    """
    upload = get_upload_with_read_access(
        upload_id=upload_id, user=user, include_others=True
    )
    if not upload.published:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='The upload should be published first.',
        )

    target_deployment_url = transfer_options.target_deployment_url

    _validate_target_deployment_url(target_deployment_url)
    _check_upload_not_processing(upload)
    _check_external_deployment_status(target_deployment_url)

    upload.publish_externally(
        target_deployment_url=target_deployment_url,
        auth_token=transfer_options.auth_token,
        embargo_length=transfer_options.embargo_length,
    )
    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


def _validate_target_deployment_url(url: str):
    # urlparse never raises, but may return empty parts
    parsed = urlparse(url)
    # Check scheme
    if parsed.scheme not in ('http', 'https'):
        raise HTTPException(
            status_code=422, detail='URL must start with http:// or https://'
        )

    # Check host
    if not parsed.netloc:
        raise HTTPException(
            status_code=422, detail='URL must contain a valid host (e.g., example.com)'
        )

    # Check that path ends with /api
    if not parsed.path.endswith('/api'):
        raise HTTPException(status_code=422, detail="URL path must end with '/api'")


async def _get_files_if_provided(
    tmp_dir_prefix: str,
    request: Request,
    file: list[UploadFile],
    local_path: str,
    file_name: str,
    user: User,
) -> tuple[list[str], list[str], None | int]:
    """
    If the user provides one or more files with the api call, load and save them to a temporary
    folder (or, if method 0 is used, just "forward" the file path). The method thus needs to identify
    which file transfer method was used (0 - 2), and save the data to disk (if method is 1 or 2).

    Returns a list of os paths to the resulting files and method (0-2), or ([], None) if no file
    data was provided with the api call.
    """
    # Determine the source data stream
    sources: list[tuple[Any, str]] = []  # List of tuples (source, filename)
    if local_path:
        # Method 0: Local file - only for admin use.
        if not user.is_admin:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail=strip("""
                You are not authorized to access this path.
                """),
            )
        if not os.path.exists(local_path):
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='The specified local_path cannot be found.',
            )
        method = 0
    elif file:
        # Method 1: Data provided as formdata
        method = 1

        async def _async_reader(_f):
            try:
                while _data := await _f.read(io.DEFAULT_BUFFER_SIZE):
                    yield _data
            except Exception as _e:
                raise _e
            finally:
                await _f.close()

        sources = [
            (_async_reader(multipart_file), unquote(multipart_file.filename))
            for multipart_file in file
        ]
    else:
        # Method 2: Data has to be sent streamed in the body
        method = 2
        sources = [(request.stream(), file_name or 'NO NAME')]

    no_file_name_info_provided = not file_name

    for _, source_file_name in sources:
        if not is_safe_basename(source_file_name):
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST, detail='Bad file name provided.'
            )

    # Forward the file path (if method == 0) or save the file(s)
    if method == 0:
        is_file = os.path.isfile(local_path)
        # Single file
        if is_file:
            upload_paths = [local_path]
            upload_folders = ['']
        # Folder
        else:
            upload_paths = []
            upload_folders = []
            for root, _, filepaths in os.walk(local_path):
                for uploaded_file in filepaths:
                    file_path = os.path.abspath(os.path.join(root, uploaded_file))
                    folder = os.path.relpath(root, local_path)
                    if folder == '.':
                        folder = ''
                    upload_paths.append(file_path)
                    upload_folders.append(folder)
    else:
        tmp_dir = files.mkdtemp(tmp_dir_prefix)
        upload_paths = []
        uploaded_bytes = 0
        upload_folders = []
        for source_stream, source_file_name in sources:
            upload_path = os.path.join(tmp_dir, source_file_name)
            try:
                with open(upload_path, 'wb') as f:
                    uploaded_bytes = 0
                    log_interval = 1e9
                    next_log_at = log_interval
                    async for chunk in source_stream:
                        if not chunk:
                            # End of data stream
                            break
                        uploaded_bytes += len(chunk)
                        f.write(chunk)
                        if uploaded_bytes > next_log_at:
                            logger.info(
                                'large upload in progress',
                                uploaded_bytes=uploaded_bytes,
                            )
                            next_log_at += log_interval
                    logger.info(f'upload completed', uploaded_bytes={uploaded_bytes})
            except Exception as e:
                if not (isinstance(e, RuntimeError) and 'Stream consumed' in str(e)):
                    if os.path.exists(tmp_dir):
                        shutil.rmtree(tmp_dir)
                    logger.warn('IO error receiving upload data', exc_info=e)
                    raise HTTPException(
                        status.HTTP_400_BAD_REQUEST,
                        detail='Some IO went wrong, upload probably aborted/disrupted.',
                    )
            upload_paths.append(upload_path)
            upload_folders.append('')

        if not uploaded_bytes and method == 2:
            # No data was provided
            shutil.rmtree(tmp_dir)
            return [], [], None

    logger.info(f'received uploaded file(s)')
    if method == 2 and no_file_name_info_provided:
        # Only ok if uploaded file is a zip or a tar archive.
        ext = (
            '.zip'
            if zipfile.is_zipfile(upload_path)
            else '.tar'
            if tarfile.is_tarfile(upload_path)
            else None
        )
        if not ext:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='No file name provided, and the file does not look like a zip or tar file.',
            )
        # Add the correct extension
        shutil.move(upload_path, upload_path + ext)
        upload_paths = [upload_path + ext]
        upload_folders = ['']

    return upload_paths, upload_folders, method


def _query_mongodb(**kwargs):
    return Upload.objects(**kwargs)  # type: ignore


def get_role_query(roles: list[UploadRole], user: User, include_all=False) -> Q:
    """
    Create MongoDB filter query for user with given roles (default: all roles)
    """
    if not roles:
        roles = list(UploadRole)

    group_ids = MongoUserGroup.get_ids_by_user_id(user.user_id, include_all=include_all)

    role_query = Q()
    if UploadRole.main_author in roles:
        role_query |= Q(main_author=user.user_id)
    if UploadRole.coauthor in roles:
        role_query |= Q(coauthors=user.user_id) | Q(coauthor_groups__in=group_ids)
    if UploadRole.reviewer in roles:
        role_query |= Q(reviewers=user.user_id) | Q(reviewer_groups__in=group_ids)

    return role_query


def is_user_upload_viewer(upload: Upload, user: User | None) -> bool:
    """Check whether user has access to that upload."""
    if 'all' in upload.reviewer_groups:
        return True

    if user is None:
        return False

    if user.is_admin:
        return True

    if user.user_id in upload.viewers:
        return True

    group_ids = MongoUserGroup.get_ids_by_user_id(user.user_id)
    if not set(group_ids).isdisjoint(upload.viewer_groups):
        return True

    return False


def is_user_upload_writer(upload: Upload, user: User):
    if user.is_admin:
        return True

    if user.user_id in upload.writers:
        return True

    group_ids = MongoUserGroup.get_ids_by_user_id(user.user_id)
    if not set(group_ids).isdisjoint(upload.writer_groups):
        return True

    return False


def get_upload_with_read_access(
    upload_id: str, user: User | None, include_others: bool = False
) -> Upload:
    """
    Determines if the user has read access to the upload. If so, the corresponding Upload
    object is returned. If the upload does not exist, or the user has no read access to
    it, an HTTPException is raised.

    Arguments:
        upload_id: The id of the requested upload.
        user: The authenticated user, if any.
        include_others: If uploads owned by others should be included. Access to the
        uploads of other users is only granted if the upload is published and not under
        embargo.
    """
    mongodb_query = _query_mongodb(upload_id=upload_id)
    upload = mongodb_query.first()
    if upload is None:
        raise HTTPException(
            status.HTTP_404_NOT_FOUND, detail='The specified upload_id was not found.'
        )

    if is_user_upload_viewer(upload, user):
        return upload

    if not include_others:
        raise HTTPException(
            status.HTTP_403_FORBIDDEN,
            detail='You do not have access to the specified upload.',
        )

    if not upload.published:
        if user is None:
            raise HTTPException(
                status.HTTP_401_UNAUTHORIZED,
                detail='You need to log in to access the specified upload.',
            )
        raise HTTPException(
            status.HTTP_403_FORBIDDEN,
            detail='You do not have access to the specified upload.',
        )

    if upload.with_embargo:
        raise HTTPException(
            status.HTTP_403_FORBIDDEN,
            detail='You do not have access to the specified upload - published with embargo.',
        )

    return upload


def _get_upload_with_write_access(
    upload_id: str,
    user: User,
    include_published: bool = False,
    published_requires_admin: bool = True,
    include_failed_imports: bool = False,
    only_main_author: bool = False,
) -> Upload:
    """
    Determines if the user has write access to the upload. If so, the corresponding Upload
    object is returned. If the upload does not exist, or the user has no write access to
    it, an HTTPException is raised.
    """
    if not user:
        raise HTTPException(
            status.HTTP_401_UNAUTHORIZED,
            detail='User authentication required to access uploads.',
        )

    mongodb_query = _query_mongodb(upload_id=upload_id)
    if not mongodb_query.count():
        raise HTTPException(
            status.HTTP_404_NOT_FOUND, detail='The specified upload_id was not found.'
        )
    upload = mongodb_query.first()

    if not is_user_upload_writer(upload, user):
        raise HTTPException(
            status.HTTP_403_FORBIDDEN,
            detail='You do not have write access to the specified upload.',
        )

    if only_main_author and not user.is_admin and upload.main_author != user.user_id:
        raise HTTPException(
            status.HTTP_403_FORBIDDEN,
            detail='Only main author has permissions for this operation.',
        )

    if not upload.published:
        return upload

    if not include_published:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='Upload is already published, operation not possible.',
        )

    is_failed_import = (
        upload.current_process
        and upload.current_process.startswith('import_bundle')
        and upload.process_status == ProcessStatus.FAILURE
    ) or upload.last_status_message == 'Import bundle failed'
    if (
        published_requires_admin
        and not user.is_admin
        and not (is_failed_import and include_failed_imports)
    ):
        raise HTTPException(
            status.HTTP_403_FORBIDDEN,
            detail='Upload is already published, only admins can perform this operation.',
        )

    return upload


def upload_to_pydantic(
    upload: Upload, *, include_total_count: bool = True
) -> UploadProcData:
    """Converts the mongo db object to an UploadProcData object."""
    pydantic_upload = UploadProcData.model_validate(upload)
    if include_total_count:
        pydantic_upload.entries = upload.total_entries_count
    try:
        pydantic_upload.upload_files_server_path = upload.upload_files.external_os_path
    except KeyError:
        # In case the files are missing for one reason or another
        pass

    return pydantic_upload


def entry_to_pydantic(
    entry: Entry, add_es_metadata: bool = False, user=None
) -> EntryProcData:
    """
    Converts the mongo db object to an EntryProcData object, and optionally also adds metadata
    from ES
    """
    rv = EntryProcData.model_validate(entry)
    if add_es_metadata:
        # load entries's metadata from search
        metadata_entries = search(
            pagination=MetadataPagination(page_size=1),
            owner='admin' if user.is_admin else 'visible',
            user_id=user.user_id,
            query=dict(entry_id=entry.entry_id),
        )
        if len(metadata_entries.data) == 1:
            rv.entry_metadata = metadata_entries.data[0]
    return rv


def _check_upload_not_processing(upload: Upload):
    """
    Checks if the upload is processing, and raises a HTTPException (err code 400) if so.
    """
    if upload.process_running:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='The upload is currently being processed, operation not allowed.',
        )


def _perform_status_check(url: str):
    response = requests.get(url, timeout=15)
    return response


def _check_external_deployment_status(deployment_url: str):
    parsed_url = urlparse(deployment_url)

    base_url = f'{parsed_url.scheme}://{parsed_url.netloc}'
    try:
        response = _perform_status_check(f'{base_url}/-/health')
        if response.status_code != status.HTTP_200_OK:
            raise HTTPException(
                status.HTTP_400_BAD_REQUEST,
                detail='The target deployment is not available or the URL is incorrectly formatted. The target deployment URL should end with /api.',
            )
    except requests.exceptions.ConnectionError:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='The target deployment is not available. Connection refused.',
        )
    except requests.exceptions.Timeout:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail='The target deployment is not available. Timeout.',
        )
    except Exception as e:
        raise HTTPException(
            status.HTTP_400_BAD_REQUEST,
            detail=f'Failed to check external deployment health. Error: {str(e)}',
        )


@router.post(
    '/{upload_id}/action/stop-processing',
    tags=[APITag.ACTION],
    summary='Stops the processing of the upload.',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _upload_not_found, _not_authorized_to_upload, _bad_request
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
def stop_upload_processing(
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_PROCESS], allow_anonymous=False)),
    ],
):
    """
    Stops the processing of the specified upload.
    """
    upload = _get_upload_with_write_access(upload_id, user, include_published=False)

    if upload.process_status != ProcessStatus.PENDING:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='This functionality is only available when upload process state is pending.',
        )
    upload.stop_processing()

    return UploadProcDataResponse(upload_id=upload_id, data=upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/assign-doi',
    tags=[APITag.ACTION],
    summary='Assign a DOI to an upload',
    response_model=UploadProcDataResponse,
    responses=create_responses(
        _datacite_not_enabled,
        _upload_not_found,
        _bad_request,
        _not_authorized_to_upload,
        _existing_upload_with_findable_state,
        _upload_already_has_doi,
        _upload_is_empty,
        _upload_is_unpublished,
        _datacite_draft_failed,
        _db_set_upload_doi_failed,
        _datacite_publish_failed,
    ),
    response_model_exclude_unset=True,
    response_model_exclude_none=True,
)
async def assign_doi(
    upload_id: Annotated[str, Path(description='The unique id of the upload.')],
    user: Annotated[
        User,
        Depends(get_current_user([Scope.UPLOADS_ASSIGN_DOI], allow_anonymous=False)),
    ],
):
    """
    Assign a DOI at DataCite to this upload.

    Conditions:

    - The DataCite service must be enabled on this deployment.
    - The upload must be published.
    - The upload must contain at least one entry.
    - The user must be the main author of the upload.
    """

    if not config.datacite.enabled:
        raise _create_exception(*_datacite_not_enabled)

    upload = _get_upload_with_write_access(
        upload_id,
        user,
        include_published=True,
        only_main_author=True,
        published_requires_admin=False,
    )

    if upload.doi is not None:
        raise _create_exception(*_upload_already_has_doi)

    if upload.total_entries_count == 0:
        raise _create_exception(*_upload_is_empty)

    if not upload.published:
        raise _create_exception(*_upload_is_unpublished)

    try:
        doi_id = create_doi_for_upload(upload)
    except DataCiteException:
        raise _create_exception(*_datacite_draft_failed)

    try:
        upload.doi = EmbeddedDOI(id=doi_id)
        upload.save()
    except MongoEngineException:
        raise _create_exception(*_db_set_upload_doi_failed)

    try:
        publish_doi(doi_id)
    except DataCiteException:
        raise _create_exception(*_datacite_publish_failed)

    return {'upload_id': upload.upload_id, 'data': upload}
