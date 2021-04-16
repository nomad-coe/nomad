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
from datetime import datetime
from typing import List, Dict, Any, Optional
from pydantic import BaseModel, Field, validator
from fastapi import (
    APIRouter, Request, File, UploadFile, status, Depends, Path, Query as FastApiQuery,
    HTTPException)

from nomad import utils, config, files, datamodel
from nomad.processing import Upload, ProcessAlreadyRunning, FAILURE
from nomad.processing.base import PROCESS_COMPLETED
from nomad.utils import strip

from .auth import get_required_user, get_required_user_bearer_or_upload_token
from ..models import (
    BaseModel, User, Direction, Pagination, PaginationResponse)
from ..utils import parameter_dependency_from_model

router = APIRouter()
default_tag = 'uploads'

logger = utils.get_logger(__name__)


class UploadPagination(Pagination):
    @validator('order_by')
    def validate_order_by(cls, order_by):  # pylint: disable=no-self-argument
        if order_by is None:
            return 'upload_time'  # Default value
        assert order_by in ('create_time', 'published'), 'order_by must be a valid attribute'
        return order_by

    @validator('page_after_value')
    def validate_page_after_value(cls, page_after_value, values):  # pylint: disable=no-self-argument
        # Validation handled elsewhere
        return page_after_value


upload_pagination_parameters = parameter_dependency_from_model(
    'upload_pagination_parameters', UploadPagination)


class UploadData(BaseModel):
    upload_id: str = Field(
        None,
        description='The unique id for the upload.')
    name: Optional[str] = Field(
        description='The name of the upload. This can be provided during upload '
                    'using the name query parameter.')
    create_time: datetime = Field(
        None,
        description='The time of creation.')
    upload_time: datetime = Field(
        None,
        description='The time of upload.')
    upload_path: Optional[str] = Field(
        description='Path to the uploaded file on the server.')
    published: bool = Field(
        False,
        description='If this upload is already published.')
    published_to: List[str] = Field(
        None,
        description='A list of other NOMAD deployments that this upload was uploaded to already.')
    tasks: List[str] = Field()
    current_task: str = Field()
    tasks_running: bool = Field()
    tasks_status: str = Field()
    errors: List[str] = Field()
    warnings: List[str] = Field()
    complete_time: Optional[datetime] = Field()
    current_process: str = Field()
    process_running: bool = Field()
    last_status_message: str = Field(
        None,
        description='The last informative message that the processing saved about this uploads status.')

    class Config:
        orm_mode = True


class UploadDataResponse(BaseModel):
    upload_id: str = Field(None, description=strip('''
        Unique id of the upload.'''))
    data: UploadData = Field(
        None, description=strip('''
        The upload data as a dictionary.'''))


class UploadQuery(BaseModel):
    upload_id: Optional[List[str]] = Field(
        description='Search for uploads matching the given id. Multiple values can be specified.')
    upload_name: Optional[List[str]] = Field(
        description='Search for uploads matching the given name. Multiple values can be specified.')
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


upload_query_parameters = parameter_dependency_from_model(
    'upload_query_parameters', UploadQuery)


class UploadQueryResponse(BaseModel):
    query: UploadQuery = Field()
    pagination: PaginationResponse = Field()
    data: List[UploadData] = Field(
        None, description=strip('''
        The upload metadata as a list. Each item is a dictionary with the data for each
        upload.'''))


@router.get(
    '', tags=[default_tag],
    summary='List uploads of authenticated user.',
    response_model=UploadQueryResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_uploads(
        request: Request,
        query: UploadQuery = Depends(upload_query_parameters),
        pagination: UploadPagination = Depends(upload_pagination_parameters),
        user: User = Depends(get_required_user)):
    '''
    Retrieves metadata about all uploads that match the given query criteria.
    '''
    # Build query
    query_kwargs: Dict[str, Any] = {}
    query_kwargs.update(user_id=str(user.user_id))

    if query.upload_id:
        query_kwargs.update(upload_id__in=query.upload_id)

    if query.upload_name:
        query_kwargs.update(name__in=query.upload_name)

    if query.is_processing is True:
        query_kwargs.update(process_status__ne=PROCESS_COMPLETED)
    elif query.is_processing is False:
        query_kwargs.update(process_status=PROCESS_COMPLETED)

    if query.is_published is True:
        query_kwargs.update(published=True)
    elif query.is_published is False:
        query_kwargs.update(published=False)

    # Fetch data from DB
    mongodb_query = _query_mongodb(**query_kwargs)
    # Create response
    start = pagination.get_simple_index()
    end = start + pagination.page_size

    order_by = pagination.order_by
    order_by_with_sign = order_by if pagination.order == Direction.asc else '-' + order_by
    mongodb_query = mongodb_query.order_by(order_by_with_sign, 'upload_id')

    data = [_upload_to_pydantic(upload) for upload in mongodb_query[start:end]]

    pagination_response = PaginationResponse(total=mongodb_query.count(), **pagination.dict())
    pagination_response.populate_simple_index_and_urls(request)

    return UploadQueryResponse(
        query=query,
        pagination=pagination_response,
        data=data)


@router.get(
    '/{upload_id}', tags=[default_tag],
    summary='Get a specific upload',
    response_model=UploadDataResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def get_uploads_id(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload to retrieve.'),
        user: User = Depends(get_required_user)):
    '''
    Fetches a specific upload by its upload_id.
    '''
    # Get upload (or throw exception if nonexistent/no access)
    upload = _get_upload_with_read_access(upload_id, user)

    return UploadDataResponse(
        upload_id=upload_id,
        data=_upload_to_pydantic(upload))


@router.post(
    '', tags=[default_tag],
    summary='Submit a new upload',
    response_model=UploadDataResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_uploads(
        request: Request,
        file: UploadFile = File(None),
        local_path: str = None,  # Internal use/admins only
        name: str = FastApiQuery(
            None,
            description=strip('''
            Specifies the name of the upload.''')),
        publish_directly: bool = FastApiQuery(
            False,
            description=strip('''
            If the upload should be published directly. False by default.''')),
        oasis_upload_id: str = FastApiQuery(
            None,
            description=strip('''
            For oasis uploads: the upload id of the oasis system.''')),
        oasis_uploader_id: str = FastApiQuery(
            None,
            description=strip('''
            For oasis uploads: the id of the user in the oasis system who made the upload
            originally. The user must also be registered in NOMAD.''')),
        oasis_deployment_id: str = FastApiQuery(
            None,
            description=strip('''
            For oasis uploads: the deployment id.''')),
        user: User = Depends(get_required_user_bearer_or_upload_token)):
    '''
    Upload a file to the repository. Can be used to upload files via browser or other
    http clients like curl. This will also start the processing of the upload.

    There are two basic ways to upload a file: streaming the file data in the http body or
    in the multipart-formdata. Both are supported. The first method does not transfer a
    filename, so it is recommended to supply the parameter `name` in this case.

    Method 1: streaming data

        curl -X 'POST' "url" -T local_file

    Method 2: multipart-formdata

        curl -X 'POST' "url" -F file=@local_file

    Authentication is required to perform an upload. This can either be done using the
    regular bearer token, or using the simplified upload token. To use the simplified
    upload token, just specify it as a query parameter in the url, i.e.

        curl -X 'POST' ".../uploads?token=ABC.XYZ" ...

    Note, there is a limit on how many unpublished uploads a user can have. If exceeded,
    error code 400 will be returned.
    '''
    # Determine the source data stream
    src_stream = None
    if local_path:
        # Local file - only for admins
        if not user.is_admin:
            raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
                You need to be admin to use local_path as method of upload.'''))
        if not os.path.exists(local_path) or not os.path.isfile(local_path):
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                The specified local_path cannot be found or is not a file.'''))
        if not name:
            name = os.path.basename(local_path)
        src_stream = _asyncronous_file_reader(open(local_path, 'rb'))
    elif file:
        # Data provided as formdata
        if not name:
            name = file.filename
        src_stream = _asyncronous_file_reader(file)
    else:
        # Data has to be sent streamed in the body
        src_stream = request.stream()

    # TODO: allow admin users to upload in the name of other users?

    # Check upload limit
    if not user.is_admin:
        if _query_mongodb(user_id=str(user.user_id), published=False).count() >= config.services.upload_limit:  # type: ignore
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                Limit of unpublished uploads exceeded for user.'''))

    # check if allowed to perform oasis upload
    from_oasis = oasis_upload_id is not None
    if from_oasis:
        if not user.is_oasis_admin:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail='Only an oasis admin can perform an oasis upload.')
        if oasis_uploader_id is None:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='You must provide the original uploader for an oasis upload.')
        if oasis_deployment_id is None:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='You must provide the oasis deployment id for an oasis upload.')
        if publish_directly is False:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='Oasis uploads are always published directly.')
        # Switch user!
        user = datamodel.User.get(user_id=oasis_uploader_id)
        if user is None:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='The given original uploader does not exist.')
    elif oasis_uploader_id is not None or oasis_deployment_id is not None:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='For an oasis upload you must provide an oasis_upload_id.')

    # Get upload_id and path
    if from_oasis:
        upload_id = oasis_upload_id
        try:
            Upload.get(upload_id)
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='An oasis upload with the given upload_id already exists.')
        except KeyError:
            pass
    else:
        upload_id = utils.create_uuid()

    logger.info('upload created', upload_id=upload_id)

    # Read the stream and save to file
    if local_path:
        upload_path = local_path  # Use the provided path
        uploaded_bytes = os.path.getsize(local_path)
    else:
        upload_path = files.PathObject(config.fs.tmp, upload_id).os_path
        try:
            with open(upload_path, 'wb') as f:
                uploaded_bytes = 0
                log_interval = 1e9
                log_unit = 'GB'
                next_log_at = log_interval
                async for chunk in src_stream:
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
        except Exception:
            if os.path.isfile(upload_path):
                os.remove(upload_path)
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                    Some IO went wrong, download probably aborted/disrupted.'''))

    if not uploaded_bytes:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
                Empty upload - not allowed.'''))

    logger.info('received uploaded file')

    upload = Upload.create(
        upload_id=upload_id,
        user=user,
        name=name,
        upload_time=datetime.utcnow(),
        upload_path=upload_path,
        temporary=local_path != upload_path,
        publish_directly=publish_directly or from_oasis,
        from_oasis=from_oasis,
        oasis_deployment_id=oasis_deployment_id)

    upload.process_upload()

    return UploadDataResponse(
        upload_id=upload_id,
        data=_upload_to_pydantic(upload))


@router.delete(
    '/{upload_id}', tags=[default_tag],
    summary='Delete an upload',
    response_model=UploadDataResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def delete_uploads_id(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload to delete.'),
        user: User = Depends(get_required_user)):
    '''
    Delete an existing upload.

    Only uploads that are sill in staging, not already deleted, not still uploaded, and
    not currently processed, can be deleted.
    '''
    upload = _get_upload_with_write_access(upload_id, user)

    if upload.tasks_running:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            The upload is not processed yet.'''))

    try:
        upload.delete_upload()
    except ProcessAlreadyRunning:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=strip('''
            The upload is still being processed.'''))
    except Exception as e:
        logger.error('could not delete processing upload', exc_info=e)
        raise

    return UploadDataResponse(
        upload_id=upload_id,
        data=_upload_to_pydantic(upload))


@router.post(
    '/{upload_id}/action/publish', tags=[default_tag],
    summary='Publish an upload',
    response_model=UploadDataResponse,
    response_model_exclude_unset=True,
    response_model_exclude_none=True)
async def post_uploads_id_action_publish(
        upload_id: str = Path(
            ...,
            description='The unique id of the upload to publish.'),
        with_embargo: bool = FastApiQuery(
            True,
            description='If the data is published with an embargo.'),
        embargo_length: int = FastApiQuery(
            36,
            description='Length of the requested embargo in months.'),
        user: User = Depends(get_required_user)):
    '''
    Publishes an upload. The upload cannot be modified after this point, and after the
    embargo period (if any) is expired, the generated archive entries will be publicly visible.
    '''
    upload = _get_upload_with_write_access(upload_id, user)

    metadata_dict: Dict[str, Any] = {'with_embargo': with_embargo}
    if with_embargo:
        if not embargo_length or not 0 < embargo_length <= 36:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail='embargo_length needs to be between 1 and 36 months.')
        metadata_dict.update(embargo_length=embargo_length)

    if upload.tasks_running:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='The upload is not processed yet.')
    if upload.tasks_status == FAILURE:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Cannot publish an upload that failed processing.')
    if upload.processed_calcs == 0:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='Cannot publish an upload without any resulting entries.')
    try:
        upload.compress_and_set_metadata(metadata_dict)
        upload.publish_upload()
    except ProcessAlreadyRunning:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail='The upload is still/already processed.')

    return UploadDataResponse(
        upload_id=upload_id,
        data=_upload_to_pydantic(upload))


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


def _get_upload_with_read_access(upload_id, user) -> Upload:
    '''
    Determines if the specified user has read access to the specified upload. If so, the
    corresponding Upload object is returned. If the upload does not exist, or the user has
    no read access to it, a HTTPException is raised.
    '''
    if not user:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            User authentication required to access uploads.'''))
    mongodb_query = _query_mongodb(upload_id=upload_id)
    if not mongodb_query.count():
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=strip('''
            The specified upload_id was not found.'''))
    upload = mongodb_query.first()
    if user.is_admin or upload.user_id == str(user.user_id):
        # Ok, it exists and belongs to user
        return upload
    else:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            You do not have access to the specified upload.'''))


def _get_upload_with_write_access(upload_id, user) -> Upload:
    if not user:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            User authentication required to access uploads.'''))
    mongodb_query = _query_mongodb(upload_id=upload_id)
    if not mongodb_query.count():
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=strip('''
            The specified upload_id was not found.'''))
    upload = mongodb_query.first()
    if upload.user_id != str(user.user_id) and not user.is_admin:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            You do not have write access to the specified upload.'''))
    if upload.published and not user.is_admin:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail=strip('''
            Only admins can change uploads that are published.'''))
    return upload


def _upload_to_pydantic(upload: Upload) -> UploadData:
    ''' Converts the mongo db object to a dictionary. '''
    return UploadData.from_orm(upload)
