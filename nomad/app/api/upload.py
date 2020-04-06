# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

'''
The upload API of the nomad@FAIRDI APIs. Provides endpoints to upload files and
get the processing status of uploads.
'''

from typing import Dict, Any
from flask import g, request, Response
from flask_restplus import Resource, fields, abort
from datetime import datetime
from werkzeug.datastructures import FileStorage
import os.path
import os
import io
from functools import wraps

from nomad import config, utils, files, search
from nomad.processing import Upload, FAILURE
from nomad.processing import ProcessAlreadyRunning
from nomad.app import common
from nomad.app.common import RFC3339DateTime

from .api import api
from .auth import authenticate, generate_upload_token
from .common import pagination_request_parser, pagination_model, upload_route, metadata_model


ns = api.namespace(
    'uploads',
    description='Uploading data and tracing uploaded data and its processing.')

proc_model = api.model('Processing', {
    'tasks': fields.List(fields.String),
    'current_task': fields.String,
    'tasks_running': fields.Boolean,
    'tasks_status': fields.String,
    'errors': fields.List(fields.String),
    'warnings': fields.List(fields.String),
    'create_time': RFC3339DateTime,
    'complete_time': RFC3339DateTime,
    'current_process': fields.String,
    'process_running': fields.Boolean,
})

calc_metadata_model = api.inherit('CalcMetaData', metadata_model, {
    'mainfile': fields.String(description='The calculation main output file is used to identify the calculation in the upload.'),
    '_pid': fields.String(description='Assign a specific pid. It must be unique.'),
    'external_id': fields.String(description='External user provided id. Does not have to be unique necessarily.')
})

upload_metadata_model = api.inherit('UploadMetaData', metadata_model, {
    'embargo_length': fields.Integer(description='Length of the requested embargo in months.'),
    'calculations': fields.List(fields.Nested(model=calc_metadata_model, skip_none=True), description='Specific per calculation data that will override the upload data.')
})

upload_model = api.inherit('UploadProcessing', proc_model, {
    'name': fields.String(
        description='The name of the upload. This can be provided during upload '
                    'using the name query parameter.'),
    'upload_id': fields.String(
        description='The unique id for the upload.'),
    # TODO just removed during migration, where this get particularily large
    # 'metadata': fields.Nested(model=upload_metadata_model, description='Additional upload and calculation meta data.', skip_none=True),
    'upload_path': fields.String(description='The uploaded file on the server'),
    'published': fields.Boolean(description='If this upload is already published'),
    'upload_time': RFC3339DateTime(),
})

upload_list_model = api.model('UploadList', {
    'pagination': fields.Nested(model=pagination_model, skip_none=True),
    'results': fields.List(fields.Nested(model=upload_model, skip_none=True))
})

calc_model = api.inherit('UploadCalculationProcessing', proc_model, {
    'calc_id': fields.String,
    'mainfile': fields.String,
    'upload_id': fields.String,
    'parser': fields.String,
    'metadata': fields.Raw(
        attribute='_entry_metadata',
        description='The repository metadata for this entry.')
})

upload_with_calcs_model = api.inherit('UploadWithPaginatedCalculations', upload_model, {
    'processed_calcs': fields.Integer,
    'total_calcs': fields.Integer,
    'failed_calcs': fields.Integer,
    'pending_calcs': fields.Integer,
    'calcs': fields.Nested(model=api.model('UploadPaginatedCalculations', {
        'pagination': fields.Nested(model=api.inherit('UploadCalculationPagination', pagination_model, {
            'successes': fields.Integer,
            'failures': fields.Integer,
        }), skip_none=True),
        'results': fields.List(fields.Nested(model=calc_model, skip_none=True))
    }), skip_none=True)
})

upload_operation_model = api.model('UploadOperation', {
    'operation': fields.String(description='Currently publish is the only operation.'),
    'metadata': fields.Nested(model=upload_metadata_model, description='Additional upload and calculation meta data. Will replace previously given metadata.')
})


upload_metadata_parser = api.parser()
upload_metadata_parser.add_argument('name', type=str, help='An optional name for the upload.', location='args')
upload_metadata_parser.add_argument('local_path', type=str, help='Use a local file on the server.', location='args')
upload_metadata_parser.add_argument('token', type=str, help='Upload token to authenticate with curl command.', location='args')
upload_metadata_parser.add_argument('file', type=FileStorage, help='The file to upload.', location='files')

upload_list_parser = pagination_request_parser.copy()
upload_list_parser.add_argument('state', type=str, help='List uploads with given state: all, unpublished, published.', location='args')
upload_list_parser.add_argument('name', type=str, help='Filter for uploads with the given name.', location='args')


def disable_marshalling(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except DisableMarshalling as e:
            return e.un_marshalled

    return wrapper


def marshal_with(*args, **kwargs):
    '''
    A special version of the RESTPlus marshal_with decorator that allows to disable
    marshalling at runtime by raising DisableMarshalling.
    '''
    def decorator(func):
        @api.marshal_with(*args, **kwargs)
        def with_marshalling(*args, **kwargs):
            return func(*args, **kwargs)

        @wraps(with_marshalling)
        def wrapper(*args, **kwargs):
            try:
                return with_marshalling(*args, **kwargs)
            except DisableMarshalling as e:
                return e.un_marshalled

        return wrapper
    return decorator


class DisableMarshalling(Exception):
    def __init__(self, body, status, headers):
        super().__init__()
        self.un_marshalled = Response(body, status=status, headers=headers)


@ns.route('/')
class UploadListResource(Resource):
    @api.doc('get_uploads')
    @api.response(400, 'Bad parameters')
    @api.marshal_with(upload_list_model, skip_none=True, code=200, description='Uploads send')
    @api.expect(upload_list_parser)
    @authenticate(required=True)
    def get(self):
        ''' Get the list of all uploads from the authenticated user. '''
        try:
            state = request.args.get('state', 'unpublished')
            name = request.args.get('name', None)
            page = int(request.args.get('page', 1))
            per_page = int(request.args.get('per_page', 10))
        except Exception:
            abort(400, message='bad parameter types')

        try:
            assert page >= 1
            assert per_page > 0
        except AssertionError:
            abort(400, message='invalid pagination')

        query_kwargs = {}
        if state == 'published':
            query_kwargs.update(published=True)
        elif state == 'unpublished':
            query_kwargs.update(published=False)
        elif state == 'all':
            pass
        else:
            abort(400, message='bad state value %s' % state)

        if name is not None:
            query_kwargs.update(name=name)

        uploads = Upload.user_uploads(g.user, **query_kwargs)
        total = uploads.count()

        results = [
            upload
            for upload in uploads.order_by('published', '-upload_time')[(page - 1) * per_page: page * per_page]]

        return dict(
            pagination=dict(total=total, page=page, per_page=per_page),
            results=results), 200

    @api.doc('upload')
    @api.expect(upload_metadata_parser)
    @api.response(400, 'To many uploads')
    @marshal_with(upload_model, skip_none=True, code=200, description='Upload received')
    @authenticate(required=True, upload_token=True)
    def put(self):
        '''
        Upload a file and automatically create a new upload in the process.
        Can be used to upload files via browser or other http clients like curl.
        This will also start the processing of the upload.

        There are two basic ways to upload a file: multipart-formdata or simply streaming
        the file data. Both are supported. The later one does not allow to transfer a
        filename or other meta-data. If a filename is available, it will become the
        name of the upload.

        Example commands:

            curl -X put ".../nomad/api/uploads/" -F file=@local_file
            curl ".../nomad/api/uploads/" --upload-file local_file

        There is a general limit on how many unpublished uploads a user can have. Will
        return 400 if this limit is exceeded.
        '''
        # check existence of local_path if local_path is used
        local_path = request.args.get('local_path')
        if local_path:
            if not os.path.exists(local_path):
                abort(404, message='The given local_path was not found.')

        # check the upload limit
        if not g.user.is_admin:
            if Upload.user_uploads(g.user, published=False).count() >= config.services.upload_limit:
                abort(400, 'Limit of unpublished uploads exceeded for user.')

        upload_name = request.args.get('name')
        upload_id = utils.create_uuid()

        logger = common.logger.bind(upload_id=upload_id, upload_name=upload_name)
        logger.info('upload created', )

        try:
            if local_path:
                # file is already there and does not to be received
                upload_path = local_path
            elif request.mimetype in ['multipart/form-data', 'application/multipart-formdata']:
                logger.info('receive upload as multipart formdata')
                upload_path = files.PathObject(config.fs.tmp, upload_id).os_path
                # multipart formdata, e.g. with curl -X put "url" -F file=@local_file
                # might have performance issues for large files: https://github.com/pallets/flask/issues/2086
                if 'file' not in request.files:
                    abort(400, message='Bad multipart-formdata, there is no file part.')
                file = request.files['file']
                if upload_name is None or upload_name is '':
                    upload_name = file.filename

                file.save(upload_path)
            else:
                # simple streaming data in HTTP body, e.g. with curl "url" -T local_file
                logger.info('started to receive upload streaming data')
                upload_path = files.PathObject(config.fs.tmp, upload_id).os_path

                try:
                    with open(upload_path, 'wb') as f:
                        received_data = 0
                        received_last = 0
                        while True:
                            data = request.stream.read(io.DEFAULT_BUFFER_SIZE)
                            if len(data) == 0:
                                break

                            received_data += len(data)
                            received_last += len(data)
                            if received_last > 1e9:
                                received_last = 0
                                # TODO remove this logging or reduce it to debug
                                logger.info('received streaming data', size=received_data)
                            f.write(data)

                except Exception as e:
                    logger.warning('Error on streaming upload', exc_info=e)
                    abort(400, message='Some IO went wrong, download probably aborted/disrupted.')
        except Exception as e:
            if not local_path and os.path.isfile(upload_path):
                os.remove(upload_path)
            logger.info('Invalid or aborted upload')
            raise e

        logger.info('received uploaded file')

        upload = Upload.create(
            upload_id=upload_id,
            user=g.user,
            name=upload_name,
            upload_time=datetime.utcnow(),
            upload_path=upload_path,
            temporary=local_path != upload_path)

        upload.process_upload()
        logger.info('initiated processing')

        if bool(request.args.get('token', False)) and request.headers.get('Accept', '') != 'application/json':
            raise DisableMarshalling(
                '''
Thanks for uploading your data to nomad.
Go back to %s and press reload to see the progress on your upload and publish your data.

''' % config.gui_url(),
                200, {'Content-Type': 'text/plain; charset=utf-8'})

        return upload, 200


class ProxyUpload:
    def __init__(self, upload, calcs):
        self.upload = upload
        self.calcs = calcs

    def __getattr__(self, name):
        return self.upload.__getattribute__(name)


@upload_route(ns)
class UploadResource(Resource):
    @api.doc('get_upload')
    @api.response(404, 'Upload does not exist')
    @api.response(400, 'Invalid parameters')
    @api.marshal_with(upload_with_calcs_model, skip_none=True, code=200, description='Upload send')
    @api.expect(pagination_request_parser)
    @authenticate(required=True)
    def get(self, upload_id: str):
        '''
        Get an update for an existing upload.

        Will not only return the upload, but also its calculations paginated.
        Use the pagination params to determine the page.
        '''
        try:
            upload = Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != str(g.user.user_id) and not g.user.is_admin:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        try:
            page = int(request.args.get('page', 1))
            per_page = int(request.args.get('per_page', 10))
            order_by = request.args.get('order_by', None)
            order = int(str(request.args.get('order', -1)))
        except Exception:
            abort(400, message='invalid pagination or ordering')

        try:
            assert page >= 1
            assert per_page > 0
        except AssertionError:
            abort(400, message='invalid pagination')

        if order_by is not None:
            order_by = str(order_by)
            if order_by not in ['mainfile', 'tasks_status', 'parser']:
                abort(400, message='invalid order_by field %s' % order_by)

            order_by = ('-%s' if order == -1 else '+%s') % order_by

        # load upload's calcs
        calcs = list(upload.all_calcs(
            (page - 1) * per_page, page * per_page, order_by=order_by))

        calc_ids = [calc.calc_id for calc in calcs]
        search_results = {
            hit['calc_id']: hit
            for hit in search.SearchRequest().search_parameter('calc_id', calc_ids).execute_scan()}

        for calc in calcs:
            calc._entry_metadata = search_results.get(calc.calc_id)

        failed_calcs = upload.failed_calcs
        result = ProxyUpload(upload, {
            'pagination': dict(
                total=upload.total_calcs, page=page, per_page=per_page,
                successes=upload.processed_calcs - failed_calcs, failures=failed_calcs),
            'results': calcs
        })

        return result, 200

    @api.doc('delete_upload')
    @api.response(404, 'Upload does not exist')
    @api.response(401, 'Upload does not belong to authenticated user.')
    @api.response(400, 'The upload is still/already processed')
    @api.marshal_with(upload_model, skip_none=True, code=200, description='Upload deleted')
    @authenticate(required=True)
    def delete(self, upload_id: str):
        '''
        Delete an existing upload.

        Only uploads that are sill in staging, not already deleted, not still uploaded, and
        not currently processed, can be deleted.
        '''
        try:
            upload = Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != str(g.user.user_id) and not g.user.is_admin:
            abort(401, message='Upload with id %s does not belong to you.' % upload_id)

        if upload.published and not g.user.is_admin:
            abort(400, message='The upload is already published')

        if upload.tasks_running:
            abort(400, message='The upload is not processed yet')

        try:
            upload.delete_upload()
        except ProcessAlreadyRunning:
            abort(400, message='The upload is still processed')
        except Exception as e:
            common.logger.error('could not delete processing upload', exc_info=e)
            raise e

        return upload, 200

    @api.doc('exec_upload_operation')
    @api.response(404, 'Upload does not exist or not in staging')
    @api.response(400, 'Operation is not supported or the upload is still/already processed')
    @api.response(401, 'If the operation is not allowed for the current user')
    @api.marshal_with(upload_model, skip_none=True, code=200, description='Upload published successfully')
    @api.expect(upload_operation_model)
    @authenticate(required=True)
    def post(self, upload_id):
        '''
        Execute an upload operation. Available operations are ``publish`` and ``re-process``

        Publish accepts further meta data that allows to provide coauthors, comments,
        external references, etc. See the model for details. The fields that start with
        ``_underscore`` are only available for users with administrative privileges.

        Publish changes the visibility of the upload. Clients can specify the visibility
        via meta data.

        Re-process will re-process the upload and produce updated repository metadata and
        archive. Only published uploads that are not processing at the moment are allowed.
        Only for uploads where calculations have been processed with an older nomad version.
        '''
        try:
            upload = Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != str(g.user.user_id) and not g.user.is_admin:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        json_data = request.get_json()
        if json_data is None:
            json_data = {}

        operation = json_data.get('operation')

        user_metadata: Dict[str, Any] = json_data.get('metadata', {})
        metadata: Dict[str, Any] = {}
        for user_key in user_metadata:
            if user_key.startswith('_'):
                if not g.user.is_admin:
                    abort(401, message='Only admin users can use _metadata_keys.')

                key = user_key[1:]
            else:
                key = user_key

            metadata[key] = user_metadata[user_key]

        if operation == 'publish':
            if upload.tasks_running:
                abort(400, message='The upload is not processed yet')
            if upload.tasks_status == FAILURE:
                abort(400, message='Cannot publish an upload that failed processing')
            if upload.processed_calcs == 0:
                abort(400, message='Cannot publish an upload without calculations')
            try:
                upload.compress_and_set_metadata(metadata)
                upload.publish_upload()
            except ProcessAlreadyRunning:
                abort(400, message='The upload is still/already processed')

            return upload, 200
        elif operation == 're-process':
            if upload.tasks_running or upload.process_running or not upload.published:
                abort(400, message='Can only non processing, re-process published uploads')

            if len(metadata) > 0:
                abort(400, message='You can not provide metadata for re-processing')

            if len(upload.outdated_calcs) == 0:
                abort(400, message='You can only re-process uploads with at least one outdated calculation')

            upload.reset()
            upload.re_process_upload()

            return upload, 200

        abort(400, message='Unsupported operation %s.' % operation)


upload_command_model = api.model('UploadCommand', {
    'upload_url': fields.Url,
    'upload_command': fields.String,
    'upload_command_with_name': fields.String,
    'upload_progress_command': fields.String,
    'upload_command_form': fields.String,
    'upload_tar_command': fields.String
})


@ns.route('/command')
class UploadCommandResource(Resource):
    @api.doc('get_upload_command')
    @api.marshal_with(upload_command_model, code=200, description='Upload command send')
    @authenticate(required=True)
    def get(self):
        ''' Get url and example command for shell based uploads. '''
        token = generate_upload_token(g.user)
        upload_url = '%s/uploads/?token=%s' % (config.api_url(ssl=False), token)
        upload_url_with_name = upload_url + '&name=<name>'

        # upload_command = 'curl -X PUT "%s" -F file=@<local_file>' % upload_url

        # Upload via streaming data tends to work much easier, e.g. no mime type issues, etc.
        # It is also easier for the user to unterstand IMHO.
        upload_command = 'curl "%s" -T <local_file>' % upload_url

        upload_command_form = 'curl "%s" -X PUT -F file=@<local_file>' % upload_url

        upload_command_with_name = 'curl "%s" -X PUT -T <local_file>' % upload_url_with_name

        upload_progress_command = upload_command + ' | xargs echo'
        upload_tar_command = 'tar -cf - <local_folder> | curl -# -H "%s" -T - | xargs echo' % upload_url

        return dict(
            upload_url=upload_url,
            upload_command=upload_command,
            upload_command_with_name=upload_command_with_name,
            upload_progress_command=upload_progress_command,
            upload_command_form=upload_command_form,
            upload_tar_command=upload_tar_command), 200
