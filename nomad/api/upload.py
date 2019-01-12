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

"""
The upload API of the nomad@FAIRDI APIs. Provides endpoints to upload files and
get the processing status of uploads.
"""

from flask import g, request
from flask_restplus import Resource, fields, abort
from datetime import datetime
from werkzeug.datastructures import FileStorage
import os.path

from nomad import config, utils
from nomad.processing import Upload
from nomad.processing import NotAllowedDuringProcessing
from nomad.files import ArchiveBasedStagingUploadFiles, StagingUploadFiles, UploadFiles

from .app import api, with_logger
from .auth import login_really_required
from .common import pagination_request_parser, pagination_model


ns = api.namespace(
    'uploads',
    description='Uploading data and tracing uploaded data and its processing.')


proc_model = api.model('Processing', {
    'tasks': fields.List(fields.String),
    'current_task': fields.String,
    'status': fields.String,
    'completed': fields.Boolean,
    'errors': fields.List(fields.String),
    'warnings': fields.List(fields.String),
    'create_time': fields.DateTime(dt_format='iso8601'),
    'complete_time': fields.DateTime(dt_format='iso8601'),
    '_async_status': fields.String(description='Only for debugging nomad')
})

upload_model = api.inherit('UploadProcessing', proc_model, {
    'name': fields.String(
        description='The name of the upload. This can be provided during upload '
                    'using the name query parameter.'),
    'upload_id': fields.String(
        description='The unique id for the upload.'),
    'additional_metadata': fields.Arbitrary,
    'local_path': fields.String,
    'upload_time': fields.DateTime(dt_format='iso8601'),
})

calc_model = api.inherit('UploadCalculationProcessing', proc_model, {
    'calc_id': fields.String,
    'mainfile': fields.String,
    'upload_id': fields.String,
    'parser': fields.String
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
        })),
        'results': fields.List(fields.Nested(model=calc_model))
    }))
})

meta_data_model = api.model('MetaData', {
    'with_embargo': fields.Boolean(default=False, description='Data with embargo is only visible to the upload until the embargo period ended.'),
    'comment': fields.List(fields.String, description='The comment are shown in the repository for each calculation.'),
    'references': fields.List(fields.String, descriptions='References allow to link calculations to external source, e.g. URLs.'),
    'coauthors': fields.List(fields.String, description='A list of co-authors given by user_id.'),
    'shared_with': fields.List(fields.String, description='A list of users to share calculations with given by user_id.'),
    '_upload_time': fields.List(fields.DateTime(dt_format='iso8601'), description='Overrride the upload time.'),
    '_uploader': fields.List(fields.String, description='Override the uploader with the given user id.')
})

calc_meta_data_model = api.inherit('CalcMetaData', meta_data_model, {
    'mainfile': fields.String(description='The calculation main output file is used to identify the calculation in the upload.'),
    '_checksum': fields.String(description='Override the calculation checksum'),
    '_pid': fields.String(description='Assign a specific pid. It must be unique.')
})

upload_meta_data_model = api.inherit('UploadMetaData', meta_data_model, {
    'calculations': fields.List(fields.Nested(model=calc_meta_data_model), description='Specific per calculation data that will override the upload data.')
})

upload_operation_model = api.model('UploadOperation', {
    'operation': fields.String(description='Currently unstage is the only operation.'),
    'metadata': fields.Nested(model=upload_meta_data_model, description='Additional upload and calculation meta data that should be considered for the operation')
})


upload_metadata_parser = api.parser()
upload_metadata_parser.add_argument('name', type=str, help='An optional name for the upload.', location='args')
upload_metadata_parser.add_argument('local_path', type=str, help='Use a local file on the server.', location='args')
upload_metadata_parser.add_argument('file', type=FileStorage, help='The file to upload.', location='files')


@ns.route('/')
class UploadListResource(Resource):
    @api.doc('get_uploads')
    @api.marshal_list_with(upload_model, skip_none=True, code=200, description='Uploads send')
    @login_really_required
    def get(self):
        """ Get the list of all uploads from the authenticated user. """
        return [upload for upload in Upload.user_uploads(g.user)], 200

    @api.doc('upload')
    @api.marshal_with(upload_model, skip_none=True, code=200, description='Upload received')
    @api.expect(upload_metadata_parser)
    @login_really_required
    @with_logger
    def put(self, logger):
        """
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
        """
        local_path = request.args.get('local_path')
        if local_path:
            if not os.path.exists(local_path):
                abort(404, message='The given local_path was not found.')

        # create upload
        upload = Upload.create(
            user=g.user,
            name=request.args.get('name'),
            local_path=local_path)

        logger.info('upload created', upload_id=upload.upload_id)

        upload_files = ArchiveBasedStagingUploadFiles(
            upload.upload_id, create=True, local_path=local_path)

        if local_path:
            # file is already there and does not to be received
            pass
        elif request.mimetype == 'application/multipart-formdata':
            # multipart formdata, e.g. with curl -X put "url" -F file=@local_file
            # might have performance issues for large files: https://github.com/pallets/flask/issues/2086
            if 'file' in request.files:
                abort(400, message='Bad multipart-formdata, there is no file part.')
            file = request.files['file']
            if upload.name is '':
                upload.name = file.filename

            file.save(upload_files.upload_file_os_path)
        else:
            # simple streaming data in HTTP body, e.g. with curl "url" -T local_file
            try:
                with open(upload_files.upload_file_os_path, 'wb') as f:
                    while not request.stream.is_exhausted:
                        f.write(request.stream.read(1024))

            except Exception as e:
                logger.warning('Error on streaming upload', exc_info=e)
                abort(400, message='Some IO went wrong, download probably aborted/disrupted.')

        if not upload_files.is_valid:
            upload_files.delete()
            upload.delete()
            logger.info('Invalid upload')
            abort(400, message='Bad file format, excpected %s.' % ", ".join(upload_files.formats))

        logger.info('received uploaded file')
        upload.upload_time = datetime.now()
        upload.process()
        logger.info('initiated processing')

        return upload, 200


class ProxyUpload:
    def __init__(self, upload, calcs):
        self.upload = upload
        self.calcs = calcs

    def __getattr__(self, name):
        return self.upload.__getattribute__(name)


@ns.route('/<string:upload_id>')
@api.doc(params={'upload_id': 'The unique id for the requested upload.'})
class UploadResource(Resource):
    @api.doc('get_upload')
    @api.response(404, 'Upload does not exist')
    @api.marshal_with(upload_with_calcs_model, skip_none=True, code=200, description='Upload send')
    @api.expect(pagination_request_parser)
    @login_really_required
    def get(self, upload_id: str):
        """
        Get an update for an existing upload.

        Will not only return the upload, but also its calculations paginated.
        Use the pagination params to determine the page.
        """
        try:
            upload = Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != str(g.user.user_id) and not g.user.is_admin:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        try:
            page = int(request.args.get('page', 1))
            per_page = int(request.args.get('per_page', 10))
            order_by = str(request.args.get('order_by', 'mainfile'))
            order = int(str(request.args.get('order', -1)))
        except Exception:
            abort(400, message='invalid pagination or ordering')

        try:
            assert page >= 1
            assert per_page > 0
        except AssertionError:
            abort(400, message='invalid pagination')

        if order_by not in ['mainfile', 'status', 'parser']:
            abort(400, message='invalid order_by field %s' % order_by)

        order_by = ('-%s' if order == -1 else '+%s') % order_by

        calcs = upload.all_calcs((page - 1) * per_page, page * per_page, order_by)
        failed_calcs = upload.failed_calcs
        result = ProxyUpload(upload, {
            'pagination': dict(
                total=upload.total_calcs, page=page, per_page=per_page,
                successes=upload.processed_calcs - failed_calcs, failures=failed_calcs),
            'results': [calc for calc in calcs]
        })

        return result, 200

    @api.doc('delete_upload')
    @api.response(404, 'Upload does not exist')
    @api.response(401, 'Upload does not belong to authenticated user.')
    @api.response(400, 'Not allowed during processing')
    @api.marshal_with(upload_model, skip_none=True, code=200, description='Upload deleted')
    @login_really_required
    @with_logger
    def delete(self, upload_id: str, logger):
        """
        Delete an existing upload.

        Only uploads that are sill in staging, not already delete, not still uploaded, and
        not currently processed, can be deleted.
        """
        try:
            upload = Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != str(g.user.user_id) and not g.user.is_admin:
            abort(401, message='Upload with id %s does not belong to you.' % upload_id)

        with utils.lnr(logger, 'delete processing upload'):
            try:
                upload.delete()
            except NotAllowedDuringProcessing:
                abort(400, message='You must not delete an upload during processing.')

        with utils.lnr(logger, 'delete upload files'):
            upload_files = UploadFiles.get(upload_id)
            assert upload_files is not None, 'Uploads existing in staging must have files.'
            if upload_files is not None:
                assert isinstance(upload_files, StagingUploadFiles), 'Uploads in staging must have staging files.'
                upload_files.delete()

        return upload, 200

    @api.doc('exec_upload_command')
    @api.response(404, 'Upload does not exist or not in staging')
    @api.response(400, 'Operation is not supported')
    @api.response(401, 'If the operation is not allowed for the current user')
    @api.marshal_with(upload_model, skip_none=True, code=200, description='Upload unstaged successfully')
    @api.expect(upload_operation_model)
    @login_really_required
    def post(self, upload_id):
        """
        Execute an upload operation. Available operations: ``unstage``

        Unstage accepts further meta data that allows to provide coauthors, comments,
        external references, etc. See the model for details. The fields that start with
        ``_underscore`` are only available for users with administrative priviledges.

        Unstage changes the visibility of the upload. Clients can specify the visibility
        via meta data.
        """
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

        meta_data = json_data.get('meta_data', {})
        for key in meta_data:
            if key.startswith('_'):
                if not g.user.is_admin:
                    abort(401, message='Only admin users can use _meta_data_keys.')
                break

        if operation == 'unstage':
            try:
                upload.unstage(meta_data)
            except NotAllowedDuringProcessing:
                abort(400, message='You must not unstage an upload during processing.')

            return upload, 200

        abort(400, message='Unsuported operation %s.' % operation)


upload_command_model = api.model('UploadCommand', {
    'upload_url': fields.Url,
    'upload_command': fields.String
})


@ns.route('/command')
class UploadCommandResource(Resource):
    @api.doc('get_upload_command')
    @api.marshal_with(upload_command_model, code=200, description='Upload command send')
    @login_really_required
    def get(self):
        """ Get url and example command for shell based uploads. """
        upload_url = 'http://%s:%s%s/uploads/' % (
            config.services.api_host,
            config.services.api_port,
            config.services.api_base_path)

        upload_command = 'curl -H "X-Token: %s" "%s" --upload-file <local_file>' % (
            g.user.get_auth_token().decode('utf-8'), upload_url)

        return dict(upload_url=upload_url, upload_command=upload_command), 200
