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
The upload API of the nomad@FAIRDI APIs. Provides endpoints to create uploads, upload
files, and retrieve the processing status of uploads.
"""

from flask import g, request
from flask_restplus import Resource, fields, abort
from datetime import datetime

from nomad.processing import Upload as UploadProc
from nomad.processing import NotAllowedDuringProcessing
from nomad.utils import get_logger
from nomad.files import UploadFile

from .app import api, base_path
from .auth import login_really_required
from .common import pagination_request_parser, pagination_model


ns = api.namespace(
    '%s/uploads' % base_path[1:] if base_path is not '' else 'uploads',
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
        description='The unique id for the upload. Its a random uuid and '
                    'and used within nomad as long as no upload_hash is available.'),
    'upload_hash': fields.String(
        description='The unique upload hash. It is based on the uploaded content and '
                    'used within nomad to identify uploads.'
    ),
    'additional_metadata': fields.Arbitrary,
    'upload_url': fields.String,
    'upload_command': fields.String,
    'local_path': fields.String,
    'upload_time': fields.DateTime(dt_format='iso8601'),
})

calc_model = api.inherit('UploadCalculationProcessing', proc_model, {
    'archive_id': fields.String,
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

upload_operation_model = api.model('UploadOperation', {
    'operation': fields.String(description='Currently unstage is the only operation.')
})


upload_metadata_parser = api.parser()
upload_metadata_parser.add_argument('name', type=str, help='An optional name for the upload.', location='args')
upload_metadata_parser.add_argument('local_path', type=str, help='Use a local file on the server.', location='args')


@ns.route('/')
class UploadList(Resource):
    @api.marshal_list_with(upload_model, skip_none=True, code=200, description='Uploads send')
    @login_really_required
    def get(self):
        """ Get the list of all uploads from the authenticated user. """
        return [upload for upload in UploadProc.user_uploads(g.user)], 200

    @api.marshal_list_with(upload_model, skip_none=True, code=200, description='Upload received')
    @api.expect(upload_metadata_parser)
    @login_really_required
    def put(self):
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
        # create upload
        upload = UploadProc.create(
            user=g.user,
            name=request.args.get('name'),
            local_path=local_path)

        logger = get_logger(__name__, endpoint='upload', action='put', upload_id=upload.upload_id)
        logger.info('upload created')

        uploadFile = UploadFile(upload.upload_id, local_path=local_path)

        if local_path:
            pass
        elif request.mimetype == 'application/multipart-formdata':
            # multipart formdata, e.g. with curl -X put "url" -F file=@local_file
            # might have performance issues for large files: https://github.com/pallets/flask/issues/2086
            if 'file' in request.files:
                abort(400, message='Bad multipart-formdata, there is no file part.')
            file = request.files['file']
            if upload.name is '':
                upload.name = file.filename

            file.save(uploadFile.os_path)
        else:
            # simple streaming data in HTTP body, e.g. with curl "url" -T local_file
            try:
                with uploadFile.open('wb') as f:
                    while not request.stream.is_exhausted:
                        f.write(request.stream.read(1024))

            except Exception as e:
                logger.error('Error on streaming upload', exc_info=e)
                abort(400, message='Some IO went wrong, download probably aborted/disrupted.')

        if not uploadFile.is_valid:
            uploadFile.delete()
            upload.delete()
            abort(400, message='Bad file format, excpected %s.' % ", ".join(UploadFile.formats))

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
class Upload(Resource):
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
            upload = UploadProc.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != str(g.user.user_id):
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

    @api.response(404, 'Upload does not exist')
    @api.response(400, 'Not allowed during processing or when not in staging')
    @api.marshal_with(upload_model, skip_none=True, code=200, description='Upload deleted')
    @login_really_required
    def delete(self, upload_id: str):
        """
        Delete an existing upload.

        Only ``is_ready`` uploads
        can be deleted. Deleting an upload in processing is not allowed.
        """
        try:
            upload = UploadProc.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != str(g.user.user_id):
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if not upload.in_staging:
            abort(400, message='Operation not allowed, upload is not in staging.')

        try:
            upload.delete()
            return upload, 200
        except NotAllowedDuringProcessing:
            abort(400, message='You must not delete an upload during processing.')

    @api.response(404, 'Upload does not exist or is not allowed')
    @api.response(400, 'Operation is not supported')
    @api.marshal_with(upload_model, skip_none=True, code=200, description='Upload unstaged successfully')
    @api.expect(upload_operation_model)
    @login_really_required
    def post(self, upload_id):
        """
        Execute an upload operation. Available operations: ``unstage``

        Untage changes the visibility of the upload. Clients can specify, if the calcs
        should be restricted.
        """
        try:
            upload = UploadProc.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != str(g.user.user_id):
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        json_data = request.get_json()
        if json_data is None:
            json_data = {}

        operation = json_data.get('operation')
        if operation == 'unstage':
            if not upload.in_staging:
                abort(400, message='Operation not allowed, upload is not in staging.')

            try:
                upload.unstage()
            except NotAllowedDuringProcessing:
                abort(400, message='You must not unstage an upload during processing.')

            return upload, 200

        abort(400, message='Unsuported operation %s.' % operation)
