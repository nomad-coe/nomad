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

from datetime import datetime

from flask import g, request
from flask_restful import Resource, abort

from nomad.files import UploadFile
from nomad.processing import NotAllowedDuringProcessing, Upload
from nomad.utils import get_logger

from .app import api, base_path, login_really_required

"""
The upload API of the nomad@FAIRDI APIs. Provides endpoints to create uploads, upload
files, and retrieve the processing status of uploads.
"""


class UploadsRes(Resource):
    """ Uploads """
    @login_really_required
    def get(self):
        """
        Get a list of current users uploads.

        .. :quickref: upload; Get a list of current users uploads.

        **Example request**:

        .. sourcecode:: http

            GET /nomad/api/uploads HTTP/1.1
            Accept: application/json

        **Example response**:

        .. sourcecode:: http

            HTTP/1.1 200 OK
            Vary: Accept
            Content-Type: application/json

            [
                {
                    "name": "examples_vasp_6.zip",
                    "upload_id": "5b89469e0d80d40008077dbc",
                    "presigned_url": "http://minio:9000/uploads/5b89469e0d80d40008077dbc?X-Amz-Algorithm=AWS4-...",
                    "create_time": "2018-08-31T13:46:06.781000",
                    "upload_time": "2018-08-31T13:46:07.531000",
                    "is_stale": false,
                    "completed": true,
                    "status": "SUCCESS",
                    "current_task": "cleanup",
                    "tasks": ["uploading", "extracting", "parse_all", "cleanup"]
                    "errors": [],
                    "warnings": []
                }
            ]

        :resheader Content-Type: application/json
        :status 200: uploads successfully provided
        :returns: list of :class:`nomad.data.Upload`
        """
        return [upload.json_dict for upload in Upload.user_uploads(g.user)], 200

    @login_really_required
    def post(self):
        """
        Create a new upload. Creating an upload on its own wont do much, but provide
        a *presigned* upload URL. PUT a file to this URL to do the actual upload and
        initiate the processing.

        .. :quickref: upload; Create a new upload.

        **Example request**:

        .. sourcecode:: http

            POST /nomad/api/uploads HTTP/1.1
            Accept: application/json
            Content-Type: application/json

            {
                "name": "vasp_data.zip"
            }

        **Example response**:

        .. sourcecode:: http

            HTTP/1.1 200 OK
            Vary: Accept
            Content-Type: application/json

            {
                "name": "vasp_data.zip",
                "upload_id": "5b89469e0d80d40008077dbc",
                "presigned_url": "http://minio:9000/uploads/5b89469e0d80d40008077dbc?X-Amz-Algorithm=AWS4-...",
                "create_time": "2018-08-31T13:46:06.781000",
                "upload_time": "2018-08-31T13:46:07.531000",
                "is_stale": false,
                "completed": true,
                "status": "SUCCESS",
                "current_task": "cleanup",
                "tasks": ["uploading", "extracting", "parse_all", "cleanup"]
                "errors": [],
                "warnings": [],
                "calcs": [
                    {
                        "current_task": "archiving",
                        "tasks": ["parsing", "normalizing", "archiving"]
                        "status": "SUCCESS",
                        "errors": [],
                        "warnings": [],
                        "parser": "parsers/vasp",
                        "mainfile": "Si.xml"
                    }
                ]
            }

        :jsonparam string name: An optional name for the upload.
        :jsonparem string local_path: An optional path the a file that is already on the server.
            In this case, uploading a file won't be possible, the local file is processed
            immediatly as if it was uploaded.
        :reqheader Content-Type: application/json
        :resheader Content-Type: application/json
        :status 200: upload successfully created
        :returns: a new instance of :class:`nomad.data.Upload`
        """
        json_data = request.get_json()
        if json_data is None:
            json_data = {}

        upload = Upload.create(
            user=g.user,
            name=json_data.get('name'),
            local_path=json_data.get('local_path'))

        if upload.local_path is not None:
            logger = get_logger(
                __name__, endpoint='uploads', action='post', upload_id=upload.upload_id)
            logger.info('file uploaded offline')
            upload.upload_time = datetime.now()
            upload.process()
            logger.info('initiated processing')

        return upload.json_dict, 200


class UploadRes(Resource):
    """ Uploads """
    @login_really_required
    def get(self, upload_id):
        """
        Get an update on an existing upload. Will not only return the upload, but
        also its calculations paginated. Use the pagination params to determine
        the page.

        .. :quickref: upload; Get an update for an existing upload.

        **Example request**:

        .. sourcecode:: http

            GET /nomad/api/uploads/5b89469e0d80d40008077dbc HTTP/1.1
            Accept: application/json

        **Example response**:

        .. sourcecode:: http

            HTTP/1.1 200 OK
            Vary: Accept
            Content-Type: application/json

            {
                "name": "vasp_data.zip",
                "upload_id": "5b89469e0d80d40008077dbc",
                "presigned_url": "http://minio:9000/uploads/5b89469e0d80d40008077dbc?X-Amz-Algorithm=AWS4-...",
                "create_time": "2018-08-31T13:46:06.781000",
                "upload_time": "2018-08-31T13:46:07.531000",
                "is_stale": false,
                "completed": true,
                "status": "SUCCESS",
                "current_task": "cleanup",
                "tasks": ["uploading", "extracting", "parse_all", "cleanup"]
                "errors": [],
                "warnings": [],
                "calcs": {
                    "pagination": {
                        "total": 1,
                        "page": 1,
                        "per_page": 25
                    },
                    "results": [
                        {
                            "current_task": "archiving",
                            "tasks": ["parsing", "normalizing", "archiving"]
                            "status": "SUCCESS",
                            "errors": [],
                            "warnings": [],
                            "parser": "parsers/vasp",
                            "mainfile": "Si.xml"
                        }
                    ]
                }
            }

        :param string upload_id: the id for the upload
        :qparam int page: the page starting with 1
        :qparam int per_page: desired calcs per page
        :qparam str order_by: the field to sort the calcs by, use [status,mainfile]
        :resheader Content-Type: application/json
        :status 200: upload successfully updated and retrieved
        :status 404: upload with id does not exist
        :returns: the :class:`nomad.data.Upload` instance
        """
        try:
            upload = Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != g.user.email:
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
        result = upload.json_dict
        result['calcs'] = {
            'pagination': dict(
                total=upload.total_calcs, page=page, per_page=per_page,
                successes=upload.processed_calcs - failed_calcs, failures=failed_calcs),
            'results': [calc.json_dict for calc in calcs]
        }

        return result, 200

    @login_really_required
    def post(self, upload_id):
        """
        Move an upload out of the staging area. This changes the visibility of the upload.
        Clients can specify, if the calcs should be restricted.

        .. :quickref: upload; Move an upload out of the staging area.

        **Example request**:

        .. sourcecode:: http

            POST /nomad/api/uploads HTTP/1.1
            Accept: application/json
            Content-Type: application/json

            {
                "operation": "unstage"
            }


        :param string upload_id: the upload id
        :resheader Content-Type: application/json
        :status 200: upload unstaged successfully
        :status 404: upload could not be found
        :status 400: if the operation is not supported
        :returns: the upload record
        """
        try:
            upload = Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != g.user.email:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        json_data = request.get_json()
        if json_data is None:
            json_data = {}

        operation = json_data.get('operation')
        if operation == 'unstage':
            upload.unstage()
            return upload.json_dict, 200

        abort(400, message='Unsuported operation %s.' % operation)

    @login_really_required
    def delete(self, upload_id):
        """
        Deletes an existing upload. Only ``is_ready`` or ``is_stale`` uploads
        can be deleted. Deleting an upload in processing is not allowed.

        .. :quickref: upload; Delete an existing upload.

        **Example request**:

        .. sourcecode:: http

            DELETE /nomad/api/uploads/5b89469e0d80d40008077dbc HTTP/1.1
            Accept: application/json

        :param string upload_id: the id for the upload
        :resheader Content-Type: application/json
        :status 200: upload successfully deleted
        :status 400: upload cannot be deleted
        :status 404: upload with id does not exist
        :returns: the :class:`nomad.data.Upload` instance with the latest processing state
        """
        try:
            upload = Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.user_id != g.user.email:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        try:
            upload.delete()
            return upload.json_dict, 200
        except NotAllowedDuringProcessing:
            abort(400, message='You must not delete an upload during processing.')


class UploadFileRes(Resource):
    """
    Upload a file to an existing upload. Can be used to upload files via bowser
    or other http clients like curl. This will start the processing of the upload.

    There are two basic ways to upload a file: multipart-formdata or simply streaming
    the file data. Both are supported. The later one does not allow to transfer a
    filename or other meta-data. If a filename is available, it will become the
    name of the upload.

    .. :quickref: upload; Upload a file to an existing upload.

    **Curl examples for both approaches**:

    .. sourcecode:: sh

        curl -X put "/nomad/api/uploads/5b89469e0d80d40008077dbc/file" -F file=@local_file
        curl "/nomad/api/uploads/5b89469e0d80d40008077dbc/file" --upload-file local_file

    :param string upload_id: the upload_id of the upload
    :resheader Content-Type: application/json
    :status 200: upload successfully received.
    :status 404: upload with given id does not exist
    :status 400: if the fileformat is not supported or the form data is different than expected.
    :returns: the upload (see GET /uploads/<upload_id>)
    """
    @login_really_required
    def put(self, upload_id):
        logger = get_logger(__name__, endpoint='upload', action='put', upload_id=upload_id)

        try:
            upload = Upload.get(upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        if upload.upload_time is not None:
            abort(400, message='A file was already uploaded to this uploade before.')

        uploadFile = UploadFile(upload_id)

        if request.mimetype == 'application/multipart-formdata':
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
            abort(400, message='Bad file format, excpected %s.' % ", ".join(UploadFile.formats))

        logger.info('received uploaded file')
        upload.upload_time = datetime.now()
        upload.process()
        logger.info('initiated processing')

        return upload.json_dict, 200


api.add_resource(UploadsRes, '%s/uploads' % base_path)
api.add_resource(UploadRes, '%s/uploads/<string:upload_id>' % base_path)
api.add_resource(UploadFileRes, '%s/uploads/<string:upload_id>/file' % base_path)
