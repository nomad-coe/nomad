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

from werkzeug.exceptions import HTTPException
from flask import Flask, request, g, jsonify, send_file, Response
from flask_restful import Resource, Api, abort
from flask_cors import CORS
from flask_httpauth import HTTPBasicAuth
from elasticsearch.exceptions import NotFoundError
from datetime import datetime
import os.path
import zipstream
from zipfile import ZIP_DEFLATED
from contextlib import contextmanager

from nomad import config, infrastructure
from nomad.files import UploadFile, ArchiveFile, ArchiveLogFile
from nomad.utils import get_logger
from nomad.processing import Upload, NotAllowedDuringProcessing
from nomad.repo import RepoCalc
from nomad.user import User

base_path = config.services.api_base_path

app = Flask(
    __name__,
    static_url_path='%s/docs' % base_path,
    static_folder=os.path.abspath(os.path.join(os.path.dirname(__file__), '../docs/.build/html')))
CORS(app)

app.config['SECRET_KEY'] = config.services.api_secret

auth = HTTPBasicAuth()
api = Api(app)


@app.before_first_request
def setup():
    infrastructure.setup()


@auth.verify_password
def verify_password(username_or_token, password):
    # first try to authenticate by token
    user = User.verify_auth_token(username_or_token)
    if not user:
        # try to authenticate with username/password
        user = User.objects(email=username_or_token).first()
        if not user:
            g.user = None
            return True  # anonymous access

        if not user or not user.verify_password(password):
            return False

    g.user = user
    return True


def login_really_required(func):
    @auth.login_required
    def wrapper(*args, **kwargs):
        if g.user is None:
            abort(401, message='Anonymous access is forbidden, authorization required')
        else:
            return func(*args, **kwargs)
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


@app.route('/api/token')
@login_really_required
def get_auth_token():
    token = g.user.generate_auth_token(600)
    return jsonify({'token': token.decode('ascii'), 'duration': 600})


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


class RepoCalcRes(Resource):
    def get(self, upload_hash, calc_hash):
        """
        Get calculation data in repository form, which only entails the quanties shown
        in the repository. This is basically the elastic search index entry for the
        requested calculations. Calcs are references via *upload_hash*, *calc_hash*
        pairs.

        .. :quickref: repo; Get calculation data in repository form.

        **Example request**:

        .. sourcecode:: http

            GET /nomad/api/repo/W36aqCzAKxOCfIiMFsBJh3nHPb4a/7ddvtfRfZAvc3Crr7jOJ8UH0T34I HTTP/1.1
            Accept: application/json

        **Example response**:

        .. sourcecode:: http

            HTTP/1.1 200 OK
            Vary: Accept
            Content-Type: application/json

            {
                "calc_hash":"7ddvtfRfZAvc3Crr7jOJ8UH0T34I",
                "upload_time":"2018-08-30T08:41:51.771367",
                "upload_id":"5b87adb813a441000a70a968",
                "upload_hash":"W36aqCzAKxOCfIiMFsBJh3nHPb4a",
                "mainfile":"RopD3Mo8oMV_-E5bh8uW5PiiCRkH1/data/BrK_svSi/TFCC010.CAB/vasprun.xml.relax1",
                "program_name":"VASP",
                "program_version":"4.6.35  3Apr08 complex  parallel LinuxIFC",
                "chemical_composition_bulk_reduced":"BrKSi2",
                "program_basis_set_type":"plane waves",
                "atom_species":[
                    35,
                    19,
                    14,
                    14
                ],
                "system_type":"Bulk",
                "crystal_system":"orthorhombic",
                "space_group_number":47,
                "configuration_raw_gid":"sq6wTJjRKb2VTajoDLVWDxHCgyN6i",
                "XC_functional_name":"GGA_X_PBE"
            }

        :param string upload_hash: the hash of the upload (from uploaded file contents)
        :param string calc_hash: the hash of the calculation (from mainfile)
        :resheader Content-Type: application/json
        :status 200: calc successfully retrieved
        :status 404: calc with given hashes does not exist
        :returns: the repository calculation entry
        """
        try:
            return RepoCalc.get(id='%s/%s' % (upload_hash, calc_hash)).json_dict, 200
        except NotFoundError:
            abort(404, message='There is no calculation for %s/%s' % (upload_hash, calc_hash))
        except Exception as e:
            abort(500, message=str(e))


class RepoCalcsRes(Resource):
    @auth.login_required
    def get(self):
        """
        Get *'all'* calculations in repository from, paginated.

        .. :quickref: repo; Get *'all'* calculations in repository from, paginated.

        **Example request**:

        .. sourcecode:: http

            GET /nomad/api/repo?page=1&per_page=25 HTTP/1.1
            Accept: application/json

        **Example response**:

        .. sourcecode:: http

            HTTP/1.1 200 OK
            Vary: Accept
            Content-Type: application/json

            {
                "pagination":{
                    "total":1,
                    "page":1,
                    "per_page":25
                },
                "results":[
                    {
                        "calc_hash":"7ddvtfRfZAvc3Crr7jOJ8UH0T34I",
                        "upload_time":"2018-08-30T08:41:51.771367",
                        "upload_id":"5b87adb813a441000a70a968",
                        "upload_hash":"W36aqCzAKxOCfIiMFsBJh3nHPb4a",
                        "mainfile":"RopD3Mo8oMV_-E5bh8uW5PiiCRkH1/data/BrK_svSi/TFCC010.CAB/vasprun.xml.relax1",
                        "program_name":"VASP",
                        "program_version":"4.6.35  3Apr08 complex  parallel LinuxIFC",
                        "chemical_composition_bulk_reduced":"BrKSi2",
                        "program_basis_set_type":"plane waves",
                        "atom_species":[
                            35,
                            19,
                            14,
                            14
                        ],
                        "system_type":"Bulk",
                        "crystal_system":"orthorhombic",
                        "space_group_number":47,
                        "configuration_raw_gid":"sq6wTJjRKb2VTajoDLVWDxHCgyN6i",
                        "XC_functional_name":"GGA_X_PBE"
                    }
                ]
            }

        :qparam int page: the page starting with 1
        :qparam int per_page: desired calcs per page
        :qparam string owner: specifies which cals to return: all|user|staging, default is all
        :resheader Content-Type: application/json
        :status 200: calcs successfully retrieved
        :returns: a list of repository entries in ``results`` and pagination info
        """
        # TODO use argparse? bad request reponse an bad params, pagination as decorator
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 10))
        owner = request.args.get('owner', 'all')

        try:
            assert page >= 1
            assert per_page > 0
        except AssertionError:
            abort(400, message='invalid pagination')

        if owner == 'all':
            search = RepoCalc.search().query('match_all')
        elif owner == 'user':
            if g.user is None:
                abort(401, message='Authentication required for owner value user.')
            search = RepoCalc.search().query('match_all')
            search = search.filter('term', user_id=g.user.email)
        elif owner == 'staging':
            if g.user is None:
                abort(401, message='Authentication required for owner value user.')
            search = RepoCalc.search().query('match_all')
            search = search.filter('term', user_id=g.user.email).filter('term', staging=True)
        else:
            abort(400, message='Invalid owner value. Valid values are all|user|staging, default is all')

        search = search[(page - 1) * per_page: page * per_page]
        return {
            'pagination': {
                'total': search.count(),
                'page': page,
                'per_page': per_page
            },
            'results': [result.json_dict for result in search]
        }


@app.route('%s/logs/<string:upload_hash>/<string:calc_hash>' % base_path, methods=['GET'])
def get_calc_proc_log(upload_hash, calc_hash):
    """
    Get calculation processing log. Calcs are references via *upload_hash*, *calc_hash*
    pairs.

    .. :quickref: archive; Get calculation processing logs.

    **Example request**:

    .. sourcecode:: http

        GET /nomad/api/logs/W36aqCzAKxOCfIiMFsBJh3nHPb4a/7ddvtfRfZAvc3Crr7jOJ8UH0T34I HTTP/1.1
        Accept: application/json

    :param string upload_hash: the hash of the upload (from uploaded file contents)
    :param string calc_hash: the hash of the calculation (from mainfile)
    :resheader Content-Type: application/json
    :status 200: calc successfully retrieved
    :status 404: calc with given hashes does not exist
    :returns: the log data, a line by line sequence of structured logs
    """
    archive_id = '%s/%s' % (upload_hash, calc_hash)

    try:
        archive = ArchiveLogFile(archive_id)
        archive_path = archive.os_path

        rv = send_file(
            archive_path,
            mimetype='application/text',
            as_attachment=True,
            attachment_filename=os.path.basename(archive_path))

        return rv
    except KeyError:
        abort(404, message='Archive/calculation %s does not exist.' % archive_id)
    except FileNotFoundError:
        abort(404, message='Archive/calculation %s does not exist.' % archive_id)
    except Exception as e:
        logger = get_logger(
            __name__, endpoint='logs', action='get',
            upload_hash=upload_hash, calc_hash=calc_hash)
        logger.error('Exception on accessing calc proc log', exc_info=e)
        abort(500, message='Could not accessing the logs.')


@app.route('%s/archive/<string:upload_hash>/<string:calc_hash>' % base_path, methods=['GET'])
def get_calc(upload_hash, calc_hash):
    """
    Get calculation data in archive form. Calcs are references via *upload_hash*, *calc_hash*
    pairs.

    .. :quickref: archive; Get calculation data in archive form.

    **Example request**:

    .. sourcecode:: http

        GET /nomad/api/archive/W36aqCzAKxOCfIiMFsBJh3nHPb4a/7ddvtfRfZAvc3Crr7jOJ8UH0T34I HTTP/1.1
        Accept: application/json

    :param string upload_hash: the hash of the upload (from uploaded file contents)
    :param string calc_hash: the hash of the calculation (from mainfile)
    :resheader Content-Type: application/json
    :status 200: calc successfully retrieved
    :status 404: calc with given hashes does not exist
    :returns: the metainfo formated JSON data of the requested calculation
    """
    archive_id = '%s/%s' % (upload_hash, calc_hash)

    try:
        archive = ArchiveFile(archive_id)
        archive_path = archive.os_path

        rv = send_file(
            archive_path,
            mimetype='application/json',
            as_attachment=True,
            attachment_filename=os.path.basename(archive_path))

        if config.files.compress_archive:
            rv.headers['Content-Encoding'] = 'gzip'

        return rv
    except KeyError:
        abort(404, message='Archive %s does not exist.' % archive_id)
    except FileNotFoundError:
        abort(404, message='Archive %s does not exist.' % archive_id)
    except Exception as e:
        logger = get_logger(
            __name__, endpoint='archive', action='get',
            upload_hash=upload_hash, calc_hash=calc_hash)
        logger.error('Exception on accessing archive', exc_info=e)
        abort(500, message='Could not accessing the archive.')


@app.route('%s/raw/<string:upload_hash>/<string:calc_hash>' % base_path, methods=['GET'])
def get_raw(upload_hash, calc_hash):
    """
    Get calculation mainfile raw data. Calcs are references via *upload_hash*, *calc_hash*
    pairs. Returns the mainfile, unless an aux_file is specified. Aux files are stored
    in repository entries. See ``/repo`` endpoint.

    .. :quickref: repo; Get calculation raw data.

    **Example request**:

    .. sourcecode:: http

        GET /nomad/api/raw/W36aqCzAKxOCfIiMFsBJh3nHPb4a/7ddvtfRfZAvc3Crr7jOJ8UH0T34I HTTP/1.1
        Accept: application/gz

    :param string upload_hash: the hash of the upload (from uploaded file contents)
    :param string calc_hash: the hash of the calculation (from mainfile)
    :qparam str auxfile: an optional aux_file to download the respective aux file, default is mainfile
    :qparam all: set any value to get a .zip with main and aux files instead of an individual file
    :resheader Content-Type: application/json
    :status 200: calc raw data successfully retrieved
    :status 404: calc with given hashes does not exist or the given aux file does not exist
    :returns: the raw data in body
    """
    archive_id = '%s/%s' % (upload_hash, calc_hash)
    logger = get_logger(__name__, endpoint='raw', action='get', archive_id=archive_id)

    try:
        repo = RepoCalc.get(id=archive_id)
    except NotFoundError:
        abort(404, message='There is no calculation for %s/%s' % (upload_hash, calc_hash))
    except Exception as e:
        abort(500, message=str(e))

    @contextmanager
    def raw_file(filename):
        try:
            upload = Upload.get(repo.upload_id)
            upload_file = UploadFile(upload.upload_id, local_path=upload.local_path)
            the_file = upload_file.get_file(filename)
            with the_file.open() as f:
                yield f
        except KeyError:
            abort(404, message='The file %s does not exist.' % filename)
        except FileNotFoundError:
            abort(404, message='The file %s does not exist.' % filename)

    get_all = request.args.get('all', None) is not None
    if get_all:
        # retrieve the 'whole' calculation, meaning the mainfile and all aux files as
        # a .zip archive
        def generator():
            """ Stream a zip file with all files using zipstream. """
            def iterator():
                """ Replace the directory based iter of zipstream with an iter over all raw files. """
                def write(filename):
                    """ Write a raw file to the zipstream. """
                    def iter_content():
                        """ Iterate the raw file contents. """
                        with raw_file(filename) as file_object:
                            while True:
                                data = file_object.read(1024)
                                if not data:
                                    break
                                yield data
                    return dict(arcname=filename, iterable=iter_content())

                yield write(repo.mainfile)
                try:
                    for auxfile in repo.aux_files:
                        yield write(os.path.join(os.path.dirname(repo.mainfile), auxfile))
                except Exception as e:
                    logger.error('Exception while accessing auxfiles.', exc_info=e)

            zip_stream = zipstream.ZipFile(mode='w', compression=ZIP_DEFLATED)
            zip_stream.paths_to_write = iterator()

            for chunk in zip_stream:
                yield chunk

        response = Response(generator(), mimetype='application/zip')
        response.headers['Content-Disposition'] = 'attachment; filename={}'.format('%s.zip' % archive_id)
        return response
    else:
        # retrieve an individual raw file
        auxfile = request.args.get('auxfile', None)
        if auxfile:
            filename = os.path.join(os.path.dirname(repo.mainfile), auxfile)
        else:
            filename = repo.mainfile

        try:
            with raw_file(filename) as f:
                rv = send_file(
                    f,
                    mimetype='application/octet-stream',
                    as_attachment=True,
                    attachment_filename=os.path.basename(filename))
                return rv
        except HTTPException as e:
            raise e
        except Exception as e:
            logger = get_logger(
                __name__, endpoint='archive', action='get',
                upload_hash=upload_hash, calc_hash=calc_hash)
            logger.error('Exception on accessing archive', exc_info=e)
            abort(500, message='Could not accessing the archive.')


@app.route('%s/admin/<string:operation>' % base_path, methods=['POST'])
def call_admin_operation(operation):
    if operation == 'repair_uploads':
        Upload.repair_all()
    if operation == 'reset':
        infrastructure.reset()
    else:
        abort(400, message='Unknown operation %s' % operation)

    return 'done', 200


api.add_resource(UploadsRes, '%s/uploads' % base_path)
api.add_resource(UploadRes, '%s/uploads/<string:upload_id>' % base_path)
api.add_resource(UploadFileRes, '%s/uploads/<string:upload_id>/file' % base_path)
api.add_resource(RepoCalcsRes, '%s/repo' % base_path)
api.add_resource(RepoCalcRes, '%s/repo/<string:upload_hash>/<string:calc_hash>' % base_path)


if __name__ == '__main__':
    app.run(debug=True, port=8000)
