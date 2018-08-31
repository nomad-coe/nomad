from flask import Flask, request, redirect
from flask_restful import Resource, Api, abort
from flask_cors import CORS
from elasticsearch.exceptions import NotFoundError

from nomad import config, files
from nomad.utils import get_logger
from nomad.data import Calc, Upload, User, InvalidId, NotAllowedDuringProcessing, me

base_path = config.services.api_base_path

app = Flask(
    __name__,
    static_url_path='%s/docs' % base_path,
    static_folder='../docs/.build/html')
CORS(app)
api = Api(app)


class UploadsRes(Resource):
    """ Uploads """
    def get(self):
        """
        Get a list of current users uploads.

        .. :quickref: Get a list of current users uploads.

        **Example request**:

        .. sourcecode:: http

            GET /nomadxt/api/uploads HTTP/1.1
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
                    "is_ready": true,
                    "proc": {
                        "task_names": [
                            "uploading",
                            "extracting",
                            "parse_all",
                            "cleanup"
                        ],
                        "current_task_name": "cleanup",
                        "status": "SUCCESS",
                        "errors": [],
                        "warnings": [],
                        "upload_id": "5b89469e0d80d40008077dbc",
                        "upload_hash": "rMB5F-gyHT0KY22eePoTjXibK95S",
                        "calc_procs": []
                    }
                }
            ]

        :resheader Content-Type: application/json
        :status 200: uploads successfully provided
        :returns: list of :class:`nomad.data.Upload`
        """
        return [upload.json_dict for upload in Upload.user_uploads(me)], 200

    def post(self):
        """
        Create a new upload. Creating an upload on its own wont do much, but provide
        a *presigned* upload URL. PUT a file to this URL to do the actual upload and
        initiate the processing.

        .. :quickref: Create a new upload.

        **Example request**:

        .. sourcecode:: http

            POST /nomadxt/api/uploads HTTP/1.1
            Accept: application/json
            Content-Type: application/json

            {
                name: 'vasp_data.zip'
            }

        **Example response**:

        .. sourcecode:: http

            HTTP/1.1 200 OK
            Vary: Accept
            Content-Type: application/json

            [
                {
                    "name": "vasp_data.zip",
                    "upload_id": "5b89469e0d80d40008077dbc",
                    "presigned_url": "http://minio:9000/uploads/5b89469e0d80d40008077dbc?X-Amz-Algorithm=AWS4-...",
                    "create_time": "2018-08-31T13:46:06.781000",
                    "is_stale": false,
                    "is_ready": false,
                    "proc": {
                        "task_names": [
                            "uploading",
                            "extracting",
                            "parse_all",
                            "cleanup"
                        ],
                        "current_task_name": "uploading",
                        "status": "PENDING",
                        "errors": [],
                        "warnings": [],
                        "upload_id": "5b89469e0d80d40008077dbc",
                        "calc_procs": []
                    }
                }
            ]
        :jsonparam string name: An optional name for the upload.
        :reqheader Content-Type: application/json
        :resheader Content-Type: application/json
        :status 200: upload successfully created
        :returns: a new instance of :class:`nomad.data.Upload`
        """
        json_data = request.get_json()
        if json_data is None:
            json_data = {}

        return Upload.create(user=me, name=json_data.get('name', None)).json_dict, 200


class UploadRes(Resource):
    """ Uploads """
    def get(self, upload_id):
        """
        Get an update for an existing upload. If the upload will be upaded with new data from the
        processing infrastucture. You can use this endpoint to periodically pull the
        new processing state until the upload ``is_ready``.

        .. :quickref: Get an update for an existing upload.

        **Example request**:

        .. sourcecode:: http

            GET /nomadxt/api/uploads/5b89469e0d80d40008077dbc HTTP/1.1
            Accept: application/json

        **Example response**:

        .. sourcecode:: http

            HTTP/1.1 200 OK
            Vary: Accept
            Content-Type: application/json

            [
                {
                    "name": "vasp_data.zip",
                    "upload_id": "5b89469e0d80d40008077dbc",
                    "presigned_url": "http://minio:9000/uploads/5b89469e0d80d40008077dbc?X-Amz-Algorithm=AWS4-...",
                    "create_time": "2018-08-31T13:46:06.781000",
                    "upload_time": "2018-08-31T13:46:16.824000",
                    "is_stale": false,
                    "is_ready": false,
                    "proc": {
                        "task_names": [
                            "uploading",
                            "extracting",
                            "parse_all",
                            "cleanup"
                        ],
                        "current_task_name": "extracting",
                        "status": "PROGRESS",
                        "errors": [],
                        "warnings": [],
                        "upload_id": "5b89469e0d80d40008077dbc",
                        "calc_procs": []
                    }
                }
            ]
        :param string upload_id: the id for the upload
        :resheader Content-Type: application/json
        :status 200: upload successfully updated and retrieved
        :status 400: bad upload id
        :status 404: upload with id does not exist
        :returns: the :class:`nomad.data.Upload` instance
        """
        try:
            return Upload.get(upload_id=upload_id).json_dict, 200
        except InvalidId:
            abort(400, message='%s is not a valid upload id.' % upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

    def delete(self, upload_id):
        """
        Deletes an existing upload. Only ``is_ready`` or ``is_stale`` uploads
        can be deleted. Deleting an upload in processing is not allowed.

        .. :quickref: Delete an existing upload.

        **Example request**:

        .. sourcecode:: http

            DELETE /nomadxt/api/uploads/5b89469e0d80d40008077dbc HTTP/1.1
            Accept: application/json

        :param string upload_id: the id for the upload
        :resheader Content-Type: application/json
        :status 200: upload successfully deleted
        :status 400: upload cannot be deleted
        :status 404: upload with id does not exist
        :returns: the :class:`nomad.data.Upload` instance with the latest processing state
        """
        try:
            return Upload.get(upload_id=upload_id).delete().json_dict, 200
        except InvalidId:
            abort(400, message='%s is not a valid upload id.' % upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)
        except NotAllowedDuringProcessing:
            abort(400, message='You must not delete an upload during processing.')


class RepoCalcRes(Resource):
    def get(self, upload_hash, calc_hash):
        try:
            return Calc.get(id='%s/%s' % (upload_hash, calc_hash)).json_dict, 200
        except NotFoundError:
            abort(404, message='There is no calculation for %s/%s' % (upload_hash, calc_hash))
        except Exception as e:
            abort(500, message=str(e))


class RepoCalcsRes(Resource):
    def get(self):
        logger = get_logger(__name__, endpoint='repo', action='get')

        # TODO use argparse? bad request reponse an bad params, pagination as decorator
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 10))

        assert page >= 1
        assert per_page > 0

        try:
            search = Calc.search().query('match_all')
            search = search[(page - 1) * per_page: page * per_page]
            return {
                'pagination': {
                    'total': search.count(),
                    'page': page,
                    'per_page': per_page
                },
                'results': [result.json_dict for result in search]
            }
        except Exception as e:
            logger.error('Could not execute repo calcs get', exc_info=e)
            abort(500, message=str(e))


@app.route('%s/archive/<string:upload_hash>/<string:calc_hash>' % base_path, methods=['GET'])
def get_calc(upload_hash, calc_hash):
    logger = get_logger(__name__, endpoint='archive', action='get', upload_hash=upload_hash, calc_hash=calc_hash)

    archive_id = '%s/%s' % (upload_hash, calc_hash)

    try:
        url = files.external_objects_url(files.archive_url(archive_id))
        return redirect(url, 302)
    except KeyError:
        abort(404, message='Archive %s does not exist.' % archive_id)
    except Exception as e:
        logger.error('Exception on accessing archive', exc_info=e)
        abort(500, message='Could not accessing the archive.')


api.add_resource(UploadsRes, '%s/uploads' % base_path)
api.add_resource(UploadRes, '%s/uploads/<string:upload_id>' % base_path)
api.add_resource(RepoCalcsRes, '%s/repo' % base_path)
api.add_resource(RepoCalcRes, '%s/repo/<string:upload_hash>/<string:calc_hash>' % base_path)


if __name__ == '__main__':
    app.run(debug=True, port=8000)
