from typing import Tuple
from flask import Flask, request, redirect
from flask_restful import Resource, Api, abort
from datetime import datetime
import mongoengine.errors
from flask_cors import CORS
import logging
from elasticsearch.exceptions import NotFoundError

from nomad import users, files, search, config
from nomad.processing import UploadProc
from nomad.utils import get_logger

base_path = config.services.api_base_path

app = Flask(
    __name__,
    static_url_path='%s/docs' % base_path,
    static_folder='../docs/.build/html')
CORS(app)
api = Api(app)


# provid a fake user for testing
me = users.User.objects(email='me@gmail.com').first()
if me is None:
    me = users.User(email='me@gmail.com', name='Me Meyer')
    me.save()


def _external_objects_url(url):
    """ Replaces the given internal object storage url (minio) with an URL that allows
        external access. """
    port_with_colon = ''
    if config.services.objects_port > 0:
        port_with_colon = ':%d' % config.services.objects_port

    return url.replace(
        '%s:%s' % (config.minio.host, config.minio.port),
        '%s%s%s' % (config.services.objects_host, port_with_colon, config.services.objects_base_path))


def _updated_proc(upload: users.Upload) -> Tuple[UploadProc, bool]:
    is_stale = False

    if upload.proc:
        proc = UploadProc(**upload.proc)
        if proc.update_from_backend():
            upload.proc = proc
            upload.save()

        if proc.current_task_name == proc.task_names[0] and upload.upload_time is None:
            is_stale = (datetime.now() - upload.create_time).days > 1

    else:
        proc = None

    return proc, is_stale


def _render(upload: users.Upload, proc: UploadProc, is_stale: bool) -> dict:
    data = {
        'name': upload.name,
        'upload_id': upload.upload_id,
        'presigned_url': _external_objects_url(upload.presigned_url),
        'presigned_orig': upload.presigned_url,
        'create_time': upload.create_time.isoformat() if upload.create_time is not None else None,
        'upload_time': upload.upload_time.isoformat() if upload.upload_time is not None else None,
        'proc_time': upload.proc_time.isoformat() if upload.proc_time is not None else None,
        'is_stale': is_stale,
        'proc': proc
    }

    return {key: value for key, value in data.items() if value is not None}


def _update_and_render(upload: users.Upload) -> dict:
    """
    If the given upload as a processing state attached, it will attempt to update this
    state and store the results, before the upload is rendered for the client.
    """
    proc, is_stale = _updated_proc(upload)
    return _render(upload, proc, is_stale)


class Uploads(Resource):

    def get(self):
        return [_update_and_render(user) for user in users.Upload.objects()], 200

    def post(self):
        json_data = request.get_json()
        if json_data is None:
            json_data = {}

        upload = users.Upload(user=me, name=json_data.get('name'))
        upload.save()

        upload.presigned_url = files.get_presigned_upload_url(upload.upload_id)
        upload.create_time = datetime.now()
        upload.proc = UploadProc(upload.upload_id)
        upload.save()

        return _update_and_render(upload), 200


class Upload(Resource):
    def get(self, upload_id):
        try:
            upload = users.Upload.objects(id=upload_id).first()
        except mongoengine.errors.ValidationError:
            abort(400, message='%s is not a valid upload id.' % upload_id)

        if upload is None:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        return _update_and_render(upload), 200

    def delete(self, upload_id):
        try:
            upload = users.Upload.objects(id=upload_id).first()
        except mongoengine.errors.ValidationError:
            print('###')
            abort(400, message='%s is not a valid upload id.' % upload_id)

        if upload is None:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        proc, is_stale = _updated_proc(upload)
        if not (proc.ready() or is_stale or proc.current_task_name == 'uploading'):
            abort(400, message='%s has not finished processing.' % upload_id)

        logger = get_logger(__name__, upload_id=upload_id)
        with logger.lnr_error('Delete upload file'):
            try:
                files.Upload(upload.upload_id).delete()
            except KeyError:
                logger.error('Upload exist, but file does not exist.')

        if proc.upload_hash is not None:
            with logger.lnr_error('Deleting archives.'):
                files.delete_archives(proc.upload_hash)

            with logger.lnr_error('Deleting indexed calcs.'):
                search.Calc.delete_all(upload_id=proc.upload_id)

        with logger.lnr_error('Deleting user upload.'):
            upload.delete()

        return _render(upload, proc, is_stale), 200


class RepoCalc(Resource):
    @staticmethod
    def _render(data: dict):
        upload_time = data.get('upload_time', None)
        if upload_time is not None and isinstance(upload_time, datetime):
            data['upload_time'] = data['upload_time'].isoformat()

        return {key: value for key, value in data.items() if value is not None}

    def get(self, upload_hash, calc_hash):
        try:
            data = search.Calc.get(id='%s/%s' % (upload_hash, calc_hash))
        except NotFoundError:
            abort(404, message='There is no calculation for %s/%s' % (upload_hash, calc_hash))
        except Exception as e:
            abort(500, message=str(e))

        return RepoCalc._render(data.to_dict()), 200


class RepoCalcs(Resource):
    def get(self):
        # TODO use argparse? bad request reponse an bad params, pagination as decorator
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 10))

        assert page >= 1
        assert per_page > 0

        body = {
            'from': (page - 1) * per_page,
            'size': per_page,
            'query': {
                'match_all': {}
            }
        }

        try:
            results = search.Calc.search(body=body)
        except Exception as e:
            get_logger(__name__).error('Could not execute repo calcs get.', exc_info=e)
            abort(500, message=str(e))

        return {
            'pagination': {
                'total': results['hits']['total'],
                'page': page,
                'per_page': per_page
            },
            'results': [RepoCalc._render(hit['_source']) for hit in results['hits']['hits']]
        }


@app.route('%s/archive/<string:upload_hash>/<string:calc_hash>' % base_path, methods=['GET'])
def get_calc(upload_hash, calc_hash):
    archive_id = '%s/%s' % (upload_hash, calc_hash)
    logger = get_logger(__name__, archive_id=archive_id)
    try:
        url = _external_objects_url(files.archive_url(archive_id))
        return redirect(url, 302)
    except KeyError:
        abort(404, message='Archive %s does not exist.' % archive_id)
    except Exception as e:
        logger.error('Exception on accessing archive', exc_info=e)
        abort(500, message='Could not accessing the archive.')


api.add_resource(Uploads, '%s/uploads' % base_path)
api.add_resource(Upload, '%s/uploads/<string:upload_id>' % base_path)
api.add_resource(RepoCalcs, '%s/repo' % base_path)
api.add_resource(RepoCalc, '%s/repo/<string:upload_hash>/<string:calc_hash>' % base_path)


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    app.run(debug=True, port=8000)
