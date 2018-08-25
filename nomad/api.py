from flask import Flask, Response, request, redirect
from flask_restful import Resource, Api, abort
from datetime import datetime
import mongoengine.errors
from flask_cors import CORS
import logging
from elasticsearch.exceptions import NotFoundError

from nomad import users, files, search
from nomad.processing import UploadProc
from nomad.utils import get_logger

app = Flask(__name__)
CORS(app)
api = Api(app)


# provid a fake user for testing
me = users.User.objects(email='me@gmail.com').first()
if me is None:
    me = users.User(email='me@gmail.com', name='Me Meyer')
    me.save()


class Uploads(Resource):

    @staticmethod
    def _render(upload: users.Upload):
        if upload.proc:
            proc = UploadProc(**upload.proc)
            proc.update_from_backend()
        else:
            proc = None

        data = {
            'name': upload.name,
            'upload_id': upload.upload_id,
            'presigned_url': upload.presigned_url,
            'create_time': upload.create_time.isoformat() if upload.create_time is not None else None,
            'upload_time': upload.upload_time.isoformat() if upload.upload_time is not None else None,
            'proc_time': upload.proc_time.isoformat() if upload.proc_time is not None else None,
            'proc': proc
        }

        return {key: value for key, value in data.items() if value is not None}

    def get(self):
        return [Uploads._render(user) for user in users.Upload.objects()], 200

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

        return Uploads._render(upload), 200


class Upload(Resource):
    def get(self, upload_id):
        try:
            upload = users.Upload.objects(id=upload_id).first()
        except mongoengine.errors.ValidationError:
            abort(400, message='%s is not a valid upload id.' % upload_id)

        if upload is None:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

        return Uploads._render(upload), 200


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


@app.route('/archive/<string:upload_hash>/<string:calc_hash>', methods=['GET'])
def get_calc(upload_hash, calc_hash):
    archive_id = '%s/%s' % (upload_hash, calc_hash)
    logger = get_logger(__name__, archive_id=archive_id)
    try:
        url = files.archive_url(archive_id)
        return redirect(url, 302)
    except KeyError:
        abort(404, message='Archive %s does not exist.' % archive_id)
    except Exception as e:
        logger.error('Exception on accessing archive', exc_info=e)
        abort(500, message='Could not accessing the archive.')


api.add_resource(Uploads, '/uploads')
api.add_resource(Upload, '/uploads/<string:upload_id>')
api.add_resource(RepoCalcs, '/repo')
api.add_resource(RepoCalc, '/repo/<string:upload_hash>/<string:calc_hash>')


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    app.run(debug=True)
