from flask import Flask, request, redirect
from flask_restful import Resource, Api, abort
from flask_cors import CORS
from elasticsearch.exceptions import NotFoundError

from nomad import config, files
from nomad.utils import get_logger
from nomad.data import Calc, Upload, User, InvalidId, NotAllowedDuringProcessing

base_path = config.services.api_base_path

app = Flask(
    __name__,
    static_url_path='%s/docs' % base_path,
    static_folder='../docs/.build/html')
CORS(app)
api = Api(app)


# provid a fake user for testing
me = User.objects(email='me@gmail.com').first()
if me is None:
    me = User(email='me@gmail.com', name='Me Meyer')
    me.save()


class UploadsRes(Resource):

    def get(self):
        return [upload.json_dict for upload in Upload.user_uploads(me)], 200

    def post(self):
        json_data = request.get_json()
        if json_data is None:
            json_data = {}

        return Upload.create(user=me, name=json_data.get('name', None)).json_dict, 200


class UploadRes(Resource):
    def get(self, upload_id):
        try:
            return Upload.get(upload_id=upload_id).json_dict, 200
        except InvalidId:
            abort(400, message='%s is not a valid upload id.' % upload_id)
        except KeyError:
            abort(404, message='Upload with id %s does not exist.' % upload_id)

    def delete(self, upload_id):
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
