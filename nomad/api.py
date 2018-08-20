from flask import Flask
from flask_restful import Resource, Api, abort
from datetime import datetime
from threading import Thread
import logging
import mongoengine.errors

from nomad import users, files, processing

logger = logging.getLogger(__name__)

app = Flask(__name__)
api = Api(app)


# provida a fake user for testing
me = users.User.objects(email='me@gmail.com').first()
if me is None:
    me = users.User(email='me@gmail.com', name='Me Meyer')
    me.save()


class Uploads(Resource):

    @staticmethod
    def _render(upload: users.Upload):
        data = {
            'id': upload.upload_id,
            'presigned_url': upload.presigned_url,
            'create_time': upload.create_time.isoformat() if upload.create_time is not None else None,
            'upload_time': upload.upload_time.isoformat() if upload.upload_time is not None else None
        }

        return {key: value for key, value in data.items() if value is not None}

    def get(self):
        return [Uploads._render(user) for user in users.Upload.objects()], 200

    def post(self):
        upload = users.Upload(user=me)
        upload.save()

        upload.presigned_url = files.get_presigned_upload_url(upload.upload_id)
        upload.create_time = datetime.now()
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


api.add_resource(Uploads, '/uploads')
api.add_resource(Upload, '/uploads/<string:upload_id>')

if __name__ == '__main__':
    @files.upload_put_handler
    def handle_upload_put(received_upload_id: str):
        upload = users.Upload.objects(id=received_upload_id).first()
        if upload is None:
            logger.error(
                'Received upload put event on non existing upload %s.' %
                received_upload_id)
            return

        upload.upload_time = datetime.now()
        try:
            proc = processing.UploadProcessing(received_upload_id)
            proc.start()
            upload.processing = proc.to_json()
        except Exception as e:
            logger.error(
                'Unexpected exception while starting processing of upload %s.' %
                received_upload_id, exc_info=e)

        upload.save()

    def handle_uploads():
        handle_upload_put(received_upload_id='provided by decorator')

    handle_uploads_thread = Thread(target=handle_uploads)
    handle_uploads_thread.start()

    app.run(debug=True)
    handle_uploads_thread.join()
