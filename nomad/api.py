from flask import Flask
from flask_restful import Resource, Api, abort
from datetime import datetime
from threading import Thread
import mongoengine.errors

from nomad import users, files, processing
from nomad.utils import get_logger

app = Flask(__name__)
api = Api(app)


# provid a fake user for testing
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
            'upload_time': upload.upload_time.isoformat() if upload.upload_time is not None else None,
            'upload_hash': upload.upload_hash,
        }

        if upload.processing is not None:
            proc = processing.UploadProcessing.from_result_backend(upload.upload_id, upload.processing)
            data['processing'] = {
                'status': proc.status,
                'parse_specs': proc.parse_specs,
                'processing_results': proc.processing_results,
                'current_task': proc.task_name,
                'error': proc.cause.__str__()
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


def start_upload_handler(quit=False):
    """
    Starts a notification handler for uploads in a different thread. This handler
    will initiate processing for all received upload events. The processing status
    will be saved to the users db.

    Arguments:
        quit: If true, will only handling one event and stop. Otherwise run forever.
    """
    @files.upload_put_handler
    def handle_upload_put(received_upload_id: str):
        logger = get_logger(__name__, upload_id=received_upload_id)

        try:
            upload = users.Upload.objects(id=received_upload_id).first()
            if upload is None:
                logger.error('Upload does not exist')
                raise Exception()

            with logger.lnr_error('Save upload time'):
                upload.upload_time = datetime.now()
                upload.save()

            with logger.lnr_error('Start processing'):
                proc = processing.UploadProcessing(received_upload_id)
                proc.start()
                upload.processing = proc.result_tuple
                upload.save()
        except Exception:
            pass

        if quit:
            raise StopIteration

    def handle_uploads():
        handle_upload_put(received_upload_id='provided by decorator')

    handle_uploads_thread = Thread(target=handle_uploads)
    handle_uploads_thread.start()

    return handle_uploads_thread


if __name__ == '__main__':
    handle_uploads_thread = start_upload_handler()
    app.run(debug=True)
    handle_uploads_thread.join()
