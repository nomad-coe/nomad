from nomad.processing.data import Upload
import time


if __name__ == '__main__':
    suspicious = {}

    while True:
        for upload in Upload.objects(status='RUNNING', current_task='parse_all'):
            if upload.total_calcs == upload.processed_calcs:
                if upload.upload_id in suspicious:
                    del(suspicious[upload.upload_id])
                    upload.status = 'SUCCESS'
                    upload.save()
                    print('Fixed suspicious %s' % upload.upload_id)
                else:
                    print('Found suspicious %s' % upload.upload_id)
                    suspicious[upload.upload_id] = upload.upload_id
        time.sleep(1)

