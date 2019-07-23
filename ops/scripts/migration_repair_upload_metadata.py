from nomad import infrastructure
from dateutil.parser import parse

infrastructure.setup_logging()
mongo = infrastructure.setup_mongo()
calcs = mongo.fairdi_nomad_migration.calc
uploads = mongo.fairdi_nomad_migration.upload


def doit(upload):
    example = calcs.find_one({'upload_id': upload})
    user_id = example['metadata']['uploader']['id']
    upload_time = example['metadata']['upload_time']

    uploads.update_one(
        {'_id': upload}, 
        {'$set': {
            'user_id': str(user_id),
            'upload_time': parse(str(upload_time))
        }})


for upload in uploads.distinct('_id'):
    doit(upload)
