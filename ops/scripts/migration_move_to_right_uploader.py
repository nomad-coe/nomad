from nomad import infrastructure

infrastructure.setup_logging()
mongo = infrastructure.setup_mongo()
calcs = mongo.fairdi_nomad_migration.calc
uploads = mongo.fairdi_nomad_migration.upload


def doit(upload):
    example = calcs.find_one({'upload_id': upload})
    user_id = example['metadata']['uploader']['id']

    uploads.update_one({'_id': upload}, {'$set': {'user_id': user_id}})


for upload in uploads.distinct('_id'):
    doit(upload)
