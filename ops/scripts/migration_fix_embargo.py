from nomad import infrastructure
from nomad import processing
from nomad import search

infrastructure.setup_logging()
calcs = infrastructure.setup_mongo().fairdi_nomad_migration.calc
index = infrastructure.setup_mongo().coe_migration.source_calc
uploads_col = infrastructure.setup_mongo().fairdi_nomad_migration.upload
infrastructure.setup_elastic()


def check_and_fix(upload):
    example = calcs.find_one({'upload_id': upload, 'metadata.pid': {'$exists': True}})
    if example is None:
        # can happen on multi package uploads
        return

    pid = example['metadata']['pid']
    truth = index.find_one({'_id': pid})

    if truth['metadata']['with_embargo'] != example['metadata']['with_embargo']:
        u = uploads_col.find_one({'_id': upload})
        print('need to fix from user %d, %s package id %s' % (example['metadata']['uploader']['id'], upload, u['name']))



for upload in calcs.distinct('upload_id'):
    check_and_fix(upload)
