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
    pid = example['metadata']['pid']
    truth = index.find_one({'_id': pid})

    if truth['metadata']['with_embargo'] != example['metadata']['with_embargo']:
        print('need to fix %s' % upload)

        calcs.update_many(
            {'upload_id': upload, 'metadata.with_embargo': True},
            {'$set': {'changed': True, 'metadata.with_embargo': False}})

        calcs.update_many(
            {'upload_id': upload, 'metadata.with_embargo': False, 'changed': True},
            {'$unset': {'changed': 1}, '$set': {'metadata.with_embargo': True}})

        upload_proc = processing.Upload.get_by_id(upload)
        upload_with_metadata = upload_proc.to_upload_with_metadata(upload_proc.metadata)
        calcs_with_metadata = upload_with_metadata.calcs
        search.publish(calcs_with_metadata)


for upload in calcs.distinct('upload_id'):
    check_and_fix(upload)
