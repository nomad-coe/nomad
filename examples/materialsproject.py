"""
This is a brief example demonstrating the public nomad@FAIRDI API for doing operations
that are necessary during a *materials-project<->nomad* integration.
"""

# create the bravado client
nomad_url = 'http://enc-staging-nomad.esc.rzg.mpg/fairdi/nomad/mp-test/api'
user = 'phuck@lbl.gov'
password = 'password'

host = urlparse(nomad_url).netloc.split(':')[0]
http_client = RequestsClient()
http_client.set_basic_auth(host, user, password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)

# upload data
test_file = '/...'
upload = client.uploads.upload(file=test_file).response().result
while upload.tasks_running:
    upload client.uploads.get_upload(upload_id=upload.upload_id).response().result
    print('processed: %d, failures: %d' % (upload.processed_calcs, upload_failed_calcs))

if upload.tasks_status != 'SUCCESS':
    print('something went wrong')
    print('errors: %s' % str(upload.errors))
    client.uploads.delete_upload(upload_id=upload.upload_id)
    sys.exit(1)

# publish data
client.uploads.exec_upload_operation(upload_id=upload.upload_id, payload={
    'operation': 'publish',
    'metadata': {
        'comment': 'Data from materials project',
        'references': ['http://materials-project.gov'],
        # 'coauthors': ['person@lbl.gov', '...'],  this does not yet work with emails
        # 'external_id': 'a/mp/id'  this does also not work, but we could implement something like this
    }
})

# search for data
result = client.repo.get_calcs(paths='tags1/tag2').response().result
if result.pagination.total == 0:
    print('not found')
    sys.exit(1)
elif result.pagination.total > 1:
    print('my ids are not specific enough, bummer ...')
    sys.exit(1)
else:
    calc = result.results[0]

# download data
#   via api
client.raw.get(upload_id=calc.upload_id, path=calc.mainfile).response()
#   via download
#   just the 'mainfile'
url = '%s/raw/%s/%s' % (nomad_url, calc.upload_id, calc.mainfile)
#   all files
url = '%s/raw/%s/%s/*' % (nomad_url, calc.upload_id, os.path.dirname(calc.mainfile))
