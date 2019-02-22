# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
A command that runs some example operations on a working nomad@FAIRDI installation
as a final integration test.
"""

import time

from .main import cli


example_file = 'tests/data/proc/examples_vasp.zip'


@cli.command(help='Runs a few example operations as a test.')
def integrationtests():
    from .main import create_client
    client = create_client()

    print('upload with multiple code data')
    with open(example_file, 'rb') as f:
        upload = client.uploads.upload(name='integration test upload', file=f).response().result

    def get_upload():
        return client.uploads.get_upload(upload_id=upload.upload_id, per_page=100).response().result

    print('observe the upload process to be finished')
    upload = get_upload()
    while upload.tasks_running:
        time.sleep(0.3)
        upload = get_upload()

    assert upload.tasks_status == 'SUCCESS'
    total = upload.calcs.pagination.total
    assert 100 > total > 0
    assert len(upload.calcs.results) == total

    try:
        print('get repo data')
        for calc in upload.calcs.results:
            repo = client.repo.get_repo_calc(upload_id=upload.upload_id, calc_id=calc.calc_id).response().result
            repo['calc_id'] == calc.calc_id

        print('get archive data')
        for calc in upload.calcs.results:
            client.archive.get_archive_calc(upload_id=upload.upload_id, calc_id=calc.calc_id).response()

        print('get archive logs')
        for calc in upload.calcs.results:
            client.archive.get_archive_logs(upload_id=upload.upload_id, calc_id=calc.calc_id).response()

        print('perform search on data')
        search = client.repo.search(owner='staging', per_page=100).response().result
        assert search.pagination.total >= total
        assert len(search.results) <= search.pagination.total
    finally:
        print('delete the upload again')
        client.uploads.delete_upload(upload_id=upload.upload_id).response()
        while upload.process_running:
            upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result

    # TODO publish upload
    # TODO admin delete published upload -- this functionality does not yet exist
