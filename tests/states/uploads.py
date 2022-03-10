#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from nomad import infrastructure
from nomad.utils.exampledata import ExampleData
from .archives.create_archives import archive_dft_bulk


def one_entry():
    '''
    State containing one entry
    '''
    infrastructure.setup()
    main_author = infrastructure.keycloak.get_user(username='test')
    data = ExampleData(main_author=main_author)

    upload_id = 'dft_upload'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)
    entry_id = 'dft_bulk'
    data.create_entry(
        upload_id=upload_id,
        entry_id=entry_id,
        mainfile='vasp.xml',
        entry_archive=archive_dft_bulk()
    )

    data.save()


def multiple_entries():
    '''
    State containing multiple DFT entries
    '''
    infrastructure.setup()
    main_author = infrastructure.keycloak.get_user(username='test')
    data = ExampleData(main_author=main_author)

    upload_id = 'dft_upload_1'
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)

    for i in range(1, 5):
        entry_id = f'dft_bulk_{i}'
        data.create_entry(
            upload_id=upload_id,
            entry_id=entry_id,
            mainfile=f'vasp_{i}.xml',
            entry_archive=archive_dft_bulk()
        )

    data.save()
