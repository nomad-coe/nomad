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

import pytest

from nomad.app.v1.routers.metainfo import store_package_definition
from nomad.metainfo import MSection


@pytest.mark.parametrize('metainfo_data', [
    pytest.param({
        'm_def': 'nomad.metainfo.metainfo.Package',
        'name': 'test.Package',
        'section_definitions': [
            {
                'name': 'MySection'
            }
        ]
    }, id='python')
])
def test_metainfo_section_id_endpoint(metainfo_data, mongo_infra, client):
    assert MSection.from_dict(metainfo_data).m_to_dict(with_root_def=True, with_out_meta=True) == metainfo_data

    package = MSection.from_dict(metainfo_data)

    store_package_definition(package, with_root_def=True, with_out_meta=True)

    section_id = package.section_definitions[0].definition_id

    response = client.get(f'metainfo/{section_id}')
    assert response.status_code == 200
    assert response.json()['data'] == metainfo_data

    response = client.get(f'metainfo/{section_id[::-1]}')
    assert response.status_code == 404
