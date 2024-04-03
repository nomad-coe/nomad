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
from unittest.mock import MagicMock
from urllib.parse import urlparse

import pytest
import json
import pybis

from nomad.datamodel import ServerContext
from nomad.datamodel.metainfo.eln.openbis import OpenbisEntry, OpenbisImportError
from nomad.processing import Upload
from nomad.utils.exampledata import ExampleData
from nomad import utils


def mocked_login(url):
    if url is None:
        if url is None:
            raise ValueError('please provide a URL you want to connect to.')

    if not url.startswith('http'):
        url = 'https://' + url

    url_obj = urlparse(url)
    if url_obj.netloc is None or url_obj.netloc == '':
        raise ValueError(
            'please provide the url in this format: https://openbis.host.ch:8443'
        )
    if url_obj.hostname is None:
        raise ValueError('hostname is missing')
    if url_obj.scheme == 'http':
        raise ValueError('always use https!')


@pytest.mark.parametrize(
    'status_code,project_url,username,password,space_data,project_data,experiment_data',
    [
        pytest.param(
            200,
            'https://openbis-server.ethz.ch/openbis/webapp/eln-lims/',
            'test_email',
            'test_password',
            {
                'code': 'code_value',
                'frozen': False,
                'permId': 'permId_value',
            },
            {
                'code': 'code_value',
                'frozen': False,
                'permId': 'permId_value',
                'identifier': '/DEFAULT/DEFAULT',
            },
            {'code': 'code_value', 'frozen': False, 'permId': 'permId_value'},
            id='successful parsing',
        ),
        pytest.param(
            400,
            'http://openbis-server.ethz.ch/openbis/webapp/eln-lims/',
            'test_email',
            'test_password',
            {
                'code': 'code_value',
                'frozen': False,
                'permId': 'permId_value',
            },
            {
                'code': 'code_value',
                'frozen': False,
                'permId': 'permId_value',
                'identifier': '/DEFAULT/DEFAULT',
            },
            {'code': 'code_value', 'frozen': False, 'permId': 'permId_value'},
            id='incorrect project_url',
        ),
    ],
)
def test_openbis(
    mongo_function,
    monkeypatch,
    user1,
    status_code,
    project_url,
    username,
    password,
    space_data,
    project_data,
    experiment_data,
):
    logger = utils.get_logger(__name__)

    class AttrsAll:
        def __init__(self, test_data):
            self.data = test_data

        def all(self):
            return self.data

    class AttrsDescriptor:
        def __get__(self, instance, owner):
            data = instance.data if instance else {}
            return AttrsAll(data)

    class MockedSpaces:
        def __init__(self):
            self.data = space_data

        attrs = AttrsDescriptor()

        @staticmethod
        def get_projects():
            class MockedProjects:
                def __init__(self):
                    self.data = project_data

                attrs = AttrsDescriptor()

                @staticmethod
                def get_experiments():
                    class MockedExperiments:
                        def __init__(self):
                            self.data = experiment_data

                        attrs = AttrsDescriptor()

                    return [MockedExperiments()]

            return [MockedProjects()]

    def mock_openbis_response(*args, **kwargs):
        mock_response = MagicMock()
        mock_response.login.return_value = mocked_login(project_url)
        mock_response.logout.return_value = None
        mock_response.get_spaces.return_value = [MockedSpaces()]
        return mock_response

    # patching openbis attrs
    monkeypatch.setattr(pybis, 'Openbis', mock_openbis_response)

    data = ExampleData(main_author=user1)
    data.create_upload(upload_id='test_upload_id', published=False)

    data.create_entry(
        upload_id='test_upload_id',
        entry_id='test_upload_id',
        mainfile='tests/data/datamodel/metainfo/eln/material_library/example-openbis.archive.json',
    )

    data.save(with_es=False)

    upload = Upload.objects(upload_id='test_upload_id').first()
    assert upload is not None

    context = ServerContext(upload=upload)
    test_archive = data.archives['test_upload_id']
    test_archive.m_context = context
    openbis_instance = OpenbisEntry()
    openbis_instance.project_url = project_url
    openbis_instance.username = username
    openbis_instance.password = password
    test_archive.data = openbis_instance

    if status_code is 200:
        openbis_instance.normalize(test_archive, logger=logger)
        assert len(openbis_instance.spaces) is 1
        assert openbis_instance.username is None
        assert openbis_instance.password is None

        parsed_openbis_spaces = openbis_instance.spaces[0].m_to_dict()
        parsed_openbis_projects = parsed_openbis_spaces.pop('projects')[0]
        parsed_openbis_experiments = parsed_openbis_projects.pop('experiments')[0]
        assert json.dumps(parsed_openbis_spaces, sort_keys=True) == json.dumps(
            space_data, sort_keys=True
        )
        assert json.dumps(parsed_openbis_projects, sort_keys=True) == json.dumps(
            project_data, sort_keys=True
        )
        assert json.dumps(parsed_openbis_experiments, sort_keys=True) == json.dumps(
            experiment_data, sort_keys=True
        )
    if status_code is 400:
        with pytest.raises(OpenbisImportError):
            openbis_instance.normalize(test_archive, logger=logger)
