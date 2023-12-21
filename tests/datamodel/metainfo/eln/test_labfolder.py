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
import pytest
from nomad import utils
import json
import requests

from tests.normalizing.conftest import run_processing
from nomad.processing.data import Upload
from nomad.datamodel.context import ServerContext
from nomad.datamodel.metainfo.eln.labfolder import (
    LabfolderProject,
    LabfolderImportError,
)
from nomad.utils.exampledata import ExampleData
from tests.test_files import example_mainfile_contents, append_raw_files  # pylint: disable=unused-import


def test_labfolder_integration(mongo, monkeypatch, test_user):
    directory = 'tests/data/datamodel/metainfo/eln/material_library'
    mainfile = 'example-labfolder.archive.json'

    test_archive = run_processing(directory, mainfile)

    # assert archive
    assert test_archive.metadata.entry_type == 'LabfolderProject'
    assert test_archive.data.project_url is not None
    assert len(test_archive.data.entries) == 6


@pytest.mark.parametrize(
    'status_code,project_url,labfolder_email,password,element_data,response_data',
    [
        pytest.param(
            200,
            'https://labfolder.labforward.app/eln/notebook#?projectIds=1',
            'test_email',
            'test_password',
            [{'elements': [{'id': '1', 'version_id': '1', 'type': 'TEXT'}]}],
            {
                'entry_id': '1',
                'id': '1',
                'version_id': '1',
                'creation_date': '2022-08-29T09:22:30+00:00',
                'version_date': '2022-08-29T09:22:30+00:00',
                'owner_id': '1',
                'element_type': 'TEXT',
                'content': '<p>test content</p>',
            },
            id='TEXT_element',
        ),
        pytest.param(
            200,
            'https://labfolder.labforward.app/eln/notebook#?projectIds=1',
            'test_email',
            'test_password',
            [{'elements': [{'id': '1', 'version_id': '1', 'type': 'TABLE'}]}],
            {
                'entry_id': '1',
                'id': '1',
                'version_id': '1',
                'creation_date': '2022-08-29T09:22:33+00:00',
                'version_date': '2022-08-29T09:22:33+00:00',
                'owner_id': '1',
                'element_type': 'TABLE',
                'title': 'Table Title',
                'content': {'version': '14.1.3', 'nested_content': 'more content'},
            },
            id='TABLE_element',
        ),
        pytest.param(
            200,
            'https://labfolder.labforward.app/eln/notebook#?projectIds=1',
            'test_email',
            'test_password',
            [{'elements': [{'id': '1', 'version_id': '1', 'type': 'WELL_PLATE'}]}],
            {
                'entry_id': '1',
                'id': '1',
                'version_id': '1',
                'creation_date': '2022-08-29T09:22:33+00:00',
                'version_date': '2022-08-29T09:22:33+00:00',
                'owner_id': '1',
                'element_type': 'WELL_PLATE',
                'title': 'Wellplate title',
                'content': {'version': '14.1.3', 'nested_content': 'more _content'},
                'meta_data': {
                    'plate': {'size': '24'},
                    'activeLayer': {'name': 'Composite', 'type': 'COMPOSITE'},
                },
            },
            id='WELL_PLATE_element',
        ),
        pytest.param(
            200,
            'https://labfolder.labforward.app/eln/notebook#?projectIds=1',
            'test_email',
            'test_password',
            [{'elements': [{'id': '1', 'version_id': '1', 'type': 'DATA'}]}],
            {
                'entry_id': '1',
                'id': '1',
                'version_id': '1',
                'creation_date': '2022-08-29T09:22:33+00:00',
                'version_date': '2022-08-29T09:22:34+00:00',
                'owner_id': '1',
                'element_type': 'DATA',
                'data_elements': [
                    {
                        'type': 'DATA_ELEMENT_GROUP',
                        'title': 'first title',
                        'children': [
                            {
                                'type': 'SINGLE_DATA_ELEMENT',
                                'title': 'child title',
                                'value': '1',
                                'unit': 'mL',
                                'physical_quantity_id': '6',
                            },
                            {
                                'type': 'DESCRIPTIVE_DATA_ELEMENT',
                                'title': 'Water',
                                'description': 'add up to 1 liter',
                            },
                        ],
                    }
                ],
            },
            id='DATA_element',
        ),
        pytest.param(
            400,
            'https://labfolder.labforward.app/eln/notebook#?projectIds=1',
            None,
            'test_password',
            None,
            None,
            id='missing_email',
        ),
        pytest.param(
            400,
            'https://labfolder.labforward.app/eln/notebook#?projectIds=1',
            'test_email',
            None,
            None,
            None,
            id='missing_password',
        ),
        pytest.param(
            400,
            None,
            'test_email',
            'test_password',
            None,
            None,
            id='missing_project_url',
        ),
        pytest.param(
            400, 'wrong_url', 'test_email', 'test_password', None, None, id='wrong_url'
        ),
        pytest.param(
            400,
            'https://labfolder.labforward.app/eln/notebook#?projectIds=1',
            'test_email',
            'test_password',
            None,
            None,
            id='wrong_element_data',
        ),
    ],
)
def test_labfolder_detailed(
    mongo,
    monkeypatch,
    test_user,
    status_code,
    project_url,
    labfolder_email,
    password,
    element_data,
    response_data,
):
    logger = utils.get_logger(__name__)
    base_api_url = 'https://labfolder.labforward.app/api/v2'

    def mock_labfolder_json_method(*args, **kwargs):
        mock_response = MagicMock()
        mock_response.status_code = status_code
        mock_response.url = '/test_endpoint'

        if args[0] == f'{base_api_url}/entries?project_ids=1':
            mock_response.return_value = element_data
        else:
            mock_response.return_value = response_data

        mock_response.json.return_value = mock_response.return_value
        return mock_response

    def mock_labfolder_response_method(*args, **kwargs):
        mock_response = MagicMock()
        mock_response.status_code = status_code
        mock_response.json.return_value = {'token': 'test'}
        return mock_response

    # patching API calls to labfolder server
    monkeypatch.setattr(requests, 'post', mock_labfolder_response_method)
    monkeypatch.setattr(requests, 'get', mock_labfolder_json_method)

    # processing preparation
    data = ExampleData(main_author=test_user)
    data.create_upload(upload_id='test_upload_id', published=False)

    data.create_entry(
        upload_id='test_upload_id',
        entry_id='test_upload_id',
        mainfile='tests/data/datamodel/metainfo/eln/material_library/example-labfolder.archive.json',
    )

    data.save(with_es=False)

    upload = Upload.objects(upload_id='test_upload_id').first()
    assert upload is not None

    context = ServerContext(upload=upload)
    test_archive = data.archives['test_upload_id']
    test_archive.m_context = context

    # Creating labfolder instance
    labfolder_instance = LabfolderProject()
    labfolder_instance.project_url = project_url
    labfolder_instance.labfolder_email = labfolder_email
    labfolder_instance.password = password
    labfolder_instance.resync_labfolder_repository = True
    test_archive.data = labfolder_instance

    if status_code is 200:
        labfolder_instance.normalize(test_archive, logger=logger)
        assert len(labfolder_instance.entries) is 1
        assert labfolder_instance.labfolder_email is None
        assert labfolder_instance.password is None

        with test_archive.m_context.raw_file(test_archive.metadata.mainfile, 'r') as f:
            archive_data = json.load(f)['data']
            assert archive_data.get('password', None) is None
            assert archive_data.get('labfolder_email', None) is None

        parsed_data = labfolder_instance.entries[0].elements[0].m_to_dict()
        del parsed_data['m_def']
        if parsed_data['element_type'] is not 'DATA':
            assert json.dumps(parsed_data, sort_keys=True) == json.dumps(
                response_data, sort_keys=True
            )
        else:
            del parsed_data['labfolder_data']
            assert json.dumps(parsed_data, sort_keys=True) == json.dumps(
                response_data, sort_keys=True
            )
    else:
        with pytest.raises(LabfolderImportError):
            labfolder_instance.normalize(test_archive, logger=logger)
