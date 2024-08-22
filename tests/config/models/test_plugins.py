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

import tempfile
import pytest

from nomad.config.models.plugins import (
    ExampleUploadEntryPoint,
    example_upload_path_prefix,
)


def mock_plugin_package(monkeypatch, directory):
    """Used for mocking the presence of a plugin package installation
    location."""

    def mock_get_package_path(package_name):
        return directory

    monkeypatch.setattr(
        'nomad.config.models.plugins.get_package_path', mock_get_package_path
    )


@pytest.mark.parametrize(
    'config, expected_local_path',
    [
        pytest.param(
            {
                'title': 'test',
                'description': 'test',
                'category': 'test',
                'path': 'example_uploads/getting_started',
            },
            f'{example_upload_path_prefix}/nomad_test_plugin/example_uploads/getting_started',
            id='load with path',
        ),
        pytest.param(
            {
                'title': 'test',
                'description': 'test',
                'category': 'test',
                'url': 'https://nomad-lab.eu/prod/v1/docs/assets/nomad-oasis.zip',
            },
            f'{example_upload_path_prefix}/nomad_test_plugin/example_uploads/nomad-oasis.zip',
            id='load with url',
        ),
    ],
)
def test_example_upload_entry_point_valid(config, expected_local_path, monkeypatch):
    # Create tmp directory that will be used as a mocked package location.
    with tempfile.TemporaryDirectory() as tmp_dir_path:
        mock_plugin_package(monkeypatch, tmp_dir_path)
        config['plugin_package'] = 'nomad_test_plugin'
        entry_point = ExampleUploadEntryPoint(**config)
        entry_point.load()
        assert entry_point.local_path == expected_local_path


@pytest.mark.parametrize(
    'config, error',
    [
        pytest.param(
            {
                'title': 'test',
                'description': 'test',
                'category': 'test',
            },
            'Provide one of "path", "url" or "local_path".',
            id='no path, url or local_path given',
        ),
        pytest.param(
            {
                'title': 'test',
                'description': 'test',
                'category': 'test',
                'path': 'example_uploads/getting_started',
                'url': 'https://test.zip',
            },
            'Provide only "path" or "url", not both.',
            id='path and url both given',
        ),
        pytest.param(
            {
                'title': 'test',
                'description': 'test',
                'category': 'test',
                'url': 'https://test.zip',
            },
            'Could not fetch the example upload from URL: https://test.zip',
            id='cannot find url',
        ),
    ],
)
def test_example_upload_entry_point_invalid(config, error, monkeypatch):
    # Create tmp directory that will be used as a mocked package location.
    with tempfile.TemporaryDirectory() as tmp_dir_path:
        mock_plugin_package(monkeypatch, tmp_dir_path)
        config['plugin_package'] = 'nomad_test_plugin'
        with pytest.raises(Exception) as exc_info:
            entry_point = ExampleUploadEntryPoint(**config)
            entry_point.load()

        assert exc_info.match(error)
