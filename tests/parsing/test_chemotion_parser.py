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
import numpy as np
from nomad.config import config
from nomad.parsing.chemotion.chemotion import (
    ChemotionParser,
    _element_type_section_mapping,
)
from nomad.files import StagingUploadFiles
from tests.parsing.test_elabftw_parser import _assert_parsed_data
from tests.processing.test_data import run_processing


def _assert_chemotion(test_archive):
    assert test_archive.data is not None
    assert test_archive.data.Collection is not None

    assert test_archive.data.Collection[0].label == 'Modification Sequence'
    assert (
        test_archive.data.Collection[0].user_id
        == '60c41de1-b83d-4487-a599-8eb310847b8a'
    )
    assert test_archive.data.Collection[0].is_locked is False

    assert len(test_archive.data.Sample) == 4
    assert test_archive.data.Sample[0].xref == {
        'cas': {'label': '554-95-0', 'value': '554-95-0'}
    }
    assert test_archive.data.Sample[1].name == 'Aqua dest.'
    assert test_archive.data.Sample[2].target_amount_value == np.float16(0.001)
    assert test_archive.data.Sample[3].target_amount_value == np.float16(0.002)

    assert len(test_archive.data.CollectionsSample) == 4
    assert len(test_archive.data.Fingerprint) == 3
    assert len(test_archive.data.Molecule) == 3
    assert test_archive.data.Molecule[1].inchikey == 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'

    for k in _element_type_section_mapping.keys():
        k = 'Reactions' if k == 'Reaction' else k
        assert k in test_archive.data.m_def.all_properties
        assert test_archive.data.m_def.all_properties[k] is not None


@pytest.mark.timeout(config.tests.default_timeout)
def test_chemotion_parser(raw_files_function, proc_infra, api_v1, user1):
    upload = run_processing(
        ('test_upload', 'tests/data/parsers/chemotion/test.zip'), user1
    )

    assert upload.total_entries_count == 2
    assert len(upload.successful_entries) == 2

    with upload.entries_metadata() as entries:
        _assert_parsed_data(
            upload.upload_id,
            entries,
            StagingUploadFiles,
            ChemotionParser(),
            _assert_chemotion,
            'export.json',
            published=False,
        )
