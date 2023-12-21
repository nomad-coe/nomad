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
import os

from nomad.datamodel.datamodel import EntryArchive
from nomad.datamodel.context import ClientContext
from nomad.utils.exampledata import ExampleData


import pytest

from tests.normalizing.conftest import run_processing, run_normalize
from nomad.datamodel.data import User


def test_substance(raw_files_function, test_user, mongo):
    directory = 'tests/data/datamodel/metainfo/eln'
    mainfile = 'test_substance.archive.yaml'
    upload_id = 'test_upload_id'
    data = ExampleData(main_author=test_user)
    data.create_upload(upload_id=upload_id, published=False)
    context = ClientContext(local_dir=directory, upload_id=upload_id)

    test_archive = data.create_entry_from_file(
        upload_id=upload_id,
        mainfile=os.path.join(directory, mainfile),
        entry_archive=EntryArchive(m_context=context),
    )

    data.save(with_es=False)

    # Check that entry type is the custom 'my_substance' class
    assert test_archive.metadata.entry_type == 'my_substance'
    # Check that api call to CAS found Lead Iodide
    assert test_archive.data.cas_number == '10101-63-0'
    # Check that the material results section was populated by the normalizer
    assert 'I' in test_archive.results.material.elements
    assert 'Pb' in test_archive.results.material.elements


def test_ensemble(raw_files_function, test_user, mongo):
    directory = 'tests/data/datamodel/metainfo/eln'
    mainfile = 'test_ensemble.archive.yaml'
    upload_id = 'test_upload_id'
    data = ExampleData(main_author=test_user)
    data.create_upload(upload_id=upload_id, published=False)
    context = ClientContext(local_dir=directory, upload_id=upload_id)

    test_archive = data.create_entry_from_file(
        upload_id=upload_id,
        mainfile=os.path.join(directory, mainfile),
        entry_archive=EntryArchive(m_context=context),
    )

    data.save(with_es=False)

    test_archive.metadata.main_author = User.m_from_dict(
        {
            'first_name': 'Sheldon',
            'last_name': 'Cooper',
            'affiliation': 'Humboldt-Universität zu Berlin',
        }
    )
    run_normalize(test_archive)

    # Check that entry type is the custom 'my_ensemble' class
    assert test_archive.metadata.entry_type == 'my_ensemble'
    # Check that sample id was generated correctly from the author metadata
    assert test_archive.data.lab_id == 'HUB_ShCo_19930101_My-ensemble'
