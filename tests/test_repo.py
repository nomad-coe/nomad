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

import pytest
from typing import Generator
from datetime import datetime
from elasticsearch import NotFoundError

from nomad.files import ArchiveFile, UploadFile
from nomad.parsing import LocalBackend
from nomad.repo import AlreadyExists, RepoCalc, key_mappings

from tests.test_files import example_file  # noqa
from tests.test_normalizing import normalized_template_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_template_example  # pylint: disable=unused-import


@pytest.fixture(scope='function')
def example_elastic_calc(normalized_template_example: LocalBackend, elastic) \
        -> Generator[RepoCalc, None, None]:

    upload_file = UploadFile('test_upload_id', local_path=example_file)
    mainfile = next(filename for filename in upload_file.filelist if 'template.json' in filename)
    auxfiles = list(upload_file.get_siblings(mainfile))

    try:
        calc = RepoCalc.get(id='test_upload_hash/test_calc_hash')
    except NotFoundError:
        pass
    else:
        calc.delete()

    entry = RepoCalc.create_from_backend(
        normalized_template_example,
        upload_hash='test_upload_hash',
        calc_hash='test_calc_hash',
        upload_id='test_upload_id',
        additional=dict(
            mainfile=mainfile,
            upload_time=datetime.now(),
            staging=True, restricted=False, user_id='me@gmail.com',
            aux_files=auxfiles),
        refresh='true')

    yield entry

    try:
        calc = RepoCalc.get(id='test_upload_hash/test_calc_hash')
    except NotFoundError:
        pass
    else:
        calc.delete()


def assert_elastic_calc(calc: RepoCalc):
    assert calc is not None
    for property in RepoCalc._doc_type.mapping:
        property = key_mappings.get(property, property)
        assert getattr(calc, property) is not None

    assert len(getattr(calc, 'aux_files')) > 0


def test_create_elastic_calc(example_elastic_calc: RepoCalc, no_warn):
    assert_elastic_calc(example_elastic_calc)
    assert RepoCalc.upload_exists(example_elastic_calc.upload_hash)

    get_result: RepoCalc = RepoCalc.get(
        id='%s/%s' % (example_elastic_calc.upload_hash, example_elastic_calc.calc_hash))
    assert_elastic_calc(get_result)


def test_create_existing_elastic_calc(
        example_elastic_calc: RepoCalc, normalized_template_example):
    try:
        RepoCalc.create_from_backend(
            normalized_template_example,
            upload_hash='test_upload_hash',
            calc_hash='test_calc_hash',
            upload_id='test_upload_id',
            additional=dict(
                mainfile='/test/mainfile',
                upload_time=datetime.now(),
                staging=True, restricted=False, user_id='me'),
            refresh='true')
        assert False
    except AlreadyExists:
        pass
    else:
        assert False


def test_delete_elastic_calc(example_elastic_calc: RepoCalc):
    example_elastic_calc.delete()

    assert not ArchiveFile('test_upload_hash/test_calc_hash').exists()
    try:
        RepoCalc.get(id='test_upload_hash/test_calc_hash')
        assert False
    except NotFoundError:
        pass
    else:
        assert False


def test_staging_elastic_calc(example_elastic_calc: RepoCalc, no_warn):
    assert RepoCalc.get(id='test_upload_hash/test_calc_hash').staging


def test_unstage_elastic_calc(example_elastic_calc: RepoCalc, no_warn):
    RepoCalc.unstage(upload_id='test_upload_id', staging=False)

    assert not RepoCalc.get(id='test_upload_hash/test_calc_hash').staging
