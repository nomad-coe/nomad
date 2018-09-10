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
import time
import logging
from elasticsearch import NotFoundError

from nomad import config
from nomad.parsing import LocalBackend
from nomad.repo import AlreadyExists, RepoCalc, key_mappings

from tests.test_normalizing import normalized_vasp_example  # pylint: disable=unused-import
from tests.test_parsing import parsed_vasp_example  # pylint: disable=unused-import
from tests.test_files import assert_not_exists


@pytest.fixture(scope='function')
def example_elastic_calc(normalized_vasp_example: LocalBackend, caplog) \
        -> Generator[RepoCalc, None, None]:
    try:
        caplog.set_level(logging.ERROR)
        RepoCalc.get(id='test_upload_hash/test_calc_hash').delete()
    except Exception:
        pass
    finally:
        caplog.set_level(logging.WARNING)

    entry = RepoCalc.create_from_backend(
        normalized_vasp_example,
        upload_hash='test_upload_hash',
        calc_hash='test_calc_hash',
        upload_id='test_upload_id',
        mainfile='/test/mainfile',
        upload_time=datetime.now(),
        staging=True, restricted=False, user_id='me')
    time.sleep(1)  # eventually consistent?

    yield entry

    try:
        caplog.set_level(logging.ERROR)
        entry.delete()
    except Exception:
        pass
    finally:
        caplog.set_level(logging.WARNING)


def assert_elastic_calc(calc: RepoCalc):
    assert calc is not None
    for property in RepoCalc._doc_type.mapping:
        property = key_mappings.get(property, property)
        assert getattr(calc, property) is not None


def test_create_elasitc_calc(example_elastic_calc: RepoCalc):
    assert_elastic_calc(example_elastic_calc)
    assert RepoCalc.upload_exists(example_elastic_calc.upload_hash)

    get_result: RepoCalc = RepoCalc.get(
        id='%s/%s' % (example_elastic_calc.upload_hash, example_elastic_calc.calc_hash))
    assert_elastic_calc(get_result)


def test_create_existing_elastic_calc(
        example_elastic_calc: RepoCalc, normalized_vasp_example, caplog):
    try:
        caplog.set_level(logging.ERROR)
        RepoCalc.create_from_backend(
            normalized_vasp_example,
            upload_hash='test_upload_hash',
            calc_hash='test_calc_hash',
            upload_id='test_upload_id',
            mainfile='/test/mainfile',
            upload_time=datetime.now(),
            staging=True, restricted=False, user_id='me')
        assert False
    except AlreadyExists:
        caplog.set_level(logging.WARNING)
        pass
    else:
        assert False


def test_delete_elastic_calc(example_elastic_calc: RepoCalc, caplog):
    example_elastic_calc.delete()

    assert_not_exists(config.files.archive_bucket, 'test_upload_hash/test_calc_hash')
    try:
        caplog.set_level(logging.ERROR)
        RepoCalc.get(id='test_upload_hash/test_calc_hash')
        assert False
    except NotFoundError:
        pass
    else:
        assert False
    finally:
        caplog.set_level(logging.WARNING)
