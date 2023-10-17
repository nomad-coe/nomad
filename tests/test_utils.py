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

import time
import json
import pytest

from nomad import utils
from nomad.utils import structlogging
from nomad.metainfo.metainfo import MSection, Quantity, SubSection
from nomad import files
from nomad.processing import Upload


def test_decode_handle_id():
    assert utils.decode_handle_id('a1') == 321
    assert utils.decode_handle_id('6i370') == 6884576
    with pytest.raises(ValueError):
        utils.decode_handle_id('zz')


def test_timer(caplog):
    with utils.timer(utils.get_logger('test_logger'), 'test measure'):
        time.sleep(0.1)

    assert json.loads(caplog.record_tuples[0][2])['event'] == 'test measure'


def test_sleep_timer():
    sleep = utils.SleepTimeBackoff(start_time=0.1)
    start = time.time()

    for _ in range(0, 3):
        sleep()

    duration = time.time() - start
    assert duration > 0.7
    assert duration < 1


def test_sanitize_logevent():
    assert structlogging.sanitize_logevent('numbers 2 and 45.2') == 'numbers X and X'
    assert structlogging.sanitize_logevent('list [2, 3.3, 10] and (273.9, .92)') == 'list L and L'
    assert structlogging.sanitize_logevent('mat [2, [3.3, 2], 10]') == 'mat M'


def test_logging(no_warn):
    utils.get_logger(__name__).info('test msg')

    received_test_event = False
    for record in no_warn.get_records(when='call'):
        assert record.levelname == 'INFO'
        data = json.loads(record.msg)
        assert 'event' in data
        assert data['event'] == 'test msg'
        received_test_event = True
    assert received_test_event


def test_common_prefix():
    assert utils.common_prefix(['aa/bb/cc', 'aa/bb/dd']) == 'aa/bb/'
    assert utils.common_prefix(['aa/bb/dc', 'aa/bb/d']) == 'aa/bb/'
    assert utils.common_prefix(['aa/b/dc', 'aa/bb/d']) == 'aa/'
    assert utils.common_prefix(['a', 'a']) == ''
    assert utils.common_prefix(['a', 'ab']) == ''
    assert utils.common_prefix(['/a', '/a']) == '/'


def test_uuid():
    uuid = utils.create_uuid()
    assert uuid is not None


def test_class_logger():
    logger = utils.ClassicLogger(__name__, test='value')
    logger.warn('hello world', test='other value')


class TestSubSection(MSection):
    name = Quantity(type=str)
    subsection_number = Quantity(type=int)


class TestSection(MSection):
    quantity_section = Quantity(type=int)
    subsection = SubSection(sub_section=TestSubSection, repeats=True)


def test_extract_section():
    test_section = TestSection(
        quantity_section=5,
        subsection=[TestSubSection(name='subsection1', subsection_number=1), TestSubSection(name='subsection2')]
    )
    # Check if quantities are properly extracted
    assert utils.extract_section(test_section, ['quantity_section']) == 5
    # Check if last section is properly extracted:
    path = ['subsection']
    assert utils.extract_section(test_section, path).name == 'subsection2'
    # Check if last section in path is not a repeated section
    assert not isinstance(utils.extract_section(test_section, path), list)
    # Check if full_list works
    assert isinstance(utils.extract_section(test_section, path, full_list=True), list)
    assert len(utils.extract_section(test_section, path, full_list=True)) == 2
    # Check if not section is reported if the last of the intermediate elements do not contain it
    path.append('subsection_number')
    assert not utils.extract_section(test_section, path)


def create_upload(upload_id, user_id, file_paths=None):
    if file_paths is None:
        file_paths = []
    upload = Upload(
        upload_id=upload_id,
        main_author=user_id)
    upload.save()
    files.StagingUploadFiles(upload_id=upload.upload_id, create=True)
    for file_path in file_paths:
        upload.staging_upload_files.add_rawfiles(file_path)
    upload.process_upload()
    upload.block_until_complete()
    return upload
