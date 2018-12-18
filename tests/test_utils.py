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

import time
import json

from nomad import utils


def test_timer(caplog):
    with utils.timer(utils.get_logger('test_logger'), 'test measure'):
        time.sleep(0.1)

    assert json.loads(caplog.record_tuples[0][2])['event'] == 'test measure'


def test_sanitize_logevent():
    assert utils.sanitize_logevent('numbers 2 and 45.2') == 'numbers X and X'
    assert utils.sanitize_logevent('list [2, 3.3, 10] and (273.9, .92)') == 'list L and L'
    assert utils.sanitize_logevent('mat [2, [3.3, 2], 10]') == 'mat M'


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
