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
import json
import random

from nomad.processing.base import Proc, process, ProcessStatus

random.seed(0)


def assert_proc(proc, current_process, process_status=ProcessStatus.SUCCESS, errors=0, warnings=0):
    assert proc.current_process == current_process
    assert proc.process_status == process_status
    assert len(proc.errors) == errors
    assert len(proc.warnings) == warnings
    assert not proc.process_running


class SimpleProc(Proc):
    @process
    def a_process(self):
        pass

    @process
    def a_process_with_arguments(self, a_string, a_int, a_float, a_bool, a_list, a_tuple, a_dict, **kwargs):
        assert a_string == 'a string'
        assert a_int == 7
        assert a_float == 3.14
        assert a_bool is True
        assert a_list == [1, 2]
        assert a_tuple == [3, 4]  # Apparently, tuples become unmarshalled as lists...
        assert a_dict == {'1': 'one', '2': 2}  # ... and dictionary keys are converted to strings.
        assert kwargs == {'kwarg': 'kwarg_value'}


@pytest.mark.parametrize('with_args', [
    pytest.param(False, id='no-args'),
    pytest.param(True, id='with-args')])
def test_simple_process(worker, mongo, no_warn, with_args):
    p = SimpleProc.create()
    if with_args:
        process = 'a_process_with_arguments'
        p.a_process_with_arguments(
            'a string', 7, 3.14, True, [1, 2], (3, 4), {1: 'one', 2: 2}, kwarg='kwarg_value')
    else:
        process = 'a_process'
        p.a_process()
    p.block_until_complete()
    assert_proc(p, process)


class FailingProc(Proc):
    @process
    def will_fail(self):
        self.fail('fail fail fail!')

    @process
    def will_fail_with_exception(self):
        _x = 1 / 0


@pytest.mark.parametrize('with_exception', [
    pytest.param(False, id='no-exception'),
    pytest.param(True, id='with-exception')])
def test_failing_process(worker, mongo, with_error, with_exception):
    p = FailingProc.create()
    if with_exception:
        event = 'process failed with exception'
        process = 'will_fail_with_exception'
        p.will_fail_with_exception()
    else:
        event = 'process failed'
        process = 'will_fail'
        p.will_fail()
    p.block_until_complete()
    assert_proc(p, process, ProcessStatus.FAILURE, errors=1)

    has_log = False
    for record in with_error.get_records(when='call'):
        if record.levelname == 'ERROR':
            has_log = True
            assert json.loads(record.msg)['event'] == event
    assert has_log


class ProcTwice(Proc):
    @process
    def process(self):
        pass


def test_process_twice(worker, mongo, no_warn):
    p = ProcTwice.create()
    p.process()
    p.block_until_complete()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'process', ProcessStatus.SUCCESS)
