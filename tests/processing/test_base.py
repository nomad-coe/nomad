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
from typing import List, Any

from mongoengine import StringField, IntField, ListField

from nomad.processing.base import Proc, process, ProcessStatus

random.seed(0)


_fail = 'FAIL'
_join_count = 0


def assert_proc(proc, current_process, process_status=ProcessStatus.SUCCESS, errors=0, warnings=0):
    assert proc.current_process == current_process
    assert proc.process_status == process_status
    assert len(proc.errors) == errors
    assert len(proc.warnings) == warnings
    assert not proc.process_running


class SimpleProc(Proc):
    @process()
    def a_process(self):
        pass

    @process()
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
    @process()
    def will_fail_with_exception(self):
        _ = 1 / 0


def test_failing_process(worker, mongo, with_error):
    p = FailingProc.create()

    event = 'process failed with exception'
    process = 'will_fail_with_exception'
    p.will_fail_with_exception()
    p.block_until_complete()
    assert_proc(p, process, ProcessStatus.FAILURE, errors=1)

    has_log = False
    for record in with_error.get_records(when='call'):
        if record.levelname == 'ERROR':
            has_log = True
            assert json.loads(record.msg)['event'] == event
    assert has_log


class ProcTwice(Proc):
    @process()
    def process(self):
        pass


def test_process_twice(worker, mongo, no_warn):
    p = ProcTwice.create()
    p.process()
    p.block_until_complete()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'process', ProcessStatus.SUCCESS)


class ChildProc(Proc):
    child_id = StringField(primary_key=True)
    parent_id = StringField(required=True)

    @classmethod
    def get(cls, id: str) -> 'ChildProc':
        return cls.get_by_id(id, 'child_id')

    @process(is_child=True)
    def child_proc(self, succeed: bool):
        assert succeed, 'failing child'

    def parent(self):
        return ParentProc.get(self.parent_id)


class ParentProc(Proc):
    id_field = 'parent_id'
    parent_id = StringField(primary_key=True)
    # State variables
    current_slot = IntField()
    join_args = ListField()

    @classmethod
    def get(cls, id: str) -> 'ParentProc':
        return cls.get_by_id(id, 'parent_id')

    @process()
    def spawn(self, fail_spawn: bool = False, child_args: List[bool] = [], join_args: List[Any] = []):
        self.join_args = join_args
        self.current_slot = 0
        child_count = 0
        for succeed in child_args:
            child = ChildProc.create(child_id=str(child_count), parent_id=self.parent_id)
            child.child_proc(succeed)
            child_count += 1
        assert not fail_spawn, 'failing in spawn'
        return ProcessStatus.WAITING_FOR_RESULT

    def child_cls(self):
        return ChildProc

    def join(self):
        global _join_count
        _join_count += 1
        if self.join_args:
            if self.current_slot < len(self.join_args):
                # One more join arg to process
                join_arg = self.join_args[self.current_slot]
                self.current_slot += 1
                if join_arg == _fail:
                    assert False, 'failing in join'
                else:
                    if join_arg is not None:
                        # Start up another child
                        child = ChildProc.create(
                            child_id=f'rejoin{self.current_slot}', parent_id=self.parent_id)
                        child.child_proc(join_arg)
                    return ProcessStatus.WAITING_FOR_RESULT


@pytest.mark.parametrize('spawn_kwargs', [
    pytest.param(dict(), id='no-children'),
    pytest.param(dict(child_args=[True] * 20), id='20-successful'),
    pytest.param(dict(child_args=[True, False]), id='one-succ-one-fail'),
    pytest.param(dict(fail_spawn=True, child_args=[True, False]), id='fail-spawn'),
    pytest.param(dict(child_args=[True, False], join_args=[_fail]), id='join-fail'),
    pytest.param(dict(child_args=[True, False], join_args=[True, False, False]), id='join-multiple'),
    pytest.param(dict(child_args=[True, False], join_args=[True, None, None]), id='join-multiple-no-children'),
    pytest.param(dict(child_args=[True, False], join_args=[True, False, _fail]), id='join-multiple-then-fail')])
def test_parent_child(worker, mongo, spawn_kwargs):
    global _join_count
    _join_count = 0
    child_args = spawn_kwargs.get('child_args', [])
    join_args = spawn_kwargs.get('join_args', [])
    fail_spawn = spawn_kwargs.get('fail_spawn', False)
    fail_join = _fail in join_args
    parent = ParentProc.create(parent_id='the_parent')
    parent.spawn(**spawn_kwargs)
    parent.block_until_complete()
    for i, succeed in enumerate(child_args):
        child = ChildProc.get(str(i))
        if child.process_running:
            assert fail_spawn  # Otherwise they should all be done
            child.block_until_complete()
        expected_child_status = ProcessStatus.SUCCESS if succeed else ProcessStatus.FAILURE
        assert child.process_status == expected_child_status
    for i, join_arg in enumerate(join_args):
        if join_arg in (True, False):
            child = ChildProc.get(f'rejoin{i + 1}')
            expected_child_status = ProcessStatus.SUCCESS if join_arg else ProcessStatus.FAILURE
            assert child.process_status == expected_child_status
        else:
            try:
                ChildProc.get(f'rejoin{i + 1}')
                assert False, 'Child should not exist'
            except KeyError:
                pass
    expected_join_count = 0 if fail_spawn else 1 + len(join_args)
    if fail_join:
        expected_join_count -= 1
    assert _join_count == expected_join_count
    expected_parent_status = ProcessStatus.FAILURE if fail_spawn or fail_join else ProcessStatus.SUCCESS
    assert parent.process_status == expected_parent_status
