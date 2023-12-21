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
import time
import threading
from typing import List, Any, Union

from mongoengine import StringField, IntField, ListField

from nomad.processing.base import (
    Proc,
    ProcessAlreadyRunning,
    process,
    process_local,
    ProcessStatus,
)

random.seed(0)


fail = 'FAIL'
events: List[str] = []


@pytest.fixture(scope='function')
def reset_events():
    events.clear()


def assert_proc(
    proc, current_process, process_status=ProcessStatus.SUCCESS, errors=0, warnings=0
):
    assert proc.current_process == current_process
    assert proc.process_status == process_status
    assert len(proc.errors) == errors
    assert len(proc.warnings) == warnings
    assert not proc.process_running


def assert_events(expected_events: List[Union[str, List[str]]]):
    ind = 0
    for expected in expected_events:
        if isinstance(expected, str):
            # str -> expect a specific event
            assert ind <= len(events), f'Not enough events, expecting {expected}'
            assert expected == events[ind]
            ind += 1
        elif isinstance(expected, list):
            # list -> expecting a number of events, in any order
            while expected:
                assert ind <= len(
                    events
                ), f'Not enough events, expecting one of {expected}'
                event = events[ind]
                ind += 1
                assert (
                    event in expected
                ), f'Unexpected event: {event}, expecting one of {expected}'
                expected.remove(event)  # type: ignore
        else:
            assert False, 'Bad value in expected_events'
    assert ind >= len(events), 'Too many events'


class SimpleProc(Proc):
    @process()
    def a_process(self):
        pass

    @process()
    def a_process_with_arguments(
        self, a_string, a_int, a_float, a_bool, a_list, a_tuple, a_dict, **kwargs
    ):
        assert a_string == 'a string'
        assert a_int == 7
        assert a_float == 3.14
        assert a_bool is True
        assert a_list == [1, 2]
        assert a_tuple == [3, 4]  # Mongodb turns tuples to lists
        assert a_dict == {'1': 'one', '2': 2}
        assert kwargs == {'kwarg': 'kwarg_value'}


@pytest.mark.parametrize(
    'with_args', [pytest.param(False, id='no-args'), pytest.param(True, id='with-args')]
)
def test_simple_process(worker, mongo_function, no_warn, with_args):
    p = SimpleProc.create()
    if with_args:
        process = 'a_process_with_arguments'
        p.a_process_with_arguments(
            'a string',
            7,
            3.14,
            True,
            [1, 2],
            (3, 4),
            {'1': 'one', '2': 2},
            kwarg='kwarg_value',
        )
    else:
        process = 'a_process'
        p.a_process()
    p.block_until_complete()
    assert_proc(p, process)


class FailingProc(Proc):
    @process()
    def will_fail_with_exception(self):
        _ = 1 / 0


def test_failing_process(worker, mongo_function, with_error):
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


def test_process_twice(worker, mongo_function, no_warn):
    p = ProcTwice.create()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'process', ProcessStatus.SUCCESS)
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
    def child_proc(self, succeed: bool, new_child_id: str = None, delay=0.1):
        time.sleep(delay)
        if new_child_id:
            # For testing adding new children during processing
            events.append(f'{self.child_id}:child_proc:add_child')
            new_child = ChildProc.create(
                child_id=new_child_id, parent_id=self.parent_id
            )
            new_child.child_proc(True)
        if not succeed:
            events.append(f'{self.child_id}:child_proc:fail')
            assert False, 'failing child'
        events.append(f'{self.child_id}:child_proc:succ')

    def parent(self):
        return ParentProc.get(self.parent_id)


class ParentProc(Proc):
    id_field = 'parent_id'
    parent_id = StringField(primary_key=True)
    # State variables
    current_slot = IntField()
    join_args = ListField()

    @classmethod
    def get(cls, parent_id: str) -> 'ParentProc':
        return cls.get_by_id(parent_id, 'parent_id')

    @process()
    def spawn(
        self,
        fail_spawn: bool = False,
        suffix: str = '',
        delay=0.1,
        child_args: List[Any] = [],
        join_args: List[Any] = [],
    ):
        """
        Arguments:
            fail_spawn: if we should fail after we have spawned the child processes
            child_args: For each value in the list, spawn a child with the provided args
            join_args: list of parameters controlling the behaviour when joining. `fail` mean
                fail the join. A boolean or list of booleans result in new children spawned.
        """
        events.append(f'{self.parent_id}:spawn:start{suffix}')
        self.join_args = join_args
        self.current_slot = 0
        child_count = 0
        for child_arg in child_args:
            if not isinstance(child_arg, list):
                child_arg = [child_arg]
            child = ChildProc.create(
                child_id=str(child_count), parent_id=self.parent_id
            )
            child.child_proc(*child_arg)
            child_count += 1
        if fail_spawn:
            events.append(f'{self.parent_id}:spawn:fail{suffix}')
            assert False, 'failing in spawn'
        time.sleep(delay)
        events.append(f'{self.parent_id}:spawn:waiting{suffix}')
        return ProcessStatus.WAITING_FOR_RESULT

    def child_cls(self):
        return ChildProc

    def join(self):
        if self.join_args:
            if self.current_slot < len(self.join_args):
                # One more join arg to process
                join_arg = self.join_args[self.current_slot]
                self.current_slot += 1
                if join_arg == fail:
                    events.append(f'{self.parent_id}:join:fail')
                    assert False, 'failing in join'
                else:
                    if isinstance(join_arg, bool):
                        join_arg = [join_arg]
                    for i, succeed in enumerate(join_arg):
                        # Start up another child
                        child = ChildProc.create(
                            child_id=f'rejoin{self.current_slot}.{i}',
                            parent_id=self.parent_id,
                        )
                        child.child_proc(succeed)
                    events.append(f'{self.parent_id}:join:waiting')
                    return ProcessStatus.WAITING_FOR_RESULT
        events.append(f'{self.parent_id}:join:succ')

    @process(is_blocking=True)
    def blocking(self, delay=0.1):
        events.append(f'{self.parent_id}:blocking:start')
        time.sleep(delay)
        events.append(f'{self.parent_id}:blocking:succ')

    @process(is_blocking=False)
    def non_blocking(self, delay=0.1):
        events.append(f'{self.parent_id}:non_blocking:start')
        time.sleep(delay)
        events.append(f'{self.parent_id}:non_blocking:succ')

    @process_local
    def local_process(self, delay=0.1, fail=False):
        events.append(f'{self.parent_id}:local_process:start')
        time.sleep(delay)
        if fail:
            events.append(f'{self.parent_id}:local_process:fail')
            assert False, 'Failing local process'
        events.append(f'{self.parent_id}:local_process:succ')


@pytest.mark.parametrize(
    'spawn_kwargs, expected_events',
    [
        pytest.param(
            dict(),
            ['p:spawn:start', 'p:spawn:waiting', 'p:join:succ'],
            id='no-children',
        ),
        pytest.param(
            dict(child_args=[True] * 20),
            [
                'p:spawn:start',
                ['p:spawn:waiting'] + [f'{ci}:child_proc:succ' for ci in range(20)],
                'p:join:succ',
            ],
            id='20-successful',
        ),
        pytest.param(
            dict(child_args=[True, False]),
            [
                'p:spawn:start',
                ['p:spawn:waiting', '0:child_proc:succ', '1:child_proc:fail'],
                'p:join:succ',
            ],
            id='one-succ-one-fail',
        ),
        pytest.param(
            dict(fail_spawn=True, child_args=[True, False]),
            [
                'p:spawn:start',
                ['p:spawn:fail', '0:child_proc:succ', '1:child_proc:fail'],
            ],
            id='fail-spawn',
        ),
        pytest.param(
            dict(child_args=[True, False], join_args=[fail]),
            [
                'p:spawn:start',
                ['p:spawn:waiting', '0:child_proc:succ', '1:child_proc:fail'],
                'p:join:fail',
            ],
            id='join-fail',
        ),
        pytest.param(
            dict(child_args=[True, False], join_args=[True, [True, False], False]),
            [
                'p:spawn:start',
                ['p:spawn:waiting', '0:child_proc:succ', '1:child_proc:fail'],
                'p:join:waiting',
                'rejoin1.0:child_proc:succ',
                'p:join:waiting',
                ['rejoin2.0:child_proc:succ', 'rejoin2.1:child_proc:fail'],
                'p:join:waiting',
                'rejoin3.0:child_proc:fail',
                'p:join:succ',
            ],
            id='join-multiple',
        ),
        pytest.param(
            dict(child_args=[True, False], join_args=[True, [], []]),
            [
                'p:spawn:start',
                ['p:spawn:waiting', '0:child_proc:succ', '1:child_proc:fail'],
                'p:join:waiting',
                'rejoin1.0:child_proc:succ',
                'p:join:waiting',
                'p:join:waiting',
                'p:join:succ',
            ],
            id='join-multiple-no-children',
        ),
        pytest.param(
            dict(child_args=[True, False], join_args=[True, False, fail]),
            [
                'p:spawn:start',
                ['p:spawn:waiting', '0:child_proc:succ', '1:child_proc:fail'],
                'p:join:waiting',
                'rejoin1.0:child_proc:succ',
                'p:join:waiting',
                'rejoin2.0:child_proc:fail',
                'p:join:fail',
            ],
            id='join-multiple-then-fail',
        ),
        pytest.param(
            dict(child_args=[[True, 'new_child']], join_args=[]),
            [
                'p:spawn:start',
                [
                    'p:spawn:waiting',
                    '0:child_proc:add_child',
                    '0:child_proc:succ',
                    'new_child:child_proc:succ',
                ],
                'p:join:succ',
            ],
            id='add-child-while-processing',
        ),
    ],
)
def test_parent_child(
    worker, mongo_function, reset_events, spawn_kwargs, expected_events
):
    child_args = spawn_kwargs.get('child_args', [])
    join_args = spawn_kwargs.get('join_args', [])
    fail_spawn = spawn_kwargs.get('fail_spawn', False)
    parent = ParentProc.create(parent_id='p')
    parent.spawn(**spawn_kwargs)
    parent.block_until_complete()
    for i, succeed in enumerate(child_args):
        child = ChildProc.get(str(i))
        if child.process_running:
            assert fail_spawn  # Otherwise they should all be done
            child.block_until_complete()
        expected_child_status = (
            ProcessStatus.SUCCESS if succeed else ProcessStatus.FAILURE
        )
        assert child.process_status == expected_child_status
    for i, join_arg in enumerate(join_args):
        if join_arg != fail:
            if isinstance(join_arg, bool):
                join_arg = [join_arg]
            for i2, succeed in enumerate(join_arg):
                child = ChildProc.get(f'rejoin{i + 1}.{i2}')
                assert (
                    child.process_status == ProcessStatus.SUCCESS
                    if succeed
                    else ProcessStatus.FAILURE
                )
    expected_parent_status = (
        ProcessStatus.FAILURE
        if fail_spawn or fail in join_args
        else ProcessStatus.SUCCESS
    )
    assert parent.process_status == expected_parent_status
    assert_events(expected_events)


def test_queueing(worker, mongo_function, reset_events):
    p = ParentProc.create(parent_id='p')
    expected_events = []
    # Schedule 20 calls
    for i in range(20):
        # Spawn, with "long" delay for the first process so everything will be queued up
        p.spawn(suffix=f':{i}', child_args=[True], delay=1.0 if i == 0 else 0.1)
        expected_events.extend(
            [
                f'p:spawn:start:{i}',
                [f'p:spawn:waiting:{i}', '0:child_proc:succ'],
                'p:join:succ',
            ]
        )
    assert len(p.queue) >= 19  # The first process may have started
    p.block_until_complete()
    assert p.process_status == ProcessStatus.SUCCESS
    assert_events(expected_events)


def test_queueing_failure(worker, mongo_function, reset_events):
    p = ParentProc.create(parent_id='p')
    # Schedule 20 calls, the second should fail
    for i in range(20):
        # Spawn, with "long" delay for the first process so everything will be queued up
        p.spawn(
            fail_spawn=(i == 1),
            suffix=f':{i}',
            child_args=[],
            delay=1.0 if i == 0 else 0.1,
        )
    assert len(p.queue) >= 19  # The first process may have started
    p.block_until_complete()
    assert p.process_status == ProcessStatus.FAILURE
    # After the failure, the queue will be cleared
    assert_events(
        [
            'p:spawn:start:0',
            'p:spawn:waiting:0',
            'p:join:succ',
            'p:spawn:start:1',
            'p:spawn:fail:1',
        ]
    )


def test_non_blocking_then_blocking(worker, mongo_function, reset_events):
    p = ParentProc.create(parent_id='p')
    p.spawn(delay=1.0)
    p.blocking()
    assert len(p.queue) >= 1  # The first process may have started
    p.block_until_complete()
    assert p.process_status == ProcessStatus.SUCCESS
    assert_events(
        [
            'p:spawn:start',
            'p:spawn:waiting',
            'p:join:succ',
            'p:blocking:start',
            'p:blocking:succ',
        ]
    )


def test_blocking_then_non_blocking(worker, mongo_function, reset_events):
    p = ParentProc.create(parent_id='p')
    p.blocking(delay=1.0)
    with pytest.raises(ProcessAlreadyRunning):
        p.spawn()
    p.block_until_complete()
    assert p.process_status == ProcessStatus.SUCCESS
    assert_events(['p:blocking:start', 'p:blocking:succ'])


def test_local_blocked(worker, mongo_function, reset_events):
    p = ParentProc.create(parent_id='p')

    def other_call():
        time.sleep(0.5)
        p2 = ParentProc.get(parent_id='p')
        try:
            p2.local_process()
        except ProcessAlreadyRunning:
            events.append('other_call:blocked')

    a_thread = threading.Thread(target=other_call)
    a_thread.start()
    p.non_blocking(delay=1.0)
    p.block_until_complete()
    assert p.process_status == ProcessStatus.SUCCESS
    a_thread.join()
    assert_events(['p:non_blocking:start', 'other_call:blocked', 'p:non_blocking:succ'])


def test_local_blocking(worker, mongo_function, reset_events):
    p = ParentProc.create(parent_id='p')

    def other_call():
        time.sleep(0.5)
        p2 = ParentProc.get(parent_id='p')
        try:
            p2.non_blocking()
        except ProcessAlreadyRunning:
            events.append('other_call:blocked1')
        try:
            p2.local_process()
        except ProcessAlreadyRunning:
            events.append('other_call:blocked2')

    a_thread = threading.Thread(target=other_call)
    a_thread.start()
    p.local_process(delay=1.0)
    assert p.process_status == ProcessStatus.SUCCESS
    a_thread.join()
    assert_events(
        [
            'p:local_process:start',
            'other_call:blocked1',
            'other_call:blocked2',
            'p:local_process:succ',
        ]
    )


def test_local_failed(worker, mongo_function, reset_events):
    p = ParentProc.create(parent_id='p')
    try:
        p.local_process(fail=True)
        raise RuntimeError('local_process expected to throw an exception')
    except AssertionError:
        pass  # Expected to happen
    assert p.process_status == ProcessStatus.FAILURE
