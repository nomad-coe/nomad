import pytest
from mongoengine import ReferenceField
import time
import json
import random
import time

from nomad.processing.base import Proc, Chord, process, task, SUCCESS, FAILURE, RUNNING, PENDING

random.seed(0)


def assert_proc(proc, current_task, tasks_status=SUCCESS, errors=0, warnings=0):
    assert proc.current_task == current_task
    assert proc.tasks_status == tasks_status
    assert len(proc.errors) == errors
    assert len(proc.warnings) == warnings
    assert not proc.process_running


class Tasks(Proc):
    @task
    def a(self):
        pass

    @task
    def b(self):
        pass


class SingleTask(Proc):
    @task
    def single(self):
        pass


def test_tasks(mockmongo):
    p = Tasks.create()
    assert p.tasks == ['a', 'b']
    assert_proc(p, None, PENDING)

    p.a()
    assert_proc(p, 'a', RUNNING)

    p.b()
    assert_proc(p, 'b')

    p = SingleTask.create()
    p.single()
    assert_proc(p, 'single')


class FailTasks(Proc):
    @task
    def will_fail(self):
        self.fail('fail fail fail')


def test_fail(mockmongo, with_error):
    p = FailTasks.create()
    p.will_fail()

    assert_proc(p, 'will_fail', FAILURE, errors=1)
    has_log = False
    for record in with_error.get_records(when='call'):
        if record.levelname == 'ERROR':
            has_log = True
            assert json.loads(record.msg)['event'] == 'task failed'
    assert has_log


class SimpleProc(Proc):
    @process
    def process(self):
        self.one()
        self.two()

    @task
    def one(self):
        pass

    @task
    def two(self):
        pass


def test_simple_process(mockmongo, worker, no_warn):
    p = SimpleProc.create()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'two')


class TaskInProc(Proc):
    @process
    @task
    def process(self):
        pass


@pytest.mark.timeout(5)
def test_task_as_proc(mockmongo, worker, no_warn):
    p = TaskInProc.create()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'process')


class ProcInProc(Proc):
    @process
    @task
    def one(self):
        self.two()

    @process
    @task
    def two(self):
        pass


def test_fail_on_proc_in_proc(mockmongo, worker):
    p = ProcInProc.create()
    p.one()
    p.block_until_complete()
    assert_proc(p, 'one', FAILURE, 1)


class ParentProc(Chord):

    @process
    @task
    def spawn_children(self):
        count = 23
        for _ in range(0, count):
            ChildProc.create(parent=self).process()

        self.spwaned_childred(count)

    @task
    def join(self):
        pass


class ChildProc(Proc):
    parent = ReferenceField(ParentProc)

    @process
    @task
    def process(self):
        time.sleep(random.uniform(0, 0.1))
        self.parent.completed_child()


@pytest.mark.timeout(10)
def test_counter(mockmongo, worker, no_warn):
    p = ParentProc.create()
    p.spawn_children()
    p.block_until_complete()

    p = ParentProc.get(p.id)
    assert_proc(p, 'join')
    # TODO there seems to be a bug, that makes this fail from time to time.
    # assert p.joined
