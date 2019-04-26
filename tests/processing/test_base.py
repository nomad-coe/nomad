import pytest
import json
import random

from nomad.processing.base import Proc, process, task, SUCCESS, FAILURE, RUNNING, PENDING

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


def test_tasks(mongo):
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


def test_fail(mongo, with_error):
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


def test_simple_process(worker, mongo, no_warn):
    p = SimpleProc.create()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'two')


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
    assert_proc(p, None)


class TaskInProc(Proc):
    @process
    @task
    def process(self):
        pass


@pytest.mark.timeout(5)
def test_task_as_proc(worker, mongo, no_warn):
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


def test_fail_on_proc_in_proc(worker, mongo):
    p = ProcInProc.create()
    p.one()
    p.block_until_complete()
    assert_proc(p, 'one', FAILURE, 1)
