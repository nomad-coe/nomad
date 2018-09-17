import pytest
from mongoengine import connect, IntField, ReferenceField
from mongoengine.connection import disconnect
import time
import logging
import json

from nomad import config
from nomad.processing.base import Proc, process, task, SUCCESS, FAILURE, RUNNING, PENDING


def assert_proc(proc, current_task, status=SUCCESS, errors=0, warnings=0):
    assert proc.current_task == current_task
    assert proc.status == status
    assert len(proc.errors) == errors
    assert len(proc.warnings) == warnings


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


def test_tasks(mongomock):
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


def test_fail(caplog):
    caplog.set_level(logging.CRITICAL, logger='nomad.processing.base')
    p = FailTasks.create()
    p.will_fail()

    assert_proc(p, 'will_fail', FAILURE, errors=1)
    has_log = False
    for record in caplog.records:
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


def test_simple_process(worker):
    p = SimpleProc.create()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'two')


class TaskInProc(Proc):
    @process
    @task
    def process(self):
        pass


# @pytest.mark.timeout(5)
def test_task_as_proc(worker):
    p = TaskInProc.create()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'process')


class ParentProc(Proc):
    children = IntField(default=0)

    @process
    @task
    def spawn_children(self):
        ChildProc.create(parent=self).process()

    @process
    @task
    def after_children(self):
        pass

    def on_child_complete(self):
        if self.incr_counter('children') == 1:
            self.after_children()


class ChildProc(Proc):
    parent = ReferenceField(ParentProc)

    @process
    @task
    def process(self):
        self.parent.on_child_complete()


# @pytest.mark.timeout(5)
def test_counter(worker, caplog):
    p = ParentProc.create()
    p.spawn_children()
    p.block_until_complete()
    assert_proc(p, 'after_children')

    for record in caplog.records:
        assert record.levelname not in ['WARNING', 'ERROR', 'CRITICAL']
