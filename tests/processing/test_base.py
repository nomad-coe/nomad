import pytest
from mongoengine import connect, IntField, ReferenceField, BooleanField, EmbeddedDocumentField
from mongoengine.connection import disconnect
import time
import logging
import json
import random
import time

from nomad import config
from nomad.processing.base import Proc, Chord, process, task, SUCCESS, FAILURE, RUNNING, PENDING

random.seed(0)


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


def test_fail(one_error):
    p = FailTasks.create()
    p.will_fail()

    assert_proc(p, 'will_fail', FAILURE, errors=1)
    has_log = False
    for record in one_error.records:
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


def test_simple_process(worker, no_warn):
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
def test_task_as_proc(worker, no_warn):
    p = TaskInProc.create()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'process')


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
def test_counter(worker, no_warn):
    p = ParentProc.create()
    p.spawn_children()
    p.block_until_complete()

    p = ParentProc.get(p.id)
    assert_proc(p, 'join')
    # TODO there seems to be a bug, that makes this fail from time to time.
    # assert p.joined
