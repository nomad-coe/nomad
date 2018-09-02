import pytest
from mongoengine import connect, IntField, ReferenceField
from mongoengine.connection import disconnect
import time

from nomad import config
from nomad.proc import Proc, process, task, SUCCESS, FAILURE, RUNNING, PENDING


@pytest.fixture(scope='session')
def celery_includes():
    return ['nomad.proc']


@pytest.fixture(scope='session')
def celery_config():
    return {
        'broker_url': config.celery.broker_url,
        'result_backend': config.celery.backend_url,
        'accept_content': ['json', 'pickle'],
        'task_serializer': config.celery.serializer,
        'result_serializer': config.celery.serializer
    }


@pytest.fixture(scope='function')
def mongomock(monkeypatch):
    def mock_connect(**kwargs):
        return connect('test_db', host='mongomock://localhost')

    disconnect()
    connection = mock_connect()
    monkeypatch.setattr('nomad.proc.mongo_connect', mock_connect)
    yield
    connection.drop_database('test_db')


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


def test_fail(mongomock):
    p = FailTasks.create()
    p.will_fail()

    assert_proc(p, 'will_fail', FAILURE, errors=1)


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


def test_simple_process(celery_session_worker, mongomock):
    p = SimpleProc.create()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'two')


class TaskInProc(Proc):
    @process
    @task
    def process(self):
        pass


def test_task_as_proc(celery_session_worker, mongomock):
    p = TaskInProc.create()
    p.process()
    p.block_until_complete()
    assert_proc(p, 'process')


class ParentProc(Proc):
    children = IntField(default=0)

    @process
    @task
    def spawn_children(self):
        print('## spawn')
        ChildProc.create(parent=self).process()

    @process
    @task
    def after_children(self):
        print('## after')
        pass

    def on_child_complete(self):
        print('## on complete')
        if self.incr_counter('children') == 1:
            self.after_children()


class ChildProc(Proc):
    parent = ReferenceField(ParentProc)

    @process
    @task
    def process(self):
        self.parent.on_child_complete()


def test_counter(celery_session_worker, mongomock):
    p = ParentProc.create()
    p.spawn_children()
    p.block_until_complete()
    assert_proc(p, 'after_children')
