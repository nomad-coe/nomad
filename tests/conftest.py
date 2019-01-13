import pytest
import logging
from sqlalchemy.orm import Session
from mongoengine import connect
from mongoengine.connection import disconnect

from nomad import config, infrastructure


@pytest.fixture(scope="session")
def monkeysession(request):
    from _pytest.monkeypatch import MonkeyPatch
    mpatch = MonkeyPatch()
    yield mpatch
    mpatch.undo()


@pytest.fixture(scope='session', autouse=True)
def nomad_logging():
    config.logstash = config.logstash._replace(enabled=False)
    config.console_log_level = logging.CRITICAL
    infrastructure.setup_logging()


@pytest.fixture(scope='session')
def celery_includes():
    return ['nomad.processing.base']


@pytest.fixture(scope='session')
def celery_config():
    return {
        'broker_url': config.celery.broker_url
    }


@pytest.fixture(scope='session')
def purged_app(celery_session_app):
    """
    Purges all pending tasks of the celery app before test. This is necessary to
    remove tasks from the queue that might be 'left over' from prior tests.
    """
    celery_session_app.control.purge()
    yield celery_session_app


@pytest.fixture()
def patched_celery(monkeypatch):
    # There is a bug in celery, which prevents to use the celery_worker for multiple
    # tests: https://github.com/celery/celery/issues/4088
    # The bug has a fix from Aug 2018, but it is not yet released (TODO).
    # We monkeypatch a similar solution here.
    def add_reader(self, fds, callback, *args):
        from kombu.utils.eventio import ERR, READ, poll

        if self.poller is None:
            self.poller = poll()

        return self.add(fds, callback, READ | ERR, args)

    monkeypatch.setattr('kombu.asynchronous.hub.Hub.add_reader', add_reader)
    yield


@pytest.fixture(scope='session')
def celery_inspect(purged_app):
    yield purged_app.control.inspect()


@pytest.fixture()
def worker(patched_celery, celery_inspect, celery_session_worker):
    """
    Extension of the celery_session_worker fixture that ensures a clean task queue.
    """
    yield

    # wait until there no more active tasks, to leave clean worker and queues for the next
    # test.
    while True:
        empty = True
        for value in celery_inspect.active().values():
            empty = empty and len(value) == 0
        if empty:
            break


@pytest.fixture(scope='function')
def mockmongo(monkeypatch):

    disconnect()
    connection = connect('test_db', host='mongomock://localhost')
    monkeypatch.setattr('nomad.infrastructure.setup_mongo', lambda **kwargs: None)

    yield

    connection.drop_database('test_db')


@pytest.fixture(scope='session')
def elastic():
    infrastructure.setup_elastic()
    assert infrastructure.elastic_client is not None


@pytest.fixture(scope='session')
def repository_db(monkeysession):
    infrastructure.setup_repository_db()
    assert infrastructure.repository_db_conn is not None

    # we use a transaction around the session to rollback anything that happens within
    # test execution
    trans = infrastructure.repository_db_conn.begin()
    session = Session(bind=infrastructure.repository_db_conn, autocommit=True)
    monkeysession.setattr('nomad.infrastructure.repository_db', session)
    yield infrastructure.repository_db
    trans.rollback()
    session.close()


@pytest.fixture(scope='function')
def clean_repository_db(repository_db):
    # do not wonder, this will not setback the id counters
    repository_db.execute('TRUNCATE uploads CASCADE;')
    yield repository_db


@pytest.fixture(scope='session')
def test_user(repository_db):
    from nomad import coe_repo
    return coe_repo.ensure_test_user(email='sheldon.cooper@nomad-fairdi.tests.de')


@pytest.fixture(scope='session')
def other_test_user(repository_db):
    from nomad import coe_repo
    return coe_repo.ensure_test_user(email='leonard.hofstadter@nomad-fairdi.tests.de')


@pytest.fixture(scope='session')
def admin_user(repository_db):
    from nomad import coe_repo
    return coe_repo.admin_user()


# @pytest.fixture(scope='function')
# def mocksearch(monkeypatch):
#     uploads_by_id = {}
#     by_archive_id = {}

#     def persist(calc):
#         uploads_by_id.setdefault(calc.upload_id, []).append(calc)
#         by_archive_id[calc.calc_id] = calc

#     def upload_exists(self):
#         return self.upload_id in uploads_by_id

#     def upload_delete(self):
#         upload_id = self.upload_id
#         if upload_id in uploads_by_id:
#             for calc in uploads_by_id[upload_id]:
#                 del(by_archive_id[calc.calc_id])
#             del(uploads_by_id[upload_id])

#     @property
#     def upload_calcs(self):
#         return uploads_by_id.get(self.upload_id, [])

#     monkeypatch.setattr('nomad.repo.RepoCalc.persist', persist)
#     monkeypatch.setattr('nomad.repo.RepoUpload.exists', upload_exists)
#     monkeypatch.setattr('nomad.repo.RepoUpload.delete', upload_delete)
#     monkeypatch.setattr('nomad.repo.RepoUpload.calcs', upload_calcs)
#     monkeypatch.setattr('nomad.repo.RepoUpload.unstage', lambda *args, **kwargs: None)

#     return by_archive_id


@pytest.fixture(scope='function')
def no_warn(caplog):
    yield caplog
    for record in caplog.get_records(when='call'):
        if record.levelname in ['WARNING', 'ERROR', 'CRITICAL']:
            assert False, record.msg


@pytest.fixture(scope='function')
def with_error(caplog):
    yield caplog
    count = 0
    for record in caplog.get_records(when='call'):
        if record.levelname in ['ERROR', 'CRITICAL']:
            count += 1

    assert count > 0
