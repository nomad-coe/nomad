import os
import shutil
import tempfile
import time

import elasticsearch
import elasticsearch.exceptions
import pytest

from nomad import infrastructure
from nomad.config import config

elastic_test_entries_index = 'nomad_entries_v1_test'
elastic_test_materials_index = 'nomad_materials_v1_test'

indices = [elastic_test_entries_index, elastic_test_materials_index]


@pytest.fixture(scope='session', autouse=True)
def raw_files_infra():
    parent_directory = '.volumes'
    if not os.path.isdir(parent_directory):
        os.makedirs(parent_directory, exist_ok=True)
    directory = tempfile.TemporaryDirectory(dir=parent_directory, prefix='test_fs')
    config.fs.tmp = tempfile.TemporaryDirectory(dir=directory.name, prefix='tmp').name
    config.fs.staging = tempfile.TemporaryDirectory(
        dir=directory.name, prefix='staging'
    ).name
    config.fs.public = tempfile.TemporaryDirectory(
        dir=directory.name, prefix='public'
    ).name
    config.fs.staging_external = os.path.abspath(config.fs.staging)
    config.fs.public_external = os.path.abspath(config.fs.public)
    config.fs.prefix_size = 2
    clear_raw_files()
    yield
    directory.cleanup()


@pytest.fixture(scope='module')
def raw_files_module(raw_files_infra):
    """Provides cleaned out files directory structure per module. Clears files before test."""
    clear_raw_files()


@pytest.fixture(scope='function')
def raw_files_function(raw_files_infra):
    """Provides cleaned out files directory structure per function. Clears files before test."""
    clear_raw_files()


def clear_raw_files():
    directories = [config.fs.staging, config.fs.public, config.fs.tmp]
    for directory in directories:
        try:
            shutil.rmtree(directory)
        except FileNotFoundError:
            pass

        os.makedirs(directory)


@pytest.fixture(scope='session')
def celery_includes():
    return ['nomad.processing.base']


@pytest.fixture(scope='session')
def celery_config():
    return {'broker_url': config.rabbitmq_url(), 'task_queue_max_priority': 10}


@pytest.fixture(scope='session')
def purged_app(celery_session_app):
    """
    Purges all pending tasks of the celery app before test. This is necessary to
    remove tasks from the queue that might be 'left over' from prior tests.
    """
    celery_session_app.control.purge()
    yield celery_session_app
    celery_session_app.control.purge()


@pytest.fixture(scope='session')
def celery_inspect(purged_app, pytestconfig):
    timeout = pytestconfig.getoption('celery_inspect_timeout')
    yield purged_app.control.inspect(timeout=timeout)


# It might be necessary to make this a function scoped fixture, if old tasks keep
# 'bleeding' into successive tests.
@pytest.fixture(scope='function')
def worker(mongo_function, celery_session_worker, celery_inspect):
    """Provides a clean worker (no old tasks) per function. Waits for all tasks to be completed."""
    yield

    # wait until there no more active tasks, to leave clean worker and queues for the next
    # test run.
    try:
        while True:
            empty = True
            celery_active = celery_inspect.active()
            if not celery_active:
                break
            for value in celery_active.values():
                empty = empty and len(value) == 0
            if empty:
                break
    except Exception:
        print('Exception during worker tear down.')
        import traceback

        traceback.print_exc()


@pytest.fixture(scope='session')
def mongo_infra(monkeysession):
    monkeysession.setattr('nomad.config.mongo.db_name', 'test_db')
    # disconnecting and connecting again results in an empty database with mongomock
    monkeysession.setattr('mongoengine.disconnect', lambda *args, **kwargs: None)
    return infrastructure.setup_mongo()


def clear_mongo(mongo_infra):
    # Some test cases need to reset the database connection
    infrastructure.mongo_client.drop_database('test_db')
    return infrastructure.mongo_client


@pytest.fixture(scope='module')
def mongo_module(mongo_infra):
    """Provides a cleaned mocked mongo per module."""
    return clear_mongo(mongo_infra)


@pytest.fixture(scope='function')
def mongo_function(mongo_infra):
    """Provides a cleaned mocked mongo per function."""
    return clear_mongo(mongo_infra)


@pytest.fixture(scope='session')
def elastic_infra(monkeysession):
    """Provides elastic infrastructure to the session"""
    monkeysession.setattr(
        'nomad.config.elastic.entries_index', elastic_test_entries_index
    )
    monkeysession.setattr(
        'nomad.config.elastic.materials_index', elastic_test_materials_index
    )

    # attempt to remove and recreate all indices
    return clear_elastic_infra()


def clear_elastic_infra():
    """
    Removes and re-creates all indices and mappings.
    """
    from elasticsearch_dsl import connections

    connection = connections.create_connection(
        hosts=['%s:%d' % (config.elastic.host, config.elastic.port)]
    )

    for index in indices:
        try:
            connection.indices.delete(index=index)
        except Exception:
            pass

    return infrastructure.setup_elastic()


def clear_elastic(elastic_infra):
    """
    Removes all contents from the existing indices.
    """
    try:
        for index in indices:
            retry_count = 10
            while True:
                try:
                    elastic_infra.delete_by_query(
                        index=index,
                        body=dict(query=dict(match_all={})),
                        wait_for_completion=True,
                        refresh=True,
                    )
                    break  # Success! Break the retry loop
                except elasticsearch.exceptions.ConflictError:
                    if retry_count:
                        # Sleep and try again
                        time.sleep(0.1)
                        retry_count -= 1
                    else:
                        raise

    except elasticsearch.exceptions.NotFoundError:
        # Happens if a test removed indices without recreating them.
        clear_elastic_infra()

    assert infrastructure.elastic_client is not None
    return elastic_infra


@pytest.fixture(scope='module')
def elastic_module(elastic_infra):
    """Provides a clean elastic per module. Clears elastic before test."""
    return clear_elastic(elastic_infra)


@pytest.fixture(scope='function')
def elastic_function(elastic_infra):
    """Provides a clean elastic per function. Clears elastic before test."""
    return clear_elastic(elastic_infra)


@pytest.fixture
def reset_infra(mongo_function, elastic_function):
    """Fixture that resets infrastructure after deleting db or search index."""
    yield None
