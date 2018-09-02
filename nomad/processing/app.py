# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import List
from contextlib import contextmanager
import inspect
import logging
import celery
from celery import Task
from celery.signals import after_setup_task_logger, after_setup_logger, worker_process_init
from mongoengine import Document, StringField, ListField, DateTimeField, IntField, \
    ReferenceField, connect
from pymongo import ReturnDocument
from datetime import datetime

from nomad import config, utils
import nomad.patch  # pylint: disable=unused-import


class Celery(celery.Celery):

    def mongo_connect(self):
        return connect(db=config.mongo.users_db, host=config.mongo.host)

    def __init__(self):
        if config.logstash.enabled:
            def initialize_logstash(logger=None, loglevel=logging.DEBUG, **kwargs):
                utils.add_logstash_handler(logger)
                return logger

            after_setup_task_logger.connect(initialize_logstash)
            after_setup_logger.connect(initialize_logstash)

        worker_process_init.connect(lambda **kwargs: self.mongo_connect())

        super().__init__(
            'nomad.processing', backend=config.celery.backend_url, broker=config.celery.broker_url)

        self.add_defaults(dict(
            accept_content=['json', 'pickle'],
            task_serializer=config.celery.serializer,
            result_serializer=config.celery.serializer,
        ))


app = Celery()


PENDING = 'PENDING'
RUNNING = 'RUNNING'
FAILURE = 'FAILURE'
SUCCESS = 'SUCCESS'


class InvalidId(Exception): pass


class AsyncDocumentNotRegistered(Exception): pass


class Proc(Document):
    """
    Base class for objects involved in processing and need persistent processing
    state.

    Attributes:
        current_task: the currently running or last completed task
        status: the overall status of the processing
        errors: a list of errors that happened during processing. Error fail a processing
            run
        warnings: a list of warnings that happened during processing. Warnings do not
            fail a processing run
        create_time: the time of creation (not the start of processing)
        proc_time: the time that processing completed (successfully or not)
    """

    meta = {
        'abstract': True,
    }

    tasks: List[str] = None
    """ the ordered list of tasks that comprise a processing run """

    current_task = StringField(default=None)
    status = StringField(default='CREATED')

    errors = ListField(StringField, default=[])
    warnings = ListField(StringField, default=[])

    create_time = DateTimeField(required=True)
    complete_time = DateTimeField()

    @property
    def completed(self) -> bool:
        return self.status in [SUCCESS, FAILURE]

    def __init__(self, **kwargs):
        assert self.__class__.tasks is not None and len(self.tasks) > 0, \
            """ the class attribute tasks must be overwritten with an acutal list """
        assert 'status' not in kwargs, \
            """ do not set the status manually, its managed """

        kwargs.setdefault('create_time', datetime.now())
        super().__init__(**kwargs)
        self.status = PENDING if self.current_task is None else RUNNING

    @classmethod
    def get(cls, obj_id):
        try:
            obj = cls.objects(id=obj_id).first()
        except Exception as e:
            raise InvalidId('%s is not a valid id' % obj_id)

        if obj is None:
            raise KeyError('%s with id %s does not exist' % (cls.__name__, obj_id))

        return obj

    def get_id(self):
        return self.id.__str__()

    def fail(self, error):
        assert not self.completed

        self.status = FAILURE
        self.errors = [str(error) for error in errors]
        self.complete_time = datetime.now()

        self.save()

    def warning(self, *warnings):
        assert not self.completed

        for warning in warnings:
            self.warnings.append(str(warning))

    def continue_with(self, task):
        assert task in self.tasks
        assert self.tasks.index(task) == self.tasks.index(self.current_task) + 1

        if self.status == FAILURE:
            return False

        if self.status == PENDING:
            assert self.current_task is None
            assert task == self.tasks[0]
            self.status = RUNNING

        self.current_task = task
        self.save()
        return True

    def complete(self):
        if self.status != FAILURE:
            assert self.status == RUNNING
            self.status = SUCCESS
            self.save()

    @property
    def json_dict(self) -> dict:
        """ A json serializable dictionary representation. """
        data = {
            'tasks': self.tasks,
            'current_task': self.current_task,
            'status': self.status,
            'errors': self.errors,
            'warnings': self.warnings,
            'create_time': self.create_time.isoformat() if self.create_time is not None else None,
            'complete_time': self.complete_time.isoformat() if self.complete_time is not None else None,
        }
        return {key: value for key, value in data.items() if value is not None}


def task(func):
    def wrapper(self, *args, **kwargs):
        if self.status == 'FAILURE':
            return

        self.continue_with(func.__name__)
        try:
            func(self, *args, **kwargs)
        except Exception as e:
            self.fail(e)

    return wrapper


def process(func):
    @app.task(bind=True, name=func.__name__, ignore_results=True)
    def task_func(task, cls_name, self_id):
        all_cls = AsyncDocument.__subclasses__()
        cls = next((cls for cls in all_cls if cls.__name__ == cls_name), None)
        if cls is None:
            raise AsyncDocumentNotRegistered('Document type %s not registered for async methods' % cls_name)

        self = cls.get(self_id)

        try:
            func(self)
        except Exception as e:
            self.fail(e)
            raise e

    def wrapper(self, *args, **kwargs):
        assert len(args) == 0 and len(kwargs) == 0
        self.save()

        self_id = self.get_id()
        cls_name = self.__class__.__name__

        return task_func.s(cls_name, self_id).delay()

    return wrapper


class Upload(Proc):

    data = StringField(default='Hello, World')
    processed_calcs = IntField(default=0)
    total_calcs = IntField(default=-1)

    @process
    def after_upload(self):
        self.extract()
        self.parse()

    @process
    def after_parse(self):
        self.cleanup()

    @task
    def extract(self):
        print('now extracting')

    @task
    def parse(self):
        Calc(upload=self).parse()
        Calc(upload=self).parse()
        self.total_calcs = 2
        self.save()
        self.check_calcs_complete(more_calcs=0)

    @task
    def cleanup(self):
        print('cleanup')

    def check_calcs_complete(self, more_calcs=1):
        # use a primitive but atomic pymongo call to increase the number of calcs
        updated_raw = Upload._get_collection().find_one_and_update(
            {'_id': self.id},
            {'$inc': {'processed_calcs': more_calcs}},
            return_document=ReturnDocument.AFTER)

        updated_total_calcs = updated_raw['total_calcs']
        updated_processed_calcs = updated_raw['processed_calcs']

        print('%d:%d' % (updated_processed_calcs, updated_total_calcs))
        if updated_processed_calcs == updated_total_calcs and updated_total_calcs != -1:
            self.after_parse()


class Calc(Proc):

    upload = ReferenceField(Upload, required=True)

    @process
    def parse(self):
        print('parsee')
        self.upload.check_calcs_complete()


if __name__ == '__main__':
    connect(db=config.mongo.users_db, host=config.mongo.host)
    tds = [Upload(), Upload(), Upload()]
    for td in tds:
        td.after_upload()
