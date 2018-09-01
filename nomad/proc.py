from typing import List
from mongoengine import EmbeddedDocument, Document, StringField, ListField, DateTimeField, \
    EmbeddedDocumentField, IntField, ReferenceField
from mongoengine.errors import ValidationError
from datetime import datetime
from celery import Celery
from celery.signals import after_setup_task_logger, after_setup_logger
import logging

from nomad import utils, files, config, parsing, normalizing
import nomad.patch  # pylint: disable=unused-import

if config.logstash.enabled:
    def initialize_logstash(logger=None, loglevel=logging.DEBUG, **kwargs):
        utils.add_logstash_handler(logger)
        return logger

    after_setup_task_logger.connect(initialize_logstash)
    after_setup_logger.connect(initialize_logstash)

app = Celery('nomad.proc', broker=config.celery.broker_url)


PENDING = 'PENDING'
RUNNING = 'RUNNING'
FAILURE = 'FAILURE'
SUCCESS = 'SUCCESS'


class Proc(EmbeddedDocument):
    current_task = StringField(default=None)
    tasks = ListField(StringField, required=True)
    status = StringField(default=PENDING)
    errors = ListField(StringField, default=[])
    warnings = ListField(StringField, default=[])

    create_time = DateTimeField(required=True)
    complete_time = DateTimeField()

    def __init__(self, **kwargs):
        kwargs.setdefault('create_time', datetime.now())
        super().__init__(**kwargs)
        assert self.current_task is None
        assert len(self.tasks) > 0

    @property
    def completed(self) -> bool:
        return self.status in [SUCCESS, FAILURE]

    def continue_with(self, task: str):
        assert task in self.tasks
        assert self.tasks.index(task) == self.tasks.index(self.current_task) + 1

        if self.status == PENDING:
            assert self.current_task is None
            assert task == self.tasks[0]
            self.status = RUNNING

        self.current_task = task

    def complete(self):
        assert self.current_task == self.tasks[-1]

        self.status = SUCCESS
        self.complete_time = datetime.now()

    def fail(self, *errors):
        assert not self.completed

        self.status = FAILURE
        self.errors = [str(error) for error in errors]
        self.complete_time = datetime.now()

    def warning(self, *warnings):
        assert not self.completed

        for warning in warnings:
            self.warnings.append(str(warning))

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


class InvalidId(Exception): pass


UPLOADING = 'uploading'
EXTRACTING = 'extracting'
PARSING = 'parsing'
CLEANING = 'clearning'


class BaseDocument(Document):
    meta = {
        'abstract': True,
    }

    @classmethod
    def get(cls, id: str):
        try:
            obj = cls.objects(id=id).first()
        except ValidationError:
            raise InvalidId('Invalid %s id' % cls.__name__)

        if obj is None:
            raise KeyError('%s does not exist' % cls.__name__)

        return obj


@app.task(bind=True, ignore_result=True)
def process_upload(task, upload_id):
    self = Upload.get(upload_id)
    logger = utils.get_logger(__name__, task=task.name, upload_id=upload_id)

    self.proc.continue_with(EXTRACTING)
    self.save()

    try:
        upload = files.Upload(upload_id)
        upload.open()
        logger.debug('Opened upload')
    except KeyError as e:
        logger.info('Process request for non existing upload')
        self.proc.fail(e)
        return
    except files.UploadError as e:
        logger.info('Could not open upload', error=str(e))
        self.proc.fail(e)
        return
    except Exception as e:
        logger.error('Unknown exception', exc_info=e)
        self.proc.fail(e)
        return
    finally:
        self.save()

    try:
        self.upload_hash = upload.hash()
    except files.UploadError as e:
        logger.error('Could not create upload hash', error=str(e))
        self.proc.fail(e)
        return
    finally:
        self.save()

    try:
        from nomad.data import Calc
        if Calc.upload_exists(proc.upload_hash):
            logger.info('Upload hash doublet')
            self.proc.fail('The same file was already uploaded and processed.')
            return
    except Exception as e:
        logger.error('Exception while checking upload hash', exc_info=e)
        self.proc.fail('Could not check upload hash for existing calcs.', e)
    finally:
        self.save()

    try:
        # TODO: deal with multiple possible parser specs
        for filename in upload.filelist:
            for parser in parsers:
                try:
                    if parser.is_mainfile(filename, lambda fn: upload.open_file(fn)):
                        tmp_mainfile = upload.get_path(filename)
                        self.add_calc(filename, parser.name, tmp_mainfile)
                except Exception as e:
                    logger.warning('Exception while matching pot. mainfile.', mainfile=filename)
                    self.proc.warning(
                        'Exception while matching pot. mainfile %s with parser %s.'
                        % (filename, parser.name))
        self.start_parsing()
    except Exception as e:
        logger.error('Exception while finding parse specs.', exc_info=e)
        self.proc.fail('Exception while fining parse specs.', e)
        self.save()


@app.task(bind=True, ignore_result=True)
def clean_upload(task, upload_id: str):
    logger = utils.get_logger(__name__, task=task.name, upload_id=upload_id)
    self = Upload.get(upload_id)

    self.continue_with(CLEANING)

    try:
        upload = files.Upload(upload_id)
    except KeyError as e:
        logger.warn('Upload does not exist')
        self.fail(e)
        return
    finally:
        self.save()

    try:
        upload.close()
        logger.debug('Closed upload')
    except Exception as e:
        logger.error('Could not close upload', exc_info=e)
        self.fail(e)
        self.save()
        return

    self.complete()
    self.save()


@app.task(bind=True, ignore_result=True)
def process_calc(task, calc_id: str):
    self: Calc = Calc.get(calc_id)

    upload_hash = self.upload.upload_hash
    parser, mainfile = self.parser, self.mainfile

    logger = utils.get_logger(
        __name__, task=task.name,
        upload_id=self.upload.upload_id, upload_hash=upload_hash, mainfile=mainfile)

    # parsing
    self.proc.continue_with(parser)
    self.save()
    try:
        parser_backend = parsing.parser_dict[parser].run(self.tmp_mainfile)
        if parser_backend.status[0] != 'ParseSuccess':
            error = parser_backend.status[1]
            logger.debug('Failed parsing', parser=parser, error=error)
            self.proc.fail(error)
            self.save()
            return
        logger.debug('Completed successfully', parser=parser)
    except Exception as e:
        logger.warn('Exception wile parsing', parser=parser, exc_info=e)
        self.proc.fail(e)
        self.save()
        return

    # normalization
    for normalizer in normalizing.normalizers:
        normalizer_name = normalizer.__name__
        self.proc.continue_with(normalizer_name)
        self.save()
        try:
            normalizer(parser_backend).normalize()
            if parser_backend.status[0] != 'ParseSuccess':
                error = parser_backend.status[1]
                logger.info('Failed run of %s: %s' % (normalizer, error))
                self.proc.fail(error)
                self.save()
                return
            logger.debug('Completed %s successfully' % normalizer)
        except Exception as e:
            logger.warn('Exception wile normalizing with %s' % normalizer, exc_info=e)
            self.proc.fail(e)
            self.save()
            return

    # update search
    self.proc.continue_with('archiving')
    self.save()
    try:
        self.archive(
            parser_backend,
            upload_hash=upload_hash,
            calc_hash=self.calc_hash,
            upload_id=self.upload.upload_id,
            mainfile=mainfile,
            upload_time=datetime.now())
        logger.debug('Archived successfully')
    except Exception as e:
        logger.error('Failed to archive', exc_info=e)
        self.proc.fail(e)
        self.save()
        return

    logger.debug('Completed processing')
    self.proc.success()
    self.save()
    self.upload.completed_calc()


class Calc(BaseDocument):
    upload = ReferenceField(Upload, required=True)
    parser = StringField(required=True)
    mainfile = StringField(required=True)
    calc_hash = StringField(required=True)
    tmp_mainfile = StringField(required=True)

    @property
    def archive_id(self):
        return '%s/%s' % (self.upload.upload_hash, self.calc_hash)

    def create(cls, upload: Upload, mainfile: str, parser: str, tmp_mainfile: str):
        self = Calc(
            upload=upload, mainfile=mainfile, parser=parser, tmp_mainfile=tmp_mainfile,
            calc_hash=utils.hash(mainfile))

    def process_calc(self):
        self.save()
        process_calc.s(self.upload.upload_id, self.mainfile).delay()


class Upload(BaseDocument):

    name = StringField()
    upload_hash = StringField()
    upload_time = DateTimeField()

    proc = EmbeddedDocumentField(Proc, required=True)

    presigned_url = StringField()

    calcs = IntField(default=0)
    processed_calcs = IntField(default=0)

    user = ReferenceField(User, required=True)

    meta = {
        'indexes': [
            'user'
        ]
    }

    @property
    def upload_id(self) -> str:
        return self.id.__str__()

    @property
    def is_stale(self) -> bool:
        create_time = self.proc.create_time

        return self.upload_time is None and (datetime.now() - create_time).days > 1

    def logger(self, **kwargs):
        return get_logger(__name__, upload_id=self.upload_id, cls=self.__class__.__name__, **kwargs)

    def create(cls, name: str=None):
        proc = Proc(tasks=[UPLOADING, EXTRACTING, PARSING, CLEANING])
        self = Upload(proc=proc, name=name)
        self.save()

        self.presigned_url = files.get_presigned_upload_url(self.upload_id)
        self.proc.continue_with(UPLOADING)

    def process_uploaded_file(self):
        upload_id = self.upload_id
        if upload_id is None:
            self.save()

        process_upload.s(upload_id).delay()

    def _check_for_cleaning(self):
        if self.calcs == self.processed_calcs:
            clean_upload(upload_id).delay()

    def add_calc(self, mainfile, parser, tmp_mainfile):
        calc = Calc.create(self, filenmae, parser, tmp_mainfile)
        calc.process_calc()
        self.calcs += 1

    def start_parsing(self):
        self.proc.continue_with(PARSING)
        self.save()
        self._check_for_cleaning()

    def completed_calc(self):
        self.processed_calcs += 1
        self._check_for_cleaning()

    def delete(self):
        logger = self.logger(action='delete')

        if not (self.proc.ready or self.stale or self.proc.current_task == 'UPLOADING'):
            raise NotAllowedDuringProcessing()

        with lnr(logger, 'Delete upload file'):
            try:
                files.Upload(self.upload_id).delete()
            except KeyError:
                if self.proc.current_task == UPLOADING:
                    logger.debug('Upload exist, but file does not exist. It was probably aborted and deleted.')
                else:
                    logger.debug('Upload exist, but uploaded file does not exist.')

        if self.upload_hash is not None:
            with lnr(logger, 'Deleting calcs'):
                Calc.delete_all(upload_id=self.upload_id)

        with lnr(logger, 'Deleting upload'):
            super().delete()

        return self

    @staticmethod
    def user_uploads(user: User) -> List['Upload']:
        """ Returns all uploads for the given user. Currently returns all uploads. """
        return [upload.update_proc() for upload in Upload.objects()]

    @property
    def json_dict(self) -> dict:
        """ A json serializable dictionary representation. """
        data = {
            'name': self.name,
            'upload_id': self.upload_id,
            'upload_hash': self.upload_hash,
            'presigned_url': files.external_objects_url(self.presigned_url),
            'upload_time': self.upload_time.isoformat() if self.upload_time is not None else None,
            'proc_time': self.proc_time.isoformat() if self.proc_time is not None else None,
            'stale': self.is_stale,
            'proc': self.proc.json_dict
        }
        return {key: value for key, value in data.items() if value is not None}
