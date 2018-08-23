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

from typing import List, Any, Union, cast
from celery.result import AsyncResult, result_from_tuple
import itertools

from nomad.normalizing import normalizers
from nomad.utils import DataObject
from nomad.processing.app import app


class ProcPipeline(DataObject):
    """
    Arguments:
        task_names: A list of task names in pipeline order.

    Attributes:
        current_task_name: Name of the currently running task.
        status: Aggregated status for the whole process. Celery status names as convention.
        errors: A list of potential error that caused failure.
    """
    def __init__(self, task_names: List[str], *args, **kwargs) -> None:
        super().__init__(*args)
        self.task_names: List[str] = task_names
        self.current_task_name: str = None
        self.status: str = 'PENDING'
        self.errors: List[str] = []

        self.update(kwargs)

    def continue_with(self, task_name: str) -> bool:
        """ Upadtes itself with information about the new current task. """

        assert self.status != 'SUCCESS', 'Cannot continue on completed workflow.'

        if self.status == 'FAILURE':
            return False
        else:
            self.status = 'PROGRESS'
            self.current_task_name = task_name
            return True

    def success(self) -> None:
        self.status = 'SUCCESS'

    @property
    def is_started(self) -> bool:
        """ True, if the task is started. """
        return self.status is not 'PENDING'

    def fail(self, e: Union[List[str], List[Exception], Exception, str]):
        """ Allows tasks to mark this processing as failed. All following task will do nothing. """
        raw_errors: Union[List[str], List[Exception]] = None
        if isinstance(e, list):
            raw_errors = e
        else:
            raw_errors = cast(Union[List[str], List[Exception]], [e])

        for error in raw_errors:
            if isinstance(error, str):
                self.errors.append(error)
            elif isinstance(error, Exception):
                self.errors.append(error.__str__())
            else:
                assert False, 'Unknown error'

        self.status = 'FAILURE'


class CalcProc(ProcPipeline):
    """
    Used to represent the state of an calc processing. It is used to provide
    information to the user (via api, users) and act as a state object within the
    more complex calc processing task. Keep in mind that this task might become several
    celery tasks in the future.

    Arguments:
        upload_hash: The hash that identifies the upload in the archive.
        mainfile: The path to the mainfile in the upload.
        parser_name: The name of the parser to use/used.
        tmp_mainfile: The full path to the mainfile in the local fs.

    Attributes:
        calc_hash: The mainfile hash that identifies the calc in the archive.
        archive_id: The id that identifies the archive via `upload_hash/calc_hash`.
        celery_task_id: The celery task id for the calc parse celery task.
    """
    def __init__(self, upload_hash, mainfile, parser_name, tmp_mainfile, *args, **kwargs):
        task_names = [
            [parser_name],
            [n.__name__ for n in normalizers],
            ['indexing', 'storage']
        ]

        super().__init__(task_names=list(itertools.chain(*task_names)), *args)

        self.upload_hash = upload_hash
        self.mainfile = mainfile
        self.parser_name = parser_name
        self.tmp_mainfile = tmp_mainfile

        self.calc_hash = hash(mainfile)
        self.archive_id = '%s/%s' % (self.upload_hash, self.calc_hash)

        self.celery_task_id: str = None

        self.update(kwargs)

    def update_from_backend(self):
        assert self.celery_task_id is not None

        celery_task_result = AsyncResult(self.celery_task_id, app=app)
        if celery_task_result.ready():
            self.update(celery_task_result.result)
        else:
            info = celery_task_result.info
            if info is not None:
                self.update(info)


class UploadProc(ProcPipeline):
    """
    Used to represent the state of an upload processing. It is used to provide
    information to the user (via api, users) and act as a state object that is passed
    from celery task to celery task.

    It is serializable (JSON, pickle). Iternaly stores
    :class:`~celery.results.AsyncResults` instance in serialized *tuple* form to
    keep connected to the results backend.

    Warning:
        You have to call :func:`forget` eventually to free all resources and the celery
        results backend.

        Anyhow, results will be deleted after 1 day, depending on `configuration
        <http://docs.celeryproject.org/en/latest/userguide/configuration.html#result-expires>`_.

    Arguments:
        upload_id: The id of the uploaded file in the object storage,
                   see also :mod:`nomad.files`.
        task_name: The names of all task in pipeline order.

    Attributes:
        upload_hash: The hash of the uploaded file. E.g., used for archive/repo ids.
        calc_procs: The state data for all child calc processings.
        celery_task_ids: Serialized form of the celery async_results tree for the processing.
    """
    def __init__(self, upload_id: str, task_names: List[str], *args, **kwargs) -> None:
        assert upload_id is not None
        super().__init__(task_names, *args)

        self.upload_id = upload_id
        self.upload_hash: str = None

        self.calc_procs: List[CalcProc] = []

        self.celery_task_ids: Any = None

        self.update(kwargs)

    @property
    def _celery_task_result(self) -> AsyncResult:
        """
        The celery async_result in its regular usable, but not serializable form.

        We use the tuple form to allow serialization (i.e. storage). Keep in mind
        that the sheer `task_id` is not enough, because it does not contain
        the parent tasks, i.e. result tree.
        See `third comment <https://github.com/celery/celery/issues/1328>`_
        for details.
        """
        assert self.celery_task_ids is not None

        return result_from_tuple(self.celery_task_ids, app=app)

    def update_from_backend(self) -> 'UploadProc':
        """ Consults the result backend and updates itself with the available results. """
        assert self.is_started, 'Run is not yet started.'

        celery_task_result = self._celery_task_result

        while celery_task_result is not None:
            if celery_task_result.ready():
                self.update(celery_task_result.result)
                break
            else:
                celery_task_result = celery_task_result.parent

        if self.calc_procs is not None:
            for calc_proc in self.calc_procs:
                calc_proc.update_from_backend()

        return self

    def forget(self) -> None:
        """ Forget the results of a completed run; free all resources in the results backend. """
        # TODO, this is not forgetting the parse task in the parse_all header, right?
        assert self.ready(), 'Run is not completed.'

        celery_task_result = self._celery_task_result
        while celery_task_result is not None:
            celery_task_result.forget()
            celery_task_result = celery_task_result.parent

    def ready(self) -> bool:
        """ Returns: True if the task has been executed. """
        assert self.is_started, 'Run is not yet started.'

        return self._celery_task_result.ready()

    def get(self, *args, **kwargs) -> 'UploadProc':
        """
        Blocks until the processing has finished. Forwards args, kwargs to
        *celery.result.get* for timeouts, etc.

        Returns: An upadted instance of itself with all the results.
        """
        assert self.is_started, 'Run is not yet started.'

        self._celery_task_result.get(*args, **kwargs)
        self.update_from_backend()
        return self
