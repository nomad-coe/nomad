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

"""
Processing comprises everything that is necessary to take an uploaded user file,
processes it, and store all necessary data for *repository*, *archive*, and potential
future services (e.g. *encyclopedia*).

Processing is build on top of *celery* (http://www.celeryproject.org/).
Celery provides a task-based programming model for distributed computing. It uses
a broker, e.g. a distributed task queue like *RabbitMQ* (), to distribute task requests,
and a result backend, e.g. a *Redis* database (), to access (intermediate) task results.
This combination allows us to easily distribute processing work while having
the processing state, i.e. (intermediate) results, always available.

This module is structures into our *celery app* (``app.py``), the task definitions
(``tasks.py``), classes that represent state for processing *uploads* and their
*calculations* (``state.py``), and the *handler* service that initiates processing
based on file storage notifications (``handler.py``, ``handlerdaemon.py``).

This module does not contain the functions to do the actual work. Those are encapsulated
in :py:mod:`nomad.files`, :py:mod:`nomad.search`, :py:mod:`nomad.users`,
:py:mod:`nomad.parsing`, and :py:mod:`nomad.normalizing`.

Processing app
--------------

Refer to http://www.celeryproject.org/ to learn about celery apps and workers. The
nomad celery app uses a *RabbitMQ* broker and *Redis* result backend. It uses *pickle*
for serialization of arguments and results (usually instances of the :class:`Proc` state
classes).

Processing tasks
----------------

There are *upload* processing tasks (extracting, parse_all, cleanup) and a calculation
processing task (parse). The upload processing tasks are ment to be executed in sequence.

.. figure:: proc.png
   :alt: nomad xt processing workflow

   This is the basic workflow of a nomad xt upload processing.


.. autotask:: nomad.processing.tasks.extracting_task
.. autotask:: nomad.processing.tasks.cleanup_task
.. autotask:: nomad.processing.tasks.parse_all_task
.. autotask:: nomad.processing.tasks.parse_task


Represent processing state
--------------------------

.. autoclass:: Proc
    :members:
.. autoclass:: UploadProc
    :members:
.. autoclass:: CalcProc
    :members:

Initiate processing
-------------------

.. autofunction:: start_processing
.. autofunction:: handle_uploads
.. autofunction:: handle_uploads_thread

"""

from nomad.processing.base import app, InvalidId, ProcNotRegistered
from nomad.processing.data import Upload, Calc, NotAllowedDuringProcessing
from nomad.processing.handler import handle_uploads, handle_uploads_thread
