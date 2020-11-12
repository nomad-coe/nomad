#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

'''
Processing comprises everything that is necessary to take an uploaded user file,
processes it, and store all necessary data for *repository*, *archive*, and potential
future services (e.g. *encyclopedia*).

Processing is build on top of *celery* (http://www.celeryproject.org/) and
*mongodb* (http://www.mongodb.org).
Celery provides a task-based programming model for distributed computing. It uses
a broker, e.g. a distributed task queue like *RabbitMQ* to distribute tasks. We
use mongodb to store the current state of processing in :class:`Upload` and
:class:`Calculation` documents. This combination allows us to easily distribute
processing work while having the processing state, i.e. (intermediate) results,
always available.

This module is structured into our *celery app* and abstract process base class
:class:`Proc` (``base.py``), and the concrete processing classes
:class:`Upload` and :class:`Calc` (``data.py``).

This module does not contain the functions to do the actual work. Those are encapsulated
in :py:mod:`nomad.files`, :py:mod:`nomad.repo`, :py:mod:`nomad.users`,
:py:mod:`nomad.parsing`, and :py:mod:`nomad.normalizing`.

Refer to http://www.celeryproject.org/ to learn about celery apps and workers. The
nomad celery app uses a *RabbitMQ* broker. We use celery to distribute processing load
in a cluster.

We use an abstract processing base class and document (:class:`Proc`) that provides all
necessary functions to execute a process as a series of potentially distributed steps. In
addition the processing state is persisted in mongodb using *mongoengine*. Instead of
exchanging serialized state between celery tasks, we use the mongodb documents to
exchange data. Therefore, the mongodb always contains the latest processing state.
We also don't have to deal with celery result backends and synchronizing with them.

.. autoclass:: nomad.processing.base.Proc

There are two concrete processes :class:`Upload` and :class: `Calc`. Instances of both
classes do represent the processing state, as well as the respective entity.

.. autoclass:: nomad.processing.data.Upload
    :members:
.. autoclass:: nomad.processing.data.Calc
    :members:
'''

from nomad.processing.base import app, InvalidId, ProcNotRegistered, SUCCESS, FAILURE, \
    RUNNING, PENDING, PROCESS_COMPLETED, PROCESS_RUNNING, ProcessAlreadyRunning
from nomad.processing.data import Upload, Calc
