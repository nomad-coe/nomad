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
Interface to the NOMAD-coe repository postgres database. This implementation is based on
SQLAlchemy. There are model classes that represent entries in the *users* and *session*
tables. All DB entities are implemented as classes, but most are hidden and data
can be accessed via the various relations with :class:`User`, :class:`Calc`, :class:`Upload`.

To load an entity from the database use :data:`nomad.infrastructure.repository_db`
(the SQLAlchemy session), e.g.:

.. code-block:: python

    repository_db.Query(coe_repo.Calc).filter_by(upload_id=some_id)

.. autoclass:: User
    :members:
    :undoc-members:

.. autofunction:: ensure_test_user
.. autodata:: admin_user
.. autoexception:: LoginException

.. autoclass:: UploadMetaData
    :members:
.. autoclass:: Upload
    :members:
    :undoc-members:
.. autoclass:: Calc
    :members:
    :undoc-members:
.. autoclass:: DataSet
    :members:
    :undoc-members:
"""

from .user import User, ensure_test_user, admin_user, LoginException
from .calc import Calc, DataSet
from .upload import UploadMetaData, Upload
