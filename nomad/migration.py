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
This module contains functions to read data from NOMAD coe, external sources,
other/older nomad@FAIRDI instances to mass upload it to a new nomad@FAIRDI instance.
"""

from nomad import infrastructure
from nomad.coe_repo import User


class NomadCOEMigration:
    """
    Drives a migration from the NOMAD coe repository db to nomad@FAIRDI.

    Arguments:
        **kwargs: Are used to configure the source repository db analog to
        :data:`config.repository_db`

    Attributes:
        source_repo_db: SQLAlchemy session for the source NOMAD coe repository db
    """
    def __init__(self, **kwargs):
        self.source_connection, self.source_repo_db = infrastructure.sqlalchemy_repository_db(
            readonly=True, **kwargs)

    def copy_users(self):
        """ Copy all users, keeping their ids, within a single transaction. """
        repo_db = infrastructure.repository_db

        repo_db.begin()
        for source_user in self.source_repo_db.query(User).all():
            self.source_repo_db.expunge(source_user)  # removes user from the source session
            repo_db.merge(source_user)

        repo_db.commit()
