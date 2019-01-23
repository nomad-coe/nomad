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

import os.path
import json
from mongoengine import Document, IntField, StringField, DictField
import itertools

from nomad import utils
from nomad.coe_repo import User, Calc
from nomad.coe_repo.base import CalcMetaData


def read_metadata(calc: Calc) -> dict:
    metadata = dict(
        _pid=calc.pid,
        _uploader=calc.uploader.user_id,
        _upload_time=calc.calc_metadata.added,
        _datasets=list(
            dict(id=ds.id, dois=ds.dois, name=ds.name)
            for ds in calc.all_datasets),
        with_embargo=calc.with_embargo,
        comment=calc.comment,
        references=calc.references,
        coauthors=list(user.user_id for user in calc.coauthors),
        shared_with=list(user.user_id for user in calc.shared_with))
    return {
        key: value for key, value in metadata.items()
        if value is not None and value != []
    }


class SourceCalc(Document):
    """
    Mongo document used as a calculation, upload, and metadata db and index
    build from a given source db.
    """
    pid = IntField(primary_key=True)
    mainfile = StringField()
    upload = StringField()
    metadata = DictField()

    extracted_prefix = '$EXTRACTED/'
    sites = ['/data/nomad/extracted/', '/nomad/repository/extracted/']
    prefixes = itertools.chain([extracted_prefix], sites)

    @staticmethod
    def index(source, drop: bool = False):
        if drop:
            SourceCalc.drop_collection()

        for metadata in source.query(CalcMetaData).filter(CalcMetaData.filenames.isnot(None)).yield_per(1000):
            filenames = json.loads(metadata.filenames.decode('utf-8'))
            filename = filenames[0]
            for prefix in SourceCalc.prefixes:
                filename = filename.replace(prefix, '')
            segments = [file.strip('\\') for file in filename.split('/')]

            upload = segments[0]
            mainfile = os.path.join(*segments[1:])
            pid = metadata.calc.pid

            source_calc = SourceCalc(mainfile=mainfile, pid=pid, upload=upload)
            source_calc.metadata = read_metadata(metadata.calc)
            source_calc.save()


class NomadCOEMigration:
    """
    Drives a migration from the NOMAD coe repository db to nomad@FAIRDI. It is assumed
    that this class is never used on the worker or api service. It assumes the
    default coe repo connection as a connection to the source repository db.

    Attributes:
        source: SQLAlchemy session for the source NOMAD coe repository db.
    """
    prefixes = ['$extracted', '/data/nomad/extracted']
    sites = [
        '/nomad/repository/nomad/extracted',
        '/nomad/']

    def __init__(self):
        self.logger = utils.get_logger(__name__)
        from nomad.infrastructure import repository_db
        self.source = repository_db

    def copy_users(self, target_db):
        """ Copy all users, keeping their ids, within a single transaction. """
        target_db.begin()
        for source_user in self.source.query(User).all():
            self.source.expunge(source_user)  # removes user from the source session
            target_db.merge(source_user)

        target_db.commit()

    def migrate(self, upload: str):
        """ Migrate the given upload. """
        pass

    def index(self, *args, **kwargs):
        SourceCalc.index(self.source, *args, **kwargs)
