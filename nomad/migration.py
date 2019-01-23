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

from typing import Generator, Tuple, List
import os.path
import json
from mongoengine import Document, IntField, StringField, DictField

from nomad import utils
from nomad.coe_repo import User, Calc, DataSet


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
    prefixes = [extracted_prefix] + sites

    _dataset_cache: dict = {}

    @staticmethod
    def _read_metadata(calc: Calc) -> dict:
        datasets: List[DataSet] = []
        for parent in calc.parents:
            parents = SourceCalc._dataset_cache.get(parent, None)
            if parents is None:
                parents = parent.all_datasets
                SourceCalc._dataset_cache[parent] = parents
            datasets.append(DataSet(parent))
            datasets.extend(parents)

        metadata = dict(
            _pid=calc.pid,
            _uploader=calc.uploader.user_id,
            _upload_time=calc.calc_metadata.added,
            _datasets=list(
                dict(id=ds.id, dois=ds.dois, name=ds.name)
                for ds in datasets),
            with_embargo=calc.with_embargo,
            comment=calc.comment,
            references=calc.references,
            coauthors=list(user.user_id for user in calc.coauthors),
            shared_with=list(user.user_id for user in calc.shared_with)
        )
        return {
            key: value for key, value in metadata.items()
            if value is not None and value != []
        }

    @staticmethod
    def index(source, drop: bool = False, with_metadata: bool = True, per_query: int = 100) \
            -> Generator[Tuple['SourceCalc', int], None, None]:
        """
        Creates a collection of :class:`SourceCalc` documents that represent source repo
        db entries. Each document relates a calc's (pid, mainfile, upload). Where
        upload is the 'id'/prefix of an upload directory or upload file in the source repo's
        filesystem.

        Arguments:
            source: The source db sql alchemy session
            drop: True to create a new collection, update the existing otherwise, default is False.
            with_metadata: True to also grab all user metadata and store it, default is True.

        Returns:
            yields tuples (:class:`SourceCalc`, #calcs_total)
        """
        if drop:
            SourceCalc.drop_collection()

        last_source_calc = SourceCalc.objects().order_by('pid').first()
        start_pid = last_source_calc.pid if last_source_calc is not None else 0
        source_query = source.query(Calc)
        total = source_query.count()

        while True:
            calcs = source_query.filter(Calc.coe_calc_id > start_pid).order_by(Calc.coe_calc_id).limit(per_query)
            source_calcs = []
            for calc in calcs:
                if calc.calc_metadata is None or calc.calc_metadata.filenames is None:
                    yield None, total
                    continue  # dataset case

                filenames = json.loads(calc.calc_metadata.filenames.decode('utf-8'))
                filename = filenames[0]
                for prefix in SourceCalc.prefixes:
                    filename = filename.replace(prefix, '')
                segments = [file.strip('\\') for file in filename.split('/')]

                source_calc = SourceCalc(pid=calc.pid)
                source_calc.upload = segments[0]
                source_calc.mainfile = os.path.join(*segments[1:])
                if with_metadata:
                    source_calc.metadata = SourceCalc._read_metadata(calc)
                source_calcs.append(source_calc)
                start_pid = source_calc.pid

                yield source_calc, total

            if len(source_calcs) == 0:
                break
            else:
                SourceCalc.objects.insert(source_calcs)


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
        return SourceCalc.index(self.source, *args, **kwargs)
