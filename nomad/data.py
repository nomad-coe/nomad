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
This module comprises a set of persistent document classes that hold all user related
data. These are information about users, their uploads and datasets, the associated
calculations, and files.

.. figure:: data.png
   :alt: nomad xt data model

   The main nomad xt data classes and their associations

The blue properties are exposed via our REST API.

For detailed information on :class:`nomad.processing.UploadProc` and
:class:`nomad.processing.CalcProc` refer to :py:mod:`nomad.processing`.

.. autoclass:: Calc
    :members:
.. autoclass:: Upload
    :members:
.. autoclass:: DataSet
.. autoclass:: User

"""

from typing import List
import sys
from datetime import datetime
import elasticsearch.exceptions
from elasticsearch_dsl import Document as ElasticDocument, Date, Keyword, Search, connections
from mongoengine import \
    Document, EmailField, StringField, BooleanField, DateTimeField, \
    ListField, DictField, ReferenceField, connect
import mongoengine.errors

from nomad import config, files
from nomad.processing import UploadProc, CalcProc
from nomad.parsing import LocalBackend
from nomad.utils import get_logger, lnr

logger = get_logger(__name__)

# ensure elastic and mongo connections
if 'sphinx' not in sys.modules:
    client = connections.create_connection(hosts=[config.elastic.host])
    connect(db=config.mongo.users_db, host=config.mongo.host)

key_mappings = {
    'basis_set_type': 'program_basis_set_type',
    'chemical_composition': 'chemical_composition_bulk_reduced'
}


class AlreadyExists(Exception): pass


class InvalidId(Exception): pass


class NotAllowedDuringProcessing(Exception): pass


class Calc(ElasticDocument):
    """
    Instances of this class represent calculations. This class manages the elastic
    search index entry, files, and archive for the respective calculation.
    Instances should be created directly, but via the static :func:`create_from_backend`
    function.

    The attribute list, does not include the various repository properties generated
    while parsing, including ``program_name``, ``program_version``, etc.

    Attributes:
        calc_hash: The hash that identified the calculation within an upload
        upload_hash: The hash of the upload
        upload_id: The id of the upload used to create this calculation
        mainfile: The mainfile (including path in upload) that was used to create this calc
    """
    calc_hash = Keyword()

    upload_time = Date()
    upload_id = Keyword()
    upload_hash = Keyword()
    mainfile = Keyword()

    program_name = Keyword()
    program_version = Keyword()

    chemical_composition = Keyword()
    basis_set_type = Keyword()
    atom_species = Keyword()
    system_type = Keyword()
    crystal_system = Keyword()
    space_group_number = Keyword()
    configuration_raw_gid = Keyword()
    XC_functional_name = Keyword()

    class Index:
        name = config.elastic.calc_index

    @property
    def archive_id(self) -> str:
        """ The unique id for this calculation. """
        return '%s/%s' % (self.upload_hash, self.calc_hash)

    def delete(self):
        """
        Delete this calculation and all associated data. This includes all files,
        the archive, and this search index entry.
        """
        # delete the archive
        files.delete_archive(self.archive_id)

        # delete the search index entry
        super().delete()

    @staticmethod
    def es_search(body):
        """ Perform an elasticsearch and not elasticsearch_dsl search on the Calc index. """
        return client.search(index=config.elastic.calc_index, body=body)

    @staticmethod
    def delete_all(**kwargs):
        for calc in Calc.search().query('match', **kwargs).execute():
            calc.delete()

    @staticmethod
    def create_from_backend(
            backend: LocalBackend, upload_id: str, upload_hash: str, calc_hash: str, **kwargs) \
            -> 'Calc':
        """
        Create a new calculation instance. The data from the given backend
        will be used. Additional meta-data can be given as *kwargs*. ``upload_id``,
        ``upload_hash``, and ``calc_hash`` are mandatory.
        This will create a elastic search entry and store the backend data to the
        archive.

        Arguments:
            backend: The parsing/normalizing backend that contains the calculation data.
            upload_hash: The upload hash of the originating upload.
            upload_id: The upload id of the originating upload.
            calc_hash: The upload unique hash for this calculation.
            kwargs: Additional arguments not stored in the backend.

        Raises:
            AlreadyExists: If the calculation already exists in elastic search. We use
                the elastic document lock here. The elastic document is ided via the
                ``archive_id``.
        """
        assert upload_hash is not None and calc_hash is not None and upload_id is not None
        kwargs.update(dict(upload_hash=upload_hash, calc_hash=calc_hash, upload_id=upload_id))

        calc = Calc(meta=dict(id='%s/%s' % (upload_hash, calc_hash)))

        for property in Calc._doc_type.mapping:
            property = key_mappings.get(property, property)

            if property in kwargs:
                value = kwargs[property]
            else:
                try:
                    value = backend.get_value(property, 0)
                except KeyError:
                    logger.warning(
                        'Missing property value', property=property, upload_id=upload_id,
                        upload_hash=upload_hash, calc_hash=calc_hash)
                    continue

            setattr(calc, property, value)

        # persist to elastic search
        try:
            calc.save(op_type='create')
        except Exception:
            raise AlreadyExists('Calculation %s does already exist.' % (calc.archive_id))

        # persist the archive
        with files.write_archive_json(calc.archive_id) as out:
            backend.write_json(out, pretty=True)

        return calc

    @property
    def json_dict(self):
        """ A json serializable dictionary representation. """
        data = self.to_dict()

        upload_time = data.get('upload_time', None)
        if upload_time is not None and isinstance(upload_time, datetime):
            data['upload_time'] = data['upload_time'].isoformat()

        data['archive_id'] = self.archive_id

        return {key: value for key, value in data.items() if value is not None}

    @staticmethod
    def upload_exists(upload_hash):
        """ Returns true if there are already calcs from the given upload. """
        search = Search(using=client, index=config.elastic.calc_index) \
            .query('match', upload_hash=upload_hash) \
            .execute()

        return len(search) > 0


if 'sphinx' not in sys.modules:
    try:
        Calc.init()
    except elasticsearch.exceptions.RequestError as e:
        if e.status_code == 400 and 'resource_already_exists_exception' in e.error:
            pass  # happens if two services try this at the same time
        else:
            raise e


class User(Document):
    """ Represents users in the database. """
    email = EmailField(primary=True)
    name = StringField()


class Upload(Document):
    """
    Represents uploads in the databases. Provides persistence access to the files storage,
    and processing state.

    Attributes:
        file_name: Optional user provided upload name
        upload_id: The upload id. Generated by the database.
        in_staging: True if the upload is still in staging and can be edited by the uploader.
        is_private: True if the upload and its derivitaves are only visible to the uploader.
        proc: The :class:`nomad.processing.UploadProc` that holds the processing state.
        created_time: The timestamp this upload was created.
        upload_time: The timestamp when the system realised the upload.
        proc_time: The timestamp when the processing realised finished by the system.
    """

    name = StringField(default=None)

    in_staging = BooleanField(default=True)
    is_private = BooleanField(default=False)

    presigned_url = StringField()
    upload_time = DateTimeField()
    create_time = DateTimeField()

    proc_time = DateTimeField()
    proc = DictField()

    user = ReferenceField(User, required=True)

    meta = {
        'indexes': [
            'proc.upload_hash',
            'user'
        ]
    }

    @staticmethod
    def user_uploads(user: User) -> List['Upload']:
        """ Returns all uploads for the given user. Currently returns all uploads. """
        return [upload.update_proc() for upload in Upload.objects()]

    @staticmethod
    def get(upload_id: str) -> 'Upload':
        try:
            upload = Upload.objects(id=upload_id).first()
        except mongoengine.errors.ValidationError:
            raise InvalidId('Invalid upload id')

        if upload is None:
            raise KeyError('Upload does not exist')

        upload.update_proc()
        return upload

    def logger(self, **kwargs):
        return get_logger(__name__, upload_id=self.upload_id, cls='Upload', **kwargs)

    def delete(self) -> 'Upload':
        logger = self.logger(action='delete')

        self.update_proc()
        if not (self.is_ready or self.is_stale or self._proc.current_task_name == 'uploading'):
            raise NotAllowedDuringProcessing()

        with lnr(logger, 'Delete upload file'):
            try:
                files.Upload(self.upload_id).delete()
            except KeyError:
                if self._proc.current_task_name == 'uploading':
                    logger.debug('Upload exist, but file does not exist. It was probably aborted and deleted.')
                else:
                    logger.debug('Upload exist, but uploaded file does not exist.')

        if self._proc.upload_hash is not None:
            with lnr(logger, 'Deleting calcs'):
                Calc.delete_all(upload_id=self.upload_id)

        with lnr(logger, 'Deleting upload'):
            super().delete()

        return self

    @staticmethod
    def create(user: User, name: str=None) -> 'Upload':
        """
        Creates a new upload for the given user, a user given name is optional.
        It will populate the record with a signed url and pending :class:`UploadProc`.
        The upload will be already saved to the database.
        """
        upload = Upload(user=user, name=name)
        upload.save()

        upload.presigned_url = files.get_presigned_upload_url(upload.upload_id)
        upload.create_time = datetime.now()
        upload.proc = UploadProc(upload.upload_id)
        upload.save()

        upload.update_proc()

        return upload

    @property
    def upload_id(self) -> str:
        return self.id.__str__()

    @property
    def is_stale(self) -> bool:
        proc = self._proc
        if proc.current_task_name == proc.task_names[0] and self.upload_time is None:
            return (datetime.now() - self.create_time).days > 1
        else:
            return False

    @property
    def is_ready(self) -> bool:
        return self._proc.status in ['SUCCESS', 'FAILURE']

    @property
    def _proc(self):
        """ Cast the internal mongo dict to an actual :class:`UploadProc` instance. """
        # keep the instance cached
        if '__proc' not in self.__dict__:
            self.__dict__['__proc'] = UploadProc(**self.proc)

        return self.__dict__['__proc']

    def update_proc(self) -> 'Upload':
        """ Updates this instance with information from the celery results backend. """
        if self._proc.update_from_backend():
            self.proc = self._proc
            self.save()

        return self

    @property
    def json_dict(self) -> dict:
        """ A json serializable dictionary representation. """
        data = {
            'name': self.name,
            'upload_id': self.upload_id,
            'presigned_url': files.external_objects_url(self.presigned_url),
            'create_time': self.create_time.isoformat() if self.create_time is not None else None,
            'upload_time': self.upload_time.isoformat() if self.upload_time is not None else None,
            'proc_time': self.proc_time.isoformat() if self.proc_time is not None else None,
            'is_stale': self.is_stale,
            'is_ready': self.is_ready,
            'proc': self._proc
        }
        return {key: value for key, value in data.items() if value is not None}


class DataSet(Document):
    name = StringField()
    description = StringField()
    doi = StringField()

    user = ReferenceField(User)
    calcs = ListField(StringField)

    meta = {
        'indexes': [
            'user',
            'doi',
            'calcs'
        ]
    }

# provid a fake user for testing
me = None
if 'sphinx' not in sys.modules:
    me = User.objects(email='me@gmail.com').first()
    if me is None:
        me = User(email='me@gmail.com', name='Me Meyer')
        me.save()
