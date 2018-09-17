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
This module is about maintaining the repository search index and providing all
data to the repository related parts of nomad.

We use *elasticsearch_dsl* to interface with elastic search. The class :class:`RepoCalc`
is an elasticsearch_dsl document that is used to represent repository index entries.

.. autoclass:: nomad.repo.RepoCalc
        :members:
"""

from typing import Dict, Any
import sys
from elasticsearch.exceptions import ConflictError, RequestError, ConnectionTimeout
from elasticsearch_dsl import Document as ElasticDocument, Search, Date, Keyword, Boolean, \
    connections
from datetime import datetime
import time

from nomad import config
from nomad.parsing import LocalBackend
from nomad.utils import get_logger

logger = get_logger(__name__)

# ensure elastic and mongo connections
if 'sphinx' not in sys.modules:
    client = connections.create_connection(hosts=[config.elastic.host])

key_mappings = {
    'basis_set_type': 'program_basis_set_type',
    'chemical_composition': 'chemical_composition_bulk_reduced'
}


class AlreadyExists(Exception): pass


class RepoCalc(ElasticDocument):
    """
    Elastic search document that represents a calculation. It is supposed to be a
    component of :class:`Calc`. Should only be created by its parent :class:`Calc`
    instance and only via the :func:`create_from_backend` factory method.
    """
    class Index:
        name = config.elastic.calc_index

    calc_hash = Keyword()
    mainfile = Keyword()
    upload_hash = Keyword()
    upload_id = Keyword()

    upload_time = Date()

    staging = Boolean()
    restricted = Boolean()
    user_id = Keyword()

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

    @property
    def archive_id(self) -> str:
        """ The unique id for this calculation. """
        return '%s/%s' % (self.upload_hash, self.calc_hash)

    @classmethod
    def create_from_backend(
            cls, backend: LocalBackend, additional: Dict[str, Any],
            upload_id: str, upload_hash: str, calc_hash: str,
            **kwargs) -> 'RepoCalc':
        """
        Create a new calculation instance in elastic search. The data from the given backend
        will be used. Additional meta-data can be given as *kwargs*. ``upload_id``,
        ``upload_hash``, and ``calc_hash`` are mandatory.

        Arguments:
            backend: The parsing/normalizing backend that contains the calculation data.
            additional: Additional arguments not stored in the backend. E.g. ``user_id``,
                ``staging``, ``restricted``
            upload_hash: The upload hash of the originating upload.
            upload_id: The upload id of the originating upload.
            calc_hash: The upload unique hash for this calculation.
            kwargs: Arguments are passed to elasticsearch index operation.

        Raises:
            AlreadyExists: If the calculation already exists in elastic search. We use
                the elastic document lock here. The elastic document is IDed via the
                ``archive_id``.
        """
        assert upload_hash is not None and calc_hash is not None and upload_id is not None
        additional.update(dict(upload_hash=upload_hash, calc_hash=calc_hash, upload_id=upload_id))

        # prepare the entry with all necessary properties from the backend
        calc = cls(meta=dict(id='%s/%s' % (upload_hash, calc_hash)))
        for property in cls._doc_type.mapping:
            property = key_mappings.get(property, property)

            if property in additional:
                value = additional[property]
            else:
                try:
                    value = backend.get_value(property, 0)
                except KeyError:
                    try:
                        program_name = backend.get_value('program_name', 0)
                    except KeyError:
                        program_name = 'unknown'
                    logger.warning(
                        'Missing property value', property=property, upload_id=upload_id,
                        upload_hash=upload_hash, calc_hash=calc_hash, code=program_name)
                    continue

            setattr(calc, property, value)

        # persist to elastic search
        try:
            # In practive es operation might fail due to timeout under heavy loads/
            # bad configuration. Retries with a small delay is a pragmatic solution.
            e_after_retries = None
            for _ in range(0, 2):
                try:
                    calc.save(op_type='create', **kwargs)
                    e_after_retries = None
                    break
                except ConnectionTimeout as e:
                    e_after_retries = e
                    time.sleep(1)
                except ConflictError as e:  # this should never happen, but happens
                    e_after_retries = e
                    time.sleep(1)
                else:
                    raise e
            if e_after_retries is not None:
                # if we had and exception and could not fix with retries, throw it
                raise e_after_retries  # pylint: disable=E0702
        except ConflictError:
            raise AlreadyExists('Calculation %s does already exist.' % (calc.archive_id))

        return calc

    @staticmethod
    def delete_upload(upload_id):
        """ Deletes all repo entries of the given upload. """
        RepoCalc.search().query('match', upload_id=upload_id).delete()

    @classmethod
    def unstage(cls, upload_id, staging=False):
        """ Update the staging property for all repo entries of the given upload. """
        cls.update_by_query(upload_id, {
            'inline': 'ctx._source.staging=%s' % ('true' if staging else 'false'),
            'lang': 'painless'
        })

    @classmethod
    def update_upload(cls, upload_id, **kwargs):
        """ Update all entries of given upload with keyword args. """
        for calc in RepoCalc.search().query('match', upload_id=upload_id):
            calc.update(**kwargs)

    @classmethod
    def update_by_query(cls, upload_id, script):
        """ Update all entries of a given upload via elastic script. """
        index = cls._default_index()
        doc_type = cls._doc_type.name
        conn = cls._get_connection()
        body = {
            'script': script,
            'query': {
                'match': {
                    'upload_id': upload_id
                }
            }
        }
        response = conn.update_by_query(index, doc_type=[doc_type], body=body)

    @staticmethod
    def es_search(body):
        """ Perform an elasticsearch and not elasticsearch_dsl search on the Calc index. """
        return client.search(index=config.elastic.calc_index, body=body)

    @staticmethod
    def upload_exists(upload_hash):
        """ Returns true if there are already calcs from the given upload. """
        search = Search(using=client, index=config.elastic.calc_index) \
            .query('match', upload_hash=upload_hash) \
            .execute()

        return len(search) > 0

    @property
    def json_dict(self):
        """ A json serializable dictionary representation. """
        data = self.to_dict()

        upload_time = data.get('upload_time', None)
        if upload_time is not None and isinstance(upload_time, datetime):
            data['upload_time'] = data['upload_time'].isoformat()

        data['archive_id'] = self.archive_id

        return {key: value for key, value in data.items() if value is not None}


if 'sphinx' not in sys.modules:
    try:
        RepoCalc.init()
    except RequestError as e:
        if e.status_code == 400 and 'resource_already_exists_exception' in e.error:
            pass  # happens if two services try this at the same time
        else:
            raise e
