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
This module represents calculations in elastic search.
"""

from elasticsearch.exceptions import ConflictError, ConnectionTimeout
from datetime import datetime
import time
from elasticsearch_dsl import Document, InnerDoc, Keyword, Text, Date, \
    Nested

from nomad import config, datamodel, infrastructure, datamodel


class AlreadyExists(Exception): pass


class User(InnerDoc):
    def __init__(self, user):
        super().__init__(
            id=user.user_id,
            name='%s %s' % (user.first_name, user.last_name),
            name_keyword='%s %s' % (user.first_name, user.last_name))

    id = Keyword()
    name = Text()
    name_keyword = Keyword()


class Dataset(InnerDoc):
    def __init__(self, dataset):
        super().__init__(
            id=dataset.id,
            doi=dataset.doi.value if dataset.doi is not None else None,
            name=dataset.name)

    id = Keyword()
    doi = Keyword()
    name = Keyword()


class Entry(Document):
    class Index:
        name = config.elastic.index_name

    upload_id = Keyword()
    upload_time = Date(format='epoch_millis')
    calc_id = Keyword()
    calc_hash = Keyword()
    pid = Keyword()
    mainfile = Keyword()
    files = Keyword()
    uploader = Nested(User)

    with_embargo = Keyword()
    coauthors = Nested(User)
    shared_with = Nested(User)
    comment = Text()
    references = Keyword()
    datasets = Nested(Dataset)

    formula = Keyword()
    atoms = Keyword()
    basis_set = Keyword()
    xc_functional = Keyword()
    system = Keyword()
    crystal_system = Keyword()
    spacegroup = Keyword()
    code_name = Keyword()
    code_version = Keyword()

    @classmethod
    def from_calc_with_metadata(cls, source: datamodel.CalcWithMetadata) -> 'Entry':
        return Entry(
            meta=dict(id=source.calc_id),
            upload_id=source.upload_id,
            upload_time=source.upload_time,
            calc_id=source.calc_id,
            calc_hash=source.calc_hash,
            pid=str(source.pid),
            mainfile=source.mainfile,
            files=source.files,
            uploader=User(source.uploader) if source.uploader is not None else None,

            with_embargo=source.with_embargo,
            coauthors=[User(user) for user in source.coauthors],
            shared_with=[User(user) for user in source.shared_with],
            comment=source.comment,
            references=[ref.value for ref in source.references],
            datasets=[Dataset(ds) for ds in source.datasets],

            formula=source.formula,
            atoms=list(set(source.atoms)),
            basis_set=source.basis_set,
            xc_functional=source.xc_functional,
            system=source.system,
            crystal_system=source.crystal_system,
            spacegroup=source.spacegroup,
            code_name=source.code_name,
            code_version=source.code_version)

    def persist(self, **kwargs):
        """
            Persist this entry to elastic search. Kwargs are passed to elastic search.

            Raises:
                AlreadyExists: If the calculation already exists in elastic search. We use
                    the elastic document lock here. The elastic document is IDed via the
                    ``calc_id``.
        """
        try:
            # In practive es operation might fail due to timeout under heavy loads/
            # bad configuration. Retries with a small delay is a pragmatic solution.
            e_after_retries = None
            for _ in range(0, 2):
                try:
                    self.save(op_type='create', **kwargs)
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
            raise AlreadyExists('Calculation %s/%s does already exist.' % (self.upload_id, self.calc_id))

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
        conn.update_by_query(index, doc_type=[doc_type], body=body)

    @staticmethod
    def es_search(body):
        """ Perform an elasticsearch and not elasticsearch_dsl search on the Calc index. """
        return infrastructure.elastic_client.search(index=config.elastic.index_name, body=body)

    @property
    def json_dict(self):
        """ A json serializable dictionary representation. """
        data = self.to_dict()

        upload_time = data.get('upload_time', None)
        if upload_time is not None and isinstance(upload_time, datetime):
            data['upload_time'] = data['upload_time'].isoformat()

        return {key: value for key, value in data.items() if value is not None}


# Entry.register_mapping(datamodel.CalcWithMetadata, Entry.from_calc_with_metadata)
