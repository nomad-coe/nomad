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

from elasticsearch_dsl import Document, InnerDoc, Keyword, Text, Date, \
    Object, Boolean, Search

from nomad import config, datamodel, infrastructure, datamodel, coe_repo


class AlreadyExists(Exception): pass


class User(InnerDoc):

    @classmethod
    def from_user_popo(cls, user):
        self = cls(user_id=user.id)

        if 'first_name' not in user:
            user = coe_repo.User.from_user_id(user.id).to_popo()

        self.name = '%s %s' % (user['first_name'], user['last_name'])
        self.name_keyword = '%s %s' % (user['first_name'], user['last_name'])

        return self

    user_id = Keyword()
    name = Text()
    name_keyword = Keyword()


class Dataset(InnerDoc):

    @classmethod
    def from_dataset_popo(cls, dataset):
        return cls(
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
    upload_time = Date()
    calc_id = Keyword()
    calc_hash = Keyword()
    pid = Keyword()
    mainfile = Keyword()
    files = Keyword(multi=True)
    uploader = Object(User)

    with_embargo = Boolean()
    published = Boolean()

    coauthors = Object(User)
    shared_with = Object(User)
    comment = Text()
    references = Keyword()
    datasets = Object(Dataset)

    formula = Keyword()
    atoms = Keyword(multi=True)
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
            uploader=User.from_user_popo(source.uploader) if source.uploader is not None else None,

            with_embargo=source.with_embargo,
            published=source.published,
            coauthors=[User.from_user_popo(user) for user in source.coauthors],
            shared_with=[User.from_user_popo(user) for user in source.shared_with],
            comment=source.comment,
            references=[ref.value for ref in source.references],
            datasets=[Dataset.from_dataset_popo(ds) for ds in source.datasets],

            formula=source.formula,
            atoms=list(set(source.atoms)),
            basis_set=source.basis_set,
            xc_functional=source.xc_functional,
            system=source.system,
            crystal_system=source.crystal_system,
            spacegroup=source.spacegroup,
            code_name=source.code_name,
            code_version=source.code_version)

    @classmethod
    def add_upload(cls, source: datamodel.UploadWithMetadata):
        for calc in source.calcs:
            cls.from_calc_with_metadata(calc).save()

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

    @classmethod
    def publish_upload(cls, upload: datamodel.UploadWithMetadata):
        cls.update_by_query(upload.upload_id, 'ctx._source["published"] = true')
        # TODO run update on all calcs with their new metadata

    @classmethod
    def delete_upload(cls, upload_id):
        index = cls._default_index()
        Search(index=index).query('match', upload_id=upload_id).delete()

    @staticmethod
    def es_search(body):
        """ Perform an elasticsearch and not elasticsearch_dsl search on the Calc index. """
        return infrastructure.elastic_client.search(index=config.elastic.index_name, body=body)
