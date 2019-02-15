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
    Object, Boolean, Search, Integer

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

    n_total_energies = Integer()
    n_geometries = Integer()
    geometries = Keyword(multi=True)
    quantities = Keyword(multi=True)

    @classmethod
    def from_calc_with_metadata(cls, source: datamodel.CalcWithMetadata) -> 'Entry':
        entry = Entry(meta=dict(id=source.calc_id))
        entry.update(source)
        return entry

    def update(self, source: datamodel.CalcWithMetadata) -> None:
        self.upload_id = source.upload_id
        self.upload_time = source.upload_time
        self.calc_id = source.calc_id
        self.calc_hash = source.calc_hash
        self.pid = str(source.pid)
        self.mainfile = source.mainfile
        self.files = source.files
        self.uploader = User.from_user_popo(source.uploader) if source.uploader is not None else None

        self.with_embargo = source.with_embargo
        self.published = source.published
        self.coauthors = [User.from_user_popo(user) for user in source.coauthors]
        self.shared_with = [User.from_user_popo(user) for user in source.shared_with]
        self.comment = source.comment
        self.references = [ref.value for ref in source.references]
        self.datasets = [Dataset.from_dataset_popo(ds) for ds in source.datasets]

        self.formula = source.formula
        self.atoms = list(set(source.atoms))
        self.basis_set = source.basis_set
        self.xc_functional = source.xc_functional
        self.system = source.system
        self.crystal_system = source.crystal_system
        self.spacegroup = source.spacegroup
        self.code_name = source.code_name
        self.code_version = source.code_version

        if source.backend is not None:
            quantities = set()
            geometries = set()
            n_total_energies = 0
            n_geometries = 0

            for meta_info, _, value in source.backend._delegate.results.traverse():
                quantities.add(meta_info)
                if meta_info == 'energy_total':
                    n_total_energies += 1
                if meta_info == 'section_system':
                    n_geometries += 1
                if meta_info == 'configuration_raw_gid':
                    geometries.add(value)

            self.geometries = list(geometries)
            self.quantities = list(quantities)
            self.n_total_energies = n_total_energies
            self.n_geometries = n_geometries

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
