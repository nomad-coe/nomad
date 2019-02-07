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
from elasticsearch_dsl import Document, InnerDoc, Keyword, Text, Long, Integer, Date, \
    Nested

from nomad import config, datamodel, coe_repo, infrastructure, datamodel


class AlreadyExists(Exception): pass


class UserData(InnerDoc):
    repository_open_date = Date(format='epoch_millis')
    repository_access_now = Keyword()
    repository_comment = Keyword()

    section_citation = Nested(properties=dict(
        citation_repo_id=Integer(),
        citation_value=Keyword()
    ))

    section_author_info = Nested(properties=dict(
        author_repo_id=Integer(index=True),
        author_first_name=Keyword(),
        author_last_name=Keyword(),
        author_name=Text()
    ))

    section_shared_with = Nested(properties=dict(
        shared_with_repo_id=Keyword(),
        shared_with_first_name=Keyword(),
        shared_with_last_name=Keyword(),
        shared_with_username=Keyword(),
        shared_with_name=Text()
    ))

    section_repository_dataset = Nested(properties=dict(
        dataset_checksum=Keyword(),
        dataset_pid=Keyword(),
        dataset_name=Keyword(),
        dataset_parent_pid=Keyword(),
        dataset_calc_id=Long(),
        dataset_parent_calc_id=Long(),
        section_dataset_doi=Nested(properties=dict(
            dataset_doi_name=Keyword(),
            dataset_doi_id=Long()))
    ))

    def fill_from_coe_repo(self, calc: coe_repo.Calc):
        pass


class CalcData(InnerDoc):
    repository_checksum = Keyword()
    repository_chemical_formula = Keyword()
    repository_parser_id = Keyword()
    repository_atomic_elements = Keyword(store=True)
    repository_atomic_elements_count = Integer(store=True)
    repository_basis_set_type = Keyword(store=True)
    repository_code_version = Keyword(store=True)
    repository_crystal_system = Keyword(store=True)
    repository_program_name = Keyword(store=True)
    repository_spacegroup_nr = Keyword(store=True)
    repository_system_type = Keyword(store=True)
    repository_xc_treatment = Keyword(store=True)


class Calc(InnerDoc):
    # main_file_uri = Keyword()
    # secondary_file_uris = Keyword()
    repository_filepaths = Keyword(index=False)
    # repository_archive_gid = Keyword()
    repository_calc_id = Long(store=True)
    repository_calc_pid = Keyword(store=True)
    upload_id = Long()
    upload_date = Date(format='epoch_millis')
    repository_grouping_checksum = Keyword()

    section_repository_userdata = Nested(UserData)
    section_repository_parserdata = Nested(CalcData)

    section_uploader_info = Nested(properties=dict(
        uploader_repo_id=Keyword(),
        uploader_first_name=Keyword(),
        uploader_last_name=Keyword(),
        uploader_username=Keyword(),
        uploader_name=Text()
    ))


class Entry(Document, datamodel.Entity):
    class Index:
        name = config.elastic.index_name

    calc_id = Keyword()
    upload_id = Keyword()
    section_repository_info = Nested(Calc)

    def __init__(self, upload_id: str, calc_id: str, **kwargs) -> None:
        super().__init__(meta=dict(id=calc_id), **kwargs)
        self.calc_id = calc_id
        self.upload_id = upload_id

    @classmethod
    def from_calc_with_metadata(cls, source: datamodel.CalcWithMetadata) -> 'Entry':
        target = Entry(
            upload_id=source.upload_id,
            calc_id=source.calc_id,
            section_repository_info=Calc(
                section_repository_parserdata=CalcData(
                    repository_checksum=source.calc_hash,
                    repository_chemical_formula=source.chemical_composition,
                    # repository_parser_id
                    repository_atomic_elements=source.atom_labels,
                    repository_atomic_elements_count=len(source.atom_labels),
                    repository_basis_set_type=source.basis_set_type,
                    repository_code_version=source.program_version,
                    repository_crystal_system=source.crystal_system,
                    repository_program_name=source.program_name,
                    repository_spacegroup_nr=source.space_group_number,
                    repository_system_type=source.system_type,
                    repository_xc_treatment=source.XC_functional_name
                ),
                section_repository_userdata=UserData(

                )
            )
        )
        return target

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
