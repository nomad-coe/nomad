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

from elasticsearch_dsl import Document, InnerDoc, Keyword, Text, Long, Integer, Date, \
    Nested

from nomad import config, datamodel, coe_repo


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


class Calc(InnerDoc, datamodel.Entity):
    main_file_uri = Keyword()
    secondary_file_uris = Keyword()
    repository_filepaths = Keyword(index=False)
    repository_archive_gid = Keyword()
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

    @staticmethod
    def create_from(source: datamodel.Calc):
        coe_calc = source.to(coe_repo.Calc)
        calc = Calc()
        calc.main_file_uri = coe_calc.mainfile


class Entry(Document):
    class Index:
        name = config.elastic.coe_repo_calcs_index_name

    section_repository_info = Nested(Calc)
