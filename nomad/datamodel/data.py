#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import os.path

from typing import Any
from cachetools import cached, TTLCache
from nomad.metainfo.metainfo import (
    predefined_datatypes, Category, MCategory, MSection, Quantity, Reference, MetainfoReferenceError,
    MProxy, Capitalized, Section, Datetime)
from nomad import config
from nomad.metainfo.pydantic_extension import PydanticModel
from nomad.metainfo.elasticsearch_extension import Elasticsearch, material_entry_type


class ArchiveSection(MSection):
    '''
    Base class for sections in a NOMAD archive. Provides a framework for custom
    section normalization via the `normalize` function.
    '''

    def normalize(self, archive, logger):
        '''
        Is called during entry normalization. If you overwrite this with custom
        normalization code, make sure to call `super(YourClass, self).normalize(archive, logger)`.
        Otherwise, not all normalize functions might be called for section definitions
        with multiple base-classes.

        Arguments:
            archive: The whole archive that is normalized.
            logger: The structlog logger used during normalization.
        '''
        pass


class EntryDataCategory(MCategory):
    pass


class ElnIntegrationCategory(EntryDataCategory):
    m_def = Category(label='Third-party ELN Integration', categories=[EntryDataCategory])


class BasicElnCategory(EntryDataCategory):
    m_def = Category(label='Basic ELN', categories=[EntryDataCategory])


class ElnExampleCategory(EntryDataCategory):
    m_def = Category(label='Example ELNs', categories=[EntryDataCategory])


class UseCaseElnCategory(EntryDataCategory):
    m_def = Category(label='Use-cases', categories=[EntryDataCategory])


class WorkflowsElnCategory(EntryDataCategory):
    m_def = Category(label='Workflows', categories=[EntryDataCategory])


class EntryData(ArchiveSection):
    '''
    An empty base section definition. This can be used to add new top-level sections
    to an entry.
    '''

    def normalize(self, archive, logger):
        super(EntryData, self).normalize(archive, logger)

        from nomad.datamodel.results import Results

        archive.metadata.entry_type = self.m_def.name
        if archive.metadata.entry_name is None and archive.metadata.mainfile:
            archive.metadata.entry_name = os.path.basename(archive.metadata.mainfile)

        if not archive.results:
            archive.results = Results()


class Author(MSection):
    ''' A person that is author of data in NOMAD or references by NOMAD. '''
    name = Quantity(
        type=str,
        derived=lambda user: ('%s %s' % (user.first_name, user.last_name)).strip(),
        a_elasticsearch=[
            Elasticsearch(material_entry_type, _es_field='keyword'),
            Elasticsearch(material_entry_type, mapping='text', field='text', _es_field=''),
            Elasticsearch(suggestion="default")
        ])

    first_name = Quantity(type=Capitalized)
    last_name = Quantity(type=Capitalized)
    email = Quantity(type=str)

    affiliation = Quantity(type=str)
    affiliation_address = Quantity(type=str)


class User(Author):
    ''' A NOMAD user.

    Typically a NOMAD user has a NOMAD account. The user related data is managed by
    NOMAD keycloak user-management system. Users are used to denote authors,
    reviewers, and owners of datasets.

    Args:
        user_id: The unique, persistent keycloak UUID
        username: The unique, persistent, user chosen username
        first_name: The users first name (including all other given names)
        last_name: The users last name
        affiliation: The name of the company and institutes the user identifies with
        affiliation_address: The address of the given affiliation
        created: The time the account was created
        repo_user_id: The id that was used to identify this user in the NOMAD CoE Repository
        is_admin: Bool that indicated, iff the user the use admin user
    '''

    m_def = Section(a_pydantic=PydanticModel())

    user_id = Quantity(
        type=str,
        a_elasticsearch=Elasticsearch(material_entry_type))

    username = Quantity(type=str)

    created = Quantity(type=Datetime)

    repo_user_id = Quantity(
        type=str,
        description='Optional, legacy user id from the old NOMAD CoE repository.')

    is_admin = Quantity(
        type=bool, derived=lambda user: user.user_id == config.services.admin_user_id)

    is_oasis_admin = Quantity(type=bool, default=False)

    @staticmethod
    @cached(cache=TTLCache(maxsize=2048, ttl=24 * 3600))
    def get(*args, **kwargs) -> 'User':
        from nomad import infrastructure
        return infrastructure.user_management.get_user(*args, **kwargs)  # type: ignore

    def full_user(self) -> 'User':
        ''' Returns a User object with all attributes loaded from the user management system. '''
        from nomad import infrastructure
        assert self.user_id is not None
        return infrastructure.user_management.get_user(user_id=self.user_id)  # type: ignore


class UserReference(Reference):
    '''
    Special metainfo reference type that allows to use user_ids as values. It automatically
    resolves user_ids to User objects. This is done lazily on getting the value.
    '''

    def __init__(self):
        super().__init__(User.m_def)

    def resolve(self, proxy: MProxy) -> MSection:
        return User.get(user_id=proxy.m_proxy_value)

    def serialize_type(self, type_data):
        return dict(type_kind='User', type_data=self.target_section_def.name)

    @classmethod
    def deserialize_type(cls, type_kind, type_data, section):
        if type_kind == 'User':
            return user_reference
        return None

    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        return value.user_id


user_reference = UserReference()
predefined_datatypes["User"] = user_reference


class AuthorReference(Reference):
    '''
    Special metainfo reference type that allows to use either user_ids or direct author
    information as values. It automatically resolves user_ids to User objects and author
    data into Author objects.
    '''

    def __init__(self):
        super().__init__(Author.m_def)

    def resolve(self, proxy: MProxy) -> MSection:
        proxy_value = proxy.m_proxy_value
        if isinstance(proxy_value, str):
            return User.get(user_id=proxy.m_proxy_value)
        elif isinstance(proxy_value, dict):
            return Author.m_from_dict(proxy_value)
        else:
            raise MetainfoReferenceError()

    def serialize_type(self, type_data):
        return dict(type_kind='Author', type_data=self.target_section_def.name)

    @classmethod
    def deserialize_type(cls, type_kind, type_data, section):
        if type_kind == 'Author':
            return author_reference
        return None

    def serialize(self, section: MSection, quantity_def: Quantity, value: Any) -> Any:
        if isinstance(value, User):
            return value.user_id
        elif isinstance(value, Author):
            return value.m_to_dict()
        else:
            raise MetainfoReferenceError()


author_reference = AuthorReference()
predefined_datatypes["Author"] = author_reference
