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

''' All generic entry metadata and related classes. '''

from typing import List, Any
from enum import Enum
from cachetools import cached, TTLCache
from elasticsearch_dsl import analyzer, tokenizer

from nomad import metainfo, config
from nomad.metainfo.mongoengine_extension import Mongo, MongoDocument
from nomad.datamodel.metainfo.common import FastAccess
from nomad.metainfo.pydantic_extension import PydanticModel
from nomad.metainfo.elasticsearch_extension import Elasticsearch, material_entry_type, entry_type as es_entry_type

# This is usually defined automatically when the first metainfo definition is evaluated, but
# due to the next imports requireing the m_package already, this would be too late.
m_package = metainfo.Package()

from .results import Results  # noqa
from .data import EntryData  # noqa
from .optimade import OptimadeEntry  # noqa
from .metainfo.simulation.run import Run  # noqa
from .metainfo.workflow import Workflow  # noqa
from .metainfo.measurements import Measurement  # noqa
from .metainfo.tabulartree import TabularTree  # noqa


class AuthLevel(int, Enum):
    '''
    Used to decorate fields with the authorization level required to edit them (using `a_auth_level`).
    * `none`: No authorization required
    * `coauthor`: You must be at least a coauthor of the upload to edit the field.
    * `main_author`: You must be the main author of the upload to edit the field.
    * `admin`: You must be admin to edit the field.
    '''
    none = 0
    coauthor = 1
    main_author = 2
    admin = 3


path_analyzer = analyzer(
    'path_analyzer',
    tokenizer=tokenizer('path_tokenizer', 'pattern', pattern='/'))


def PathSearch():
    return [
        Elasticsearch(_es_field='keyword'),
        Elasticsearch(
            mapping=dict(type='text', analyzer=path_analyzer.to_dict()),
            field='path', _es_field='')]


quantity_analyzer = analyzer(
    'quantity_analyzer',
    tokenizer=tokenizer('quantity_tokenizer', 'pattern', pattern='.'))


def QuantitySearch():
    return [
        Elasticsearch(material_entry_type, _es_field='keyword'),
        Elasticsearch(
            material_entry_type,
            mapping=dict(type='text', analyzer=path_analyzer.to_dict()),
            field='path', _es_field='')]


class Author(metainfo.MSection):
    ''' A person that is author of data in NOMAD or references by NOMAD. '''
    name = metainfo.Quantity(
        type=str,
        derived=lambda user: ('%s %s' % (user.first_name, user.last_name)).strip(),
        a_elasticsearch=[
            Elasticsearch(material_entry_type, _es_field='keyword'),
            Elasticsearch(material_entry_type, mapping='text', field='text', _es_field=''),
            Elasticsearch(suggestion="default")
        ])

    first_name = metainfo.Quantity(type=metainfo.Capitalized)
    last_name = metainfo.Quantity(type=metainfo.Capitalized)
    email = metainfo.Quantity(type=str)

    affiliation = metainfo.Quantity(type=str)
    affiliation_address = metainfo.Quantity(type=str)


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

    m_def = metainfo.Section(a_pydantic=PydanticModel())

    user_id = metainfo.Quantity(
        type=str,
        a_elasticsearch=Elasticsearch(material_entry_type))

    username = metainfo.Quantity(type=str)

    created = metainfo.Quantity(type=metainfo.Datetime)

    repo_user_id = metainfo.Quantity(
        type=str,
        description='Optional, legacy user id from the old NOMAD CoE repository.')

    is_admin = metainfo.Quantity(
        type=bool, derived=lambda user: user.user_id == config.services.admin_user_id)

    is_oasis_admin = metainfo.Quantity(type=bool, default=False)

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


class UserReference(metainfo.Reference):
    '''
    Special metainfo reference type that allows to use user_ids as values. It automatically
    resolves user_ids to User objects. This is done lazily on getting the value.
    '''

    def __init__(self):
        super().__init__(User.m_def)

    def resolve(self, proxy: metainfo.MProxy) -> metainfo.MSection:
        return User.get(user_id=proxy.m_proxy_value)

    def serialize(self, section: metainfo.MSection, quantity_def: metainfo.Quantity, value: Any) -> Any:
        return value.user_id


user_reference = UserReference()


class AuthorReference(metainfo.Reference):
    '''
    Special metainfo reference type that allows to use either user_ids or direct author
    information as values. It automatically resolves user_ids to User objects and author
    data into Author objects.
    '''

    def __init__(self):
        super().__init__(Author.m_def)

    def resolve(self, proxy: metainfo.MProxy) -> metainfo.MSection:
        proxy_value = proxy.m_proxy_value
        if isinstance(proxy_value, str):
            return User.get(user_id=proxy.m_proxy_value)
        elif isinstance(proxy_value, dict):
            return Author.m_from_dict(proxy_value)
        else:
            raise metainfo.MetainfoReferenceError()

    def serialize(self, section: metainfo.MSection, quantity_def: metainfo.Quantity, value: Any) -> Any:
        if isinstance(value, User):
            return value.user_id
        elif isinstance(value, Author):
            return value.m_to_dict()
        else:
            raise metainfo.MetainfoReferenceError()


author_reference = AuthorReference()


class Dataset(metainfo.MSection):
    ''' A Dataset is attached to one or many entries to form a set of data.

    Args:
        dataset_id: The unique identifier for this dataset as a string. It should be
            a randomly generated UUID, similar to other nomad ids.
        dataset_name: The human readable name of the dataset as string. The dataset name must be
            unique for the user.
        user_id: The unique user_id of the owner and creator of this dataset. The owner
            must not change after creation.
        doi: The optional Document Object Identifier (DOI) associated with this dataset.
            Nomad can register DOIs that link back to the respective representation of
            the dataset in the nomad UI. This quantity holds the string representation of
            this DOI. There is only one per dataset. The DOI is just the DOI name, not its
            full URL, e.g. "10.17172/nomad/2019.10.29-1".
        pid: The original NOMAD CoE Repository dataset PID. Old DOIs still reference
            datasets based on this id. Is not used for new datasets.
        dataset_create_time: The date when the dataset was first created.
        dataset_modified_time: The date when the dataset was last modified. An owned dataset
            can only be extended after a DOI was assigned. A foreign dataset cannot be changed
            once a DOI was assigned.
        dataset_type: The type determined if a dataset is owned, i.e. was created by
            the authors of the contained entries; or if a dataset is foreign,
            i.e. it was created by someone not necessarily related to the entries.
    '''
    m_def = metainfo.Section(a_mongo=MongoDocument(), a_pydantic=PydanticModel())

    dataset_id = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(primary_key=True),
        a_elasticsearch=Elasticsearch(material_entry_type))
    dataset_name = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(index=True),
        a_elasticsearch=[
            Elasticsearch(material_entry_type),
            Elasticsearch(suggestion="default"),
        ])
    user_id = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(index=True))
    doi = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(index=True),
        a_elasticsearch=Elasticsearch(material_entry_type))
    pid = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(index=True))
    dataset_create_time = metainfo.Quantity(
        type=metainfo.Datetime,
        a_mongo=Mongo(index=True),
        a_elasticsearch=Elasticsearch())
    dataset_modified_time = metainfo.Quantity(
        type=metainfo.Datetime,
        a_mongo=Mongo(index=True),
        a_elasticsearch=Elasticsearch())
    dataset_type = metainfo.Quantity(
        type=metainfo.MEnum('owned', 'foreign'),
        a_mongo=Mongo(index=True),
        a_elasticsearch=Elasticsearch())
    query = metainfo.Quantity(
        type=metainfo.JSON, a_mongo=Mongo())
    entries = metainfo.Quantity(
        type=str, shape=['*'], a_mongo=Mongo())


class DatasetReference(metainfo.Reference):
    '''
    Special metainfo reference type that allows to use dataset_ids as values. It automatically
    resolves dataset_ids to Dataset objects. This is done lazily on getting the value.
    '''

    def __init__(self):
        super().__init__(Dataset.m_def)

    def resolve(self, proxy: metainfo.MProxy) -> metainfo.MSection:
        return Dataset.m_def.a_mongo.get(dataset_id=proxy.m_proxy_value)

    def serialize(self, section: metainfo.MSection, quantity_def: metainfo.Quantity, value: Any) -> Any:
        if isinstance(value, metainfo.MProxy):
            return value.m_proxy_value
        else:
            return value.dataset_id


dataset_reference = DatasetReference()


class EditableUserMetadata(metainfo.MCategory):
    ''' NOMAD entry metadata quantities that can be edited by the user before or after publish. '''
    pass


class MongoUploadMetadata(metainfo.MCategory):
    ''' The field is defined on the Upload mongo document. '''
    pass


class MongoEntryMetadata(metainfo.MCategory):
    ''' The field is defined on the Entry mongo document. '''
    pass


class MongoSystemMetadata(metainfo.MCategory):
    '''
    The field is managed directly by the system/process (or derived from data managed by the
    system/process), and should never be updated from an :class:`EntryMetadata` object.
    '''
    pass


class DomainMetadata(metainfo.MCategory):
    ''' NOMAD entry quantities that are determined by the uploaded data. '''
    pass


def derive_origin(entry: 'EntryMetadata') -> str:
    if entry.external_db is not None:
        return str(entry.external_db)

    if entry.main_author:
        return entry.main_author.name

    return None


def derive_authors(entry: 'EntryMetadata') -> List[User]:
    if entry.external_db == 'EELS Data Base':
        return list(entry.entry_coauthors)

    authors: List[User] = [entry.main_author]
    if entry.coauthors:
        authors.extend(entry.coauthors)
    if entry.entry_coauthors:
        authors.extend(entry.entry_coauthors)
    return authors


class EntryMetadata(metainfo.MSection):
    '''
    Attributes:
        upload_id: The id of the upload (random UUID).
        upload_name: The user provided upload name.
        upload_create_time: The time that the upload was created
        entry_id: The unique mainfile based entry id.
        entry_hash: The raw file content based checksum/hash of this entry.
        entry_create_time: The time that the entry was created
        last_edit_time: The date and time the user metadata was last edited.
        parser_name: The NOMAD parser used for the last processing.
        mainfile: The path to the mainfile from the root directory of the uploaded files.
        mainfile_key: Key used to differentiate between different *child entries* of an entry.
            For parent entries and entries that do not have any children, the value should
            be empty.
        files: A list of all files, relative to upload.
        pid: The unique, sequentially enumerated, integer PID that was used in the legacy
            NOMAD CoE. It allows to resolve URLs of the old NOMAD CoE Repository.
        raw_id: The code specific identifier extracted from the entry's raw files by the parser,
            if supported.
        external_id: A user provided external id. Usually the id for an entry in an external
            database where the data was imported from.
        published: Indicates if the entry is published.
        publish_time: The time when the upload was published.
        with_embargo: Indicated if this entry is under an embargo. Entries with embargo are
            only visible to the main author, the upload coauthors, and the upload reviewers.
        license: A short license description (e.g. CC BY 4.0), that refers to the
            license of this entry.
        processed: Boolean indicating if this entry was successfully processed and archive
            data and entry metadata is available.
        last_processing_time: The date and time of the last processing.
        processing_errors: Errors that occured during processing.
        nomad_version: A string that describes the version of the nomad software that was
            used to do the last successful processing.
        nomad_commit: The NOMAD commit used for the last processing.
        comment: An arbitrary string with user provided information about the entry.
        references: A list of URLs for resources that are related to the entry.
        external_db: The repository or external database where the original data resides.
        origin: A short human readable description of the entries origin. Usually it is the
            handle of an external database/repository or the name of the main author.
        main_author: Id of the main author of this entry.
        coauthors: A user provided list of co-authors for the whole upload. These can view
            and edit the upload when in staging, and view it also if it is embargoed.
        entry_coauthors: Ids of all co-authors (excl. the main author and upload coauthors)
            specified on the entry level, rather than on the upload level. They are shown
            as authors of this entry alongside its main author and upload coauthors.
        reviewers: Ids of users who can review the upload which this entry belongs to. Like
            the main author and the upload coauthors, reviewers can find, see, and download
            all data from the upload and all its entries, even if it is in staging or has
            an embargo.
        datasets: Ids of all datasets that this entry appears in
    '''
    m_def = metainfo.Section(label='Metadata')

    upload_id = metainfo.Quantity(
        type=str, categories=[MongoUploadMetadata],
        description='The persistent and globally unique identifier for the upload of the entry',
        a_elasticsearch=Elasticsearch(material_entry_type, metrics=dict(n_uploads='cardinality')))

    upload_name = metainfo.Quantity(
        type=str, categories=[MongoUploadMetadata, EditableUserMetadata],
        description='The user provided upload name',
        a_elasticsearch=Elasticsearch())

    upload_create_time = metainfo.Quantity(
        type=metainfo.Datetime, categories=[MongoUploadMetadata, EditableUserMetadata],
        description='The date and time when the upload was created in nomad',
        a_auth_level=AuthLevel.admin,
        a_elasticsearch=Elasticsearch(material_entry_type))

    entry_id = metainfo.Quantity(
        type=str,
        description='A persistent and globally unique identifier for the entry',
        categories=[MongoEntryMetadata, MongoSystemMetadata],
        a_elasticsearch=Elasticsearch(material_entry_type, metrics=dict(n_entries='cardinality')))

    entry_name = metainfo.Quantity(
        type=str,
        description='A brief human readable name for the entry.',
        a_elasticsearch=[
            Elasticsearch(material_entry_type, _es_field='keyword'),
            Elasticsearch(
                material_entry_type, field='prefix',
                es_query='match_phrase_prefix', mapping='text', _es_field='')
        ])

    entry_type = metainfo.Quantity(
        type=str,
        description='The main schema definition. This is the name of the section used for data.',
        a_elasticsearch=Elasticsearch(material_entry_type))

    calc_id = metainfo.Quantity(
        type=str, description='Legacy field name, use `entry_id` instead.',
        derived=lambda entry: entry.entry_id,
        a_elasticsearch=Elasticsearch(material_entry_type))

    entry_hash = metainfo.Quantity(
        # Note: This attribute is not stored in ES
        type=str,
        description='A raw file content based checksum/hash',
        categories=[MongoEntryMetadata])

    entry_create_time = metainfo.Quantity(
        type=metainfo.Datetime, categories=[MongoEntryMetadata, MongoSystemMetadata, EditableUserMetadata],
        description='The date and time when the entry was created in nomad',
        a_auth_level=AuthLevel.admin,
        a_elasticsearch=Elasticsearch(material_entry_type))

    last_edit_time = metainfo.Quantity(
        # Note: This attribute is not stored in ES
        type=metainfo.Datetime, categories=[MongoEntryMetadata],
        description='The date and time the user metadata was last edited.')

    parser_name = metainfo.Quantity(
        type=str, categories=[MongoEntryMetadata, MongoSystemMetadata],
        description='The NOMAD parser used for the last processing',
        a_elasticsearch=Elasticsearch())

    mainfile = metainfo.Quantity(
        type=str, categories=[MongoEntryMetadata, MongoSystemMetadata],
        description='The path to the mainfile from the root directory of the uploaded files',
        a_elasticsearch=PathSearch())

    mainfile_key = metainfo.Quantity(
        type=str, categories=[MongoEntryMetadata, MongoSystemMetadata],
        description='''
            Key used to differentiate between different *child entries* of an entry.
            For parent entries and entries that do not have any children, the value should
            be empty.
        ''',
        a_elasticsearch=PathSearch())

    files = metainfo.Quantity(
        type=str, shape=['0..*'],
        description='''
            The paths to the files within the upload that belong to this entry.
            All files within the same directory as the entry's mainfile are considered the
            auxiliary files that belong to the entry.
        ''',
        a_elasticsearch=PathSearch())

    pid = metainfo.Quantity(
        type=str,
        description='''
            The unique, sequentially enumerated, integer PID that was used in the legacy
            NOMAD CoE. It allows to resolve URLs of the old NOMAD CoE Repository.
        ''',
        categories=[MongoEntryMetadata],
        a_elasticsearch=Elasticsearch(es_entry_type))

    raw_id = metainfo.Quantity(
        type=str,
        description='''
            The code specific identifier extracted from the entry's raw files by the parser,
            if supported.
        ''',
        a_elasticsearch=Elasticsearch(es_entry_type))

    external_id = metainfo.Quantity(
        type=str, categories=[MongoEntryMetadata, EditableUserMetadata],
        description='''
            A user provided external id. Usually the id for an entry in an external database
            where the data was imported from.
        ''',
        a_elasticsearch=Elasticsearch())

    published = metainfo.Quantity(
        type=bool, default=False,
        description='Indicates if the entry is published',
        categories=[MongoUploadMetadata],
        a_elasticsearch=Elasticsearch(material_entry_type))

    publish_time = metainfo.Quantity(
        type=metainfo.Datetime, categories=[MongoUploadMetadata, EditableUserMetadata],
        description='The date and time when the upload was published in nomad',
        a_auth_level=AuthLevel.admin,
        a_elasticsearch=Elasticsearch(material_entry_type))

    with_embargo = metainfo.Quantity(
        type=bool, default=False,
        categories=[MongoUploadMetadata, MongoSystemMetadata],
        description='Indicated if this entry is under an embargo',
        a_elasticsearch=Elasticsearch(material_entry_type))

    embargo_length = metainfo.Quantity(
        # Note: This attribute is not stored in ES
        type=int, categories=[MongoUploadMetadata, EditableUserMetadata],
        description='The length of the requested embargo period, in months')

    license = metainfo.Quantity(
        # Note: This attribute is not stored in ES
        type=str,
        description='''
            A short license description (e.g. CC BY 4.0), that refers to the
            license of this entry.
        ''',
        default='CC BY 4.0',
        categories=[MongoUploadMetadata, EditableUserMetadata],
        a_auth_level=AuthLevel.admin)

    processed = metainfo.Quantity(
        type=bool, default=False, categories=[MongoEntryMetadata, MongoSystemMetadata],
        description='Indicates that the entry is successfully processed.',
        a_elasticsearch=Elasticsearch())

    last_processing_time = metainfo.Quantity(
        type=metainfo.Datetime,
        description='The date and time of the last processing.',
        categories=[MongoEntryMetadata],
        a_elasticsearch=Elasticsearch())

    processing_errors = metainfo.Quantity(
        type=str, shape=['*'], description='Errors that occured during processing',
        a_elasticsearch=Elasticsearch())

    nomad_version = metainfo.Quantity(
        type=str,
        description='The NOMAD version used for the last processing',
        categories=[MongoEntryMetadata],
        a_elasticsearch=Elasticsearch())

    nomad_commit = metainfo.Quantity(
        type=str,
        description='The NOMAD commit used for the last processing',
        categories=[MongoEntryMetadata],
        a_elasticsearch=Elasticsearch())

    comment = metainfo.Quantity(
        type=str, categories=[MongoEntryMetadata, EditableUserMetadata],
        description='A user provided comment for this entry',
        a_elasticsearch=Elasticsearch(mapping='text'))

    references = metainfo.Quantity(
        type=str, shape=['0..*'], categories=[MongoEntryMetadata, EditableUserMetadata],
        description='User provided references (URLs) for this entry',
        a_elasticsearch=Elasticsearch())

    external_db = metainfo.Quantity(
        type=metainfo.MEnum('EELS Data Base', 'Materials Project', 'AFLOW', 'OQMD', 'Kyoto Phonopy Database'),
        categories=[MongoUploadMetadata, EditableUserMetadata],
        description='The repository or external database where the original data resides',
        a_elasticsearch=Elasticsearch(material_entry_type))

    origin = metainfo.Quantity(
        type=str,
        description='''
            A short human readable description of the entries origin. Usually it is the
            handle of an external database/repository or the name of the main author.
        ''',
        derived=derive_origin,
        a_elasticsearch=Elasticsearch(material_entry_type))

    main_author = metainfo.Quantity(
        type=user_reference, categories=[MongoUploadMetadata, EditableUserMetadata],
        description='The main author of the entry',
        a_auth_level=AuthLevel.admin,
        a_elasticsearch=Elasticsearch(material_entry_type))

    coauthors = metainfo.Quantity(
        # Note: This attribute is not stored in ES
        type=author_reference, shape=['0..*'], default=[], categories=[MongoUploadMetadata, EditableUserMetadata],
        description='''
            A user provided list of co-authors for the whole upload. These can view and edit the
            upload when in staging, and view it also if it is embargoed.
        ''')

    entry_coauthors = metainfo.Quantity(
        # Note: This attribute is not stored in ES
        type=author_reference, shape=['0..*'], default=[], categories=[MongoEntryMetadata],
        description='''
            A user provided list of co-authors specific for this entry. This is a legacy field,
            for new uploads, coauthors should be specified on the upload level only.
        ''')

    reviewers = metainfo.Quantity(
        # Note: This attribute is not stored in ES
        type=user_reference, shape=['0..*'], default=[], categories=[MongoUploadMetadata, EditableUserMetadata],
        description='''
            A user provided list of reviewers. Reviewers can see the whole upload, also if
            it is unpublished or embargoed
        ''')

    authors = metainfo.Quantity(
        type=author_reference, shape=['0..*'],
        description='All authors (main author and co-authors)',
        derived=derive_authors,
        a_elasticsearch=Elasticsearch(material_entry_type, metrics=dict(n_authors='cardinality')))

    writers = metainfo.Quantity(
        type=user_reference, shape=['0..*'],
        description='All writers (main author, upload coauthors)',
        derived=lambda entry: ([entry.main_author] if entry.main_author is not None else []) + entry.coauthors,
        a_elasticsearch=Elasticsearch(material_entry_type))

    viewers = metainfo.Quantity(
        type=user_reference, shape=['0..*'],
        description='All viewers (main author, upload coauthors, and reviewers)',
        derived=lambda entry: ([entry.main_author] if entry.main_author is not None else []) + entry.coauthors + entry.reviewers,
        a_elasticsearch=Elasticsearch(material_entry_type))

    datasets = metainfo.Quantity(
        type=dataset_reference, shape=['0..*'], default=[],
        categories=[MongoEntryMetadata, EditableUserMetadata],
        description='A list of user curated datasets this entry belongs to.',
        a_elasticsearch=Elasticsearch(material_entry_type))

    optimade = metainfo.SubSection(
        sub_section=OptimadeEntry,
        description='Metadata used for the optimade API.',
        a_elasticsearch=Elasticsearch(es_entry_type))

    domain = metainfo.Quantity(
        type=metainfo.MEnum('dft', 'ems'),
        description='The material science domain',
        a_elasticsearch=Elasticsearch(material_entry_type))

    n_quantities = metainfo.Quantity(
        type=int, default=0, description='Number of metainfo quantities parsed from the entry.',
        a_elasticsearch=Elasticsearch(metrics=dict(n_quantities='sum')))

    quantities = metainfo.Quantity(
        type=str, shape=['0..*'],
        description='All quantities that are used by this entry.',
        a_elasticsearch=QuantitySearch())

    sections = metainfo.Quantity(
        type=str, shape=['*'],
        description='All sections that are present in this entry.',
        a_elasticsearch=Elasticsearch(material_entry_type))

    def apply_archvie_metadata(self, archive):
        quantities = set()
        sections = set()
        n_quantities = 0

        section_paths = {}

        def get_section_path(section):
            section_path = section_paths.get(section)
            if section_path is None:
                parent = section.m_parent
                if parent:
                    parent_path = get_section_path(parent)
                    if parent_path == '':
                        section_path = section.m_parent_sub_section.name
                    else:
                        section_path = f'{parent_path}.{section.m_parent_sub_section.name}'
                else:
                    section_path = ''
                section_paths[section] = section_path
                quantities.add(section_path)

            return section_path

        for section, property_def, _ in archive.m_traverse():
            sections.add(section.m_def)

            quantity_path = f'{get_section_path(section)}.{property_def.name}'
            quantities.add(quantity_path)
            n_quantities += 1

        self.quantities = list(quantities)
        self.quantities.sort()
        self.sections = [section.qualified_name() for section in sections]
        self.sections.sort()
        self.n_quantities = n_quantities


class EntryArchive(metainfo.MSection):
    m_def = metainfo.Section(label='Entry')

    entry_id = metainfo.Quantity(
        type=str, description='The unique primary id for this entry.',
        derived=lambda entry: entry.metadata.entry_id)

    run = metainfo.SubSection(sub_section=Run, repeats=True)
    measurement = metainfo.SubSection(sub_section=Measurement, repeats=True)

    data = metainfo.SubSection(sub_section=EntryData)

    workflow = metainfo.SubSection(sub_section=Workflow, repeats=True, categories=[FastAccess])
    metadata = metainfo.SubSection(
        sub_section=EntryMetadata, categories=[FastAccess],
        a_elasticsearch=Elasticsearch())

    processing_logs = metainfo.Quantity(
        type=Any, shape=['0..*'],
        description='The processing logs for this entry as a list of structlog entries.')

    results = metainfo.SubSection(
        sub_section=Results,
        categories=[FastAccess],
        a_elasticsearch=Elasticsearch(auto_include_subsections=True))

    tabular_tree = metainfo.SubSection(sub_section=TabularTree, repeats=False)
    definitions = metainfo.SubSection(sub_section=metainfo.Package)

    def normalize(self, archive, logger):
        if not archive.metadata.entry_type:
            if archive.definitions is not None:
                archive.metadata.entry_type = 'Schema'

    def m_update_from_dict(self, dct) -> None:
        super().m_update_from_dict(dct)
        if self.definitions is not None:
            self.definitions.archive = self


m_package.__init_metainfo__()
