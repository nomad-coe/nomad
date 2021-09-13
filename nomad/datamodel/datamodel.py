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

from typing import Any
from cachetools import cached, TTLCache
from elasticsearch_dsl import Keyword, Text, analyzer, tokenizer
import ase.data

from nomad import metainfo, config
from nomad.metainfo.search_extension import Search
from nomad.metainfo.elastic_extension import ElasticDocument
from nomad.metainfo.mongoengine_extension import Mongo, MongoDocument
from nomad.datamodel.metainfo.common_dft import FastAccess
from nomad.metainfo.pydantic_extension import PydanticModel

from .dft import DFTMetadata
from .ems import EMSMetadata
from .qcms import QCMSMetadata

# This is usually defined automatically when the first metainfo definition is evaluated, but
# due to the next imports requireing the m_package already, this would be too late.
m_package = metainfo.Package()

from .encyclopedia import EncyclopediaMetadata  # noqa
from .metainfo.common_dft import Run, Workflow  # noqa
from .metainfo.common_experimental import Experiment  # noqa
from .metainfo.common_qcms import QuantumCMS  # noqa


def _only_atoms(atoms):
    numbers = [ase.data.atomic_numbers[atom] for atom in atoms]
    only_atoms = [ase.data.chemical_symbols[number] for number in sorted(numbers)]
    return ''.join(only_atoms)


path_analyzer = analyzer(
    'path_analyzer',
    tokenizer=tokenizer('path_tokenizer', 'pattern', pattern='/'))


class Author(metainfo.MSection):
    ''' A person that is author of data in NOMAD or references by NOMAD. '''
    name = metainfo.Quantity(
        type=str,
        derived=lambda user: ('%s %s' % (user.first_name, user.last_name)).strip(),
        a_search=Search(mapping=Text(fields={'keyword': Keyword()})))

    first_name = metainfo.Quantity(type=metainfo.Capitalized)
    last_name = metainfo.Quantity(type=metainfo.Capitalized)
    email = metainfo.Quantity(
        type=str,
        a_elastic=dict(mapping=Keyword),  # TODO remove?
        a_search=Search())

    affiliation = metainfo.Quantity(type=str)
    affiliation_address = metainfo.Quantity(type=str)


class User(Author):
    ''' A NOMAD user.

    Typically a NOMAD user has a NOMAD account. The user related data is managed by
    NOMAD keycloak user-management system. Users are used to denote uploaders, authors,
    people to shared data with embargo with, and owners of datasets.

    Args:
        user_id: The unique, persistent keycloak UUID
        username: The unique, persistent, user chosen username
        first_name: The users first name (including all other given names)
        last_name: The users last name
        affiliation: The name of the company and institutes the user identifies with
        affiliation_address: The address of the given affiliation
        create: The time the account was created
        repo_user_id: The id that was used to identify this user in the NOMAD CoE Repository
        is_admin: Bool that indicated, iff the user the use admin user
    '''

    m_def = metainfo.Section(a_pydantic=PydanticModel())

    user_id = metainfo.Quantity(
        type=str,
        a_search=Search())

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
        return infrastructure.keycloak.get_user(*args, **kwargs)  # type: ignore

    def full_user(self) -> 'User':
        ''' Returns a User object with all attributes loaded from the user management system. '''
        from nomad import infrastructure
        assert self.user_id is not None
        return infrastructure.keycloak.get_user(user_id=self.user_id)  # type: ignore


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
        name: The human readable name of the dataset as string. The dataset name must be
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
        created: The date when the dataset was first created.
        modified: The date when the dataset was last modified. An owned dataset can only
            be extended after a DOI was assigned. A foreign dataset cannot be changed
            once a DOI was assigned.
        dataset_type: The type determined if a dataset is owned, i.e. was created by
            the uploader/owner of the contained entries; or if a dataset is foreign,
            i.e. it was created by someone not necessarily related to the entries.
    '''
    m_def = metainfo.Section(a_mongo=MongoDocument(), a_pydantic=PydanticModel())

    dataset_id = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(primary_key=True),
        a_search=Search())
    name = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(index=True),
        a_search=Search())
    user_id = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(index=True))
    doi = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(index=True),
        a_search=Search())
    pid = metainfo.Quantity(
        type=str,
        a_mongo=Mongo(index=True))
    created = metainfo.Quantity(
        type=metainfo.Datetime,
        a_mongo=Mongo(index=True),
        a_search=Search())
    modified = metainfo.Quantity(
        type=metainfo.Datetime,
        a_mongo=Mongo(index=True),
        a_search=Search())
    dataset_type = metainfo.Quantity(
        type=metainfo.MEnum('owned', 'foreign'),
        a_mongo=Mongo(index=True),
        a_search=Search())
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


class UserProvidableMetadata(metainfo.MCategory):
    ''' NOMAD entry metadata quantities that can be determined by the user, e.g. via nomad.yaml. '''


class EditableUserMetadata(metainfo.MCategory):
    ''' NOMAD entry metadata quantities that can be edited by the user after publish. '''
    m_def = metainfo.Category(categories=[UserProvidableMetadata])


class OasisMetadata(metainfo.MCategory):
    ''' NOMAD entry metadata quantities that can be provided by an OASIS. '''
    m_def = metainfo.Category(categories=[EditableUserMetadata])


class MongoMetadata(metainfo.MCategory):
    ''' NOMAD entry quantities that are stored in mongodb and not necessarely in the archive. '''
    pass


class DomainMetadata(metainfo.MCategory):
    ''' NOMAD entry quantities that are determined by the uploaded data. '''
    pass


def derive_origin(entry):
    if entry.external_db is not None:
        return str(entry.external_db)

    if entry.uploader:
        return entry.uploader.name

    return None


def derive_authors(entry):
    uploaders = []
    if entry.uploader is not None and entry.external_db is None:
        uploaders = [entry.uploader]
    return uploaders + entry.coauthors


class EntryMetadata(metainfo.MSection):
    '''
    Attributes:
        upload_id: The ``upload_id`` of the calculations upload (random UUID).
        calc_id: The unique mainfile based calculation id.
        calc_hash: The raw file content based checksum/hash of this calculation.
        pid: The unique persistent id of this calculation.
        mainfile: The upload relative mainfile path.
        domain: Must be the key for a registered domain. This determines which actual
            subclass is instantiated.

        files: A list of all files, relative to upload.
        upload_time: The time when the calc was uploaded.
        uploader: An object describing the uploading user, has at least ``user_id``
        processed: Boolean indicating if this calc was successfully processed and archive
            data and calc metadata is available.
        last_processing: A datatime with the time of the last successful processing.
        nomad_version: A string that describes the version of the nomad software that was
            used to do the last successful processing.

        comment: An arbitrary string with user provided information about the entry.
        references: A list of URLs for resources that are related to the entry.
        uploader: Id of the uploader of this entry.
        coauthors: Ids of all co-authors (excl. the uploader) of this entry. Co-authors are
            shown as authors of this entry alongside its uploader.
        shared_with: Ids of all users that this entry is shared with. These users can find,
            see, and download all data for this entry, even if it is in staging or
            has an embargo.
        with_embargo: Entries with embargo are only visible to the uploader, the admin
            user, and users the entry is shared with (see shared_with).
        upload_time: The time that this entry was uploaded
        datasets: Ids of all datasets that this entry appears in
    '''
    m_def = metainfo.Section(
        a_elastic=ElasticDocument(index_name=config.elastic.index_name, id=lambda x: x.calc_id))

    upload_id = metainfo.Quantity(
        type=str,
        description='The persistent and globally unique identifier for the upload of the entry',
        a_search=Search(
            many_or='append', group='uploads_grouped', metric_name='uploads', metric='cardinality'))

    calc_id = metainfo.Quantity(
        type=str,
        description='A persistent and globally unique identifier for the entry',
        categories=[OasisMetadata],
        a_search=Search(many_or='append'))

    calc_hash = metainfo.Quantity(
        type=str,
        description='A raw file content based checksum/hash',
        categories=[MongoMetadata],
        a_search=Search(
            many_or='append', metric_name='unique_entries', metric='cardinality'))

    mainfile = metainfo.Quantity(
        type=str,
        description='The path to the mainfile from the root directory of the uploaded files',
        a_search=[
            Search(
                description='Search within the mainfile path.',
                mapping=Text(multi=True, analyzer=path_analyzer, fields={'keyword': Keyword()}),
                many_or='append', search_field='mainfile.keyword'),
            Search(
                description='Search for the exact mainfile.',
                many_and='append', name='mainfile_path', search_field='mainfile.keyword')])

    files = metainfo.Quantity(
        type=str, shape=['0..*'],
        description='''
        The paths to the files within the upload that belong to this entry.
        All files within the same directory as the entry's mainfile are considered the
        auxiliary files that belong to the entry.
        ''',
        a_search=[
            Search(
                description='Search within the paths.', name='path',
                mapping=Text(
                    multi=True, analyzer=path_analyzer, fields={'keyword': Keyword()})
            ),
            Search(
                description='Search for exact paths.',
                many_or='append', name='files', search_field='files.keyword')])

    pid = metainfo.Quantity(
        type=str,
        description='''
        The unique, sequentially enumerated, integer PID that was used in the legacy
        NOMAD CoE. It allows to resolve URLs of the old NOMAD CoE Repository.''',
        categories=[MongoMetadata],
        a_search=Search(many_or='append'))

    raw_id = metainfo.Quantity(
        type=str,
        description='''
        The code specific identifier extracted from the entrie's raw files if such an
        identifier is supported by the underlying code
        ''',
        categories=[MongoMetadata, UserProvidableMetadata],
        a_search=Search(many_or='append'))

    domain = metainfo.Quantity(
        type=metainfo.MEnum('dft', 'ems', 'qcms'),
        description='The material science domain',
        categories=[MongoMetadata, UserProvidableMetadata],
        a_search=Search())

    published = metainfo.Quantity(
        type=bool, default=False,
        description='Indicates if the entry is published',
        categories=[MongoMetadata, OasisMetadata],
        a_search=Search())

    processed = metainfo.Quantity(
        type=bool, default=False,
        description='Indicates that the entry is successfully processed.',
        categories=[MongoMetadata],
        a_search=Search())

    last_processing = metainfo.Quantity(
        type=metainfo.Datetime,
        description='The datetime of the last processing',
        categories=[MongoMetadata],
        a_search=Search())

    processing_errors = metainfo.Quantity(
        type=str, shape=['*'], description='Errors that occured during processing',
        a_search=Search(many_and='append'))

    nomad_version = metainfo.Quantity(
        type=str,
        description='The NOMAD version used for the last processing',
        categories=[MongoMetadata],
        a_search=Search(many_or='append'))
    nomad_commit = metainfo.Quantity(
        type=str,
        description='The NOMAD commit used for the last processing',
        categories=[MongoMetadata],
        a_search=Search(many_or='append'))
    parser_name = metainfo.Quantity(
        type=str,
        description='The NOMAD parser used for the last processing',
        a_search=Search(many_or='append'))

    comment = metainfo.Quantity(
        type=str, categories=[MongoMetadata, EditableUserMetadata],
        description='A user provided comment for this entry',
        a_search=Search(mapping=Text()))

    references = metainfo.Quantity(
        type=str, shape=['0..*'], categories=[MongoMetadata, EditableUserMetadata],
        description='User provided references (URLs) for this entry',
        a_search=Search())

    external_db = metainfo.Quantity(
        type=metainfo.MEnum('EELSDB', 'Materials Project', 'AFLOW', 'OQMD'), categories=[MongoMetadata, UserProvidableMetadata],
        description='The repository or external database where the original data resides',
        a_search=Search())

    uploader = metainfo.Quantity(
        type=user_reference, categories=[MongoMetadata],
        description='The uploader of the entry',
        a_flask=dict(admin_only=True, verify=User),
        a_search=[
            Search(
                description='The full name of the authors for exact searches',
                metric_name='uploaders', metric='cardinality',
                many_or='append', search_field='uploader.name.keyword',
                statistic_size=10,
                statistic_order='_count'),
            Search(
                name='uploader_id', search_field='uploader.user_id',
                description='The full name of the authors',)
        ])

    origin = metainfo.Quantity(
        type=str,
        description='''
            A short human readable description of the entries origin. Usually it is the
            handle of an external database/repository or the name of the uploader.
        ''',
        derived=derive_origin,
        a_search=Search(statistic_size=10, statistic_order='_count'))

    coauthors = metainfo.Quantity(
        type=author_reference, shape=['0..*'], default=[], categories=[MongoMetadata, EditableUserMetadata],
        description='A user provided list of co-authors',
        a_flask=dict(verify=User))

    authors = metainfo.Quantity(
        type=author_reference, shape=['0..*'],
        description='All authors (uploader and co-authors)',
        derived=derive_authors,
        a_search=Search(
            description='The full name of the authors for exact searches',
            metric='cardinality',
            many_or='append', search_field='authors.name.keyword'))

    shared_with = metainfo.Quantity(
        type=user_reference, shape=['0..*'], default=[], categories=[MongoMetadata, EditableUserMetadata],
        description='A user provided list of userts to share the entry with',
        a_flask=dict(verify=User))

    owners = metainfo.Quantity(
        type=user_reference, shape=['0..*'],
        description='All owner (uploader and shared with users)',
        derived=lambda entry: ([entry.uploader] if entry.uploader is not None else []) + entry.shared_with,
        a_search=Search(
            description='The full name of the owners for exact searches',
            many_or='append', search_field='owners.name.keyword'))

    license = metainfo.Quantity(
        type=str,
        description='''
            A short license description (e.g. CC BY 4.0), that refers to the
            license of this entry.
        ''',
        default='CC BY 4.0',
        categories=[MongoMetadata, EditableUserMetadata])

    with_embargo = metainfo.Quantity(
        type=bool, default=False, categories=[MongoMetadata, EditableUserMetadata],
        description='Indicated if this entry is under an embargo',
        a_search=Search())

    upload_time = metainfo.Quantity(
        type=metainfo.Datetime, categories=[MongoMetadata, OasisMetadata],
        description='The date and time this entry was uploaded to nomad',
        a_flask=dict(admin_only=True),
        a_search=Search(order_default=True))

    upload_name = metainfo.Quantity(
        type=str, categories=[MongoMetadata],
        description='The user provided upload name',
        a_search=Search(many_or='append'))

    datasets = metainfo.Quantity(
        type=dataset_reference, shape=['0..*'], default=[],
        categories=[MongoMetadata, EditableUserMetadata],
        description='A list of user curated datasets this entry belongs to.',
        a_flask=dict(verify=Dataset),
        a_search=[
            Search(
                search_field='datasets.name', many_or='append',
                description='A list of user curated datasets this entry belongs to for exact name search'),
            Search(
                name='dataset_id', search_field='datasets.dataset_id', many_or='append',
                group='datasets_grouped',
                metric='cardinality', metric_name='datasets',
                description='A list of user curated datasets this entry belongs to for exact name search')])

    external_id = metainfo.Quantity(
        type=str, categories=[MongoMetadata, UserProvidableMetadata],
        description='''
        A user provided external id. Usually the id for an entry in an external database
        where the data was imported from.''',
        a_search=Search(many_or='split'))

    last_edit = metainfo.Quantity(
        type=metainfo.Datetime, categories=[MongoMetadata, OasisMetadata],
        description='The date and time the user metadata was edited last',
        a_search=Search())

    formula = metainfo.Quantity(
        type=str, categories=[DomainMetadata],
        description='A (reduced) chemical formula',
        a_search=Search())

    atoms = metainfo.Quantity(
        type=str, shape=['n_atoms'], default=[], categories=[DomainMetadata],
        description='The atom labels of all atoms of the entry\'s material',
        a_search=Search(
            many_and='append', statistic_size=len(ase.data.chemical_symbols)))

    only_atoms = metainfo.Quantity(
        type=str, categories=[DomainMetadata],
        description='The atom labels concatenated in order-number order',
        derived=lambda entry: _only_atoms(entry.atoms),
        a_search=Search(many_and='append', derived=_only_atoms))

    n_atoms = metainfo.Quantity(
        type=int, categories=[DomainMetadata], default=0,
        description='The number of atoms in the entry\'s material',
        a_search=Search())

    ems = metainfo.SubSection(sub_section=EMSMetadata, a_search=Search())
    dft = metainfo.SubSection(sub_section=DFTMetadata, a_search=Search(), categories=[FastAccess])
    qcms = metainfo.SubSection(sub_section=QCMSMetadata, a_search=Search())
    encyclopedia = metainfo.SubSection(sub_section=EncyclopediaMetadata, categories=[FastAccess], a_search=Search())

    def apply_user_metadata(self, metadata: dict):
        ''' Applies a user provided metadata dict to this calc. '''
        self.m_update(**metadata)

    def apply_domain_metadata(self, archive):
        ''' Used to apply metadata that is related to the domain. '''
        assert self.domain is not None, 'all entries must have a domain'
        domain_sub_section_def = self.m_def.all_sub_sections.get(self.domain)
        domain_section_def = domain_sub_section_def.sub_section
        assert domain_section_def is not None, 'unknown domain %s' % self.domain

        # add domain section if not already there
        domain_section = self.m_get_sub_section(domain_sub_section_def, -1)
        if domain_section is None:
            domain_section = self.m_create(domain_section_def.section_cls)

        domain_section.apply_domain_metadata(archive)


class EntryArchive(metainfo.MSection):

    section_run = metainfo.SubSection(sub_section=Run, repeats=True)
    section_experiment = metainfo.SubSection(sub_section=Experiment)
    section_quantum_cms = metainfo.SubSection(sub_section=QuantumCMS)
    section_workflow = metainfo.SubSection(sub_section=Workflow, categories=[FastAccess])
    section_metadata = metainfo.SubSection(sub_section=EntryMetadata, categories=[FastAccess])

    processing_logs = metainfo.Quantity(
        type=Any, shape=['0..*'],
        description='The processing logs for this entry as a list of structlog entries.')


# preemptively create the elasticsearch document definition, which populates metrics and
# search quantities in the search_extension
EntryMetadata.m_def.a_elastic.document
