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
This duplicates functionality for .base.py. It represents first pieces of a transition
towards using the new metainfo system for all repository metadata.
"""
from typing import Dict
from cachetools import cached, TTLCache
from elasticsearch_dsl import Keyword

from nomad import metainfo, config
import nomad.metainfo.mongoengine


class User(metainfo.MSection):
    """ A NOMAD user.

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
    """

    user_id = metainfo.Quantity(type=str, a_me=dict(primary_key=True))
    name = metainfo.Quantity(
        type=str,
        derived=lambda user: ('%s %s' % (user.first_name, user.last_name)).strip())
    first_name = metainfo.Quantity(type=str)
    last_name = metainfo.Quantity(type=str)
    email = metainfo.Quantity(
        type=str, a_me=dict(index=True), a_elastic=dict(mapping=Keyword))
    username = metainfo.Quantity(type=str)
    affiliation = metainfo.Quantity(type=str)
    affiliation_address = metainfo.Quantity(type=str)
    created = metainfo.Quantity(type=metainfo.Datetime)
    repo_user_id = metainfo.Quantity(type=str)
    is_admin = metainfo.Quantity(
        type=bool, derived=lambda user: user.user_id == config.services.admin_user_id)

    @staticmethod
    @cached(cache=TTLCache(maxsize=2048, ttl=24 * 3600))
    def get(*args, **kwargs) -> 'User':
        from nomad import infrastructure
        return infrastructure.keycloak.get_user(*args, **kwargs)  # type: ignore

    @staticmethod
    @cached(cache=TTLCache(maxsize=1, ttl=24 * 3600))
    def repo_users() -> Dict[str, 'User']:
        from nomad import infrastructure
        return {
            str(user.repo_user_id): user
            for user in infrastructure.keycloak.search_user()
            if user.repo_user_id is not None
        }


class Dataset(metainfo.MSection):
    """ A Dataset is attached to one or many entries to form a set of data.

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
        legacy_id: The original NOMAD CoE Repository dataset PID. Old DOIs still reference
            datasets based on this id. Is not used for new datasets.
    """
    dataset_id = metainfo.Quantity(
        type=str,
        a_me=dict(primary_key=True))
    name = metainfo.Quantity(
        type=str,
        a_me=dict(index=True))
    user_id = metainfo.Quantity(
        type=str,
        a_me=dict(index=True))
    doi = metainfo.Quantity(
        type=str,
        a_me=dict(index=True))
    legacy_id = metainfo.Quantity(
        type=str,
        a_me=dict(index=True))


class UserMetadata(metainfo.MSection):
    """ NOMAD entry quantities that are given by the user or determined by user actions.

    Args:
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
    """

    comment = metainfo.Quantity(type=str)
    references = metainfo.Quantity(type=str, shape=['0..*'])
    uploader = metainfo.Quantity(type=str, a_flask=dict(admin_only=True, verify=User))
    coauthors = metainfo.Quantity(type=str, shape=['0..*'], a_flask=dict(verify=User))
    shared_with = metainfo.Quantity(type=str, shape=['0..*'], a_flask=dict(verify=User))
    with_embargo = metainfo.Quantity(type=bool)
    upload_time = metainfo.Quantity(type=metainfo.Datetime, a_flask=dict(admin_only=True))
    datasets = metainfo.Quantity(type=str, shape=['0..*'], a_flask=dict(verify=Dataset))


nomad.metainfo.mongoengine.init_section(User)
nomad.metainfo.mongoengine.init_section(Dataset)
