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

"""
This module contains all functions necessary to manage DOI via datacite.org and its
MDS API (https://support.datacite.org/docs/mds-api-guide).
"""

import datetime

from mongoengine import DateTimeField, Document, EmbeddedDocument, StringField

from nomad.datacite import (
    create_attributes_from_args,
    create_doi,
    delete_doi,
    generate_target_url,
    generate_unique_doi_name,
    publish_doi,
)
from nomad.datamodel import User


class DOIException(Exception):
    """DOI-related errors, including errors with DataCite."""

    pass


# Collection name d_o_i (auto-generated)
# Only used for Dataset
class DOI(Document):
    doi = StringField(primary_key=True)
    url = StringField()
    metadata_url = StringField()  # unnecessary
    doi_url = StringField()  # unnecessary
    state = StringField()
    create_time = DateTimeField()
    metadata_xml = StringField()  # unnecessary

    @staticmethod
    def create() -> 'DOI':
        """Creates a unique DOI with the NOMAD DOI prefix for a dataset or upload."""
        doi = DOI()
        doi.create_time = datetime.datetime.now(datetime.timezone.utc)
        doi.doi = generate_unique_doi_name()
        doi.state = 'created'
        doi.url = generate_target_url(doi.doi, 'dataset')
        doi.save()

        return doi

    def create_draft(self, title: str, publicationYear: int, user: User):
        attributes = create_attributes_from_args(title, publicationYear, user)
        attributes.doi = self.doi
        attributes.url = self.url
        create_doi(attributes)

        self.state = 'draft'
        self.save()

    def delete(self, *args, **kwargs):
        """Deletes the DOI. Only allowed for drafts."""
        assert self.state == 'draft', 'can only delete drafts'
        delete_doi(self.doi)

        super().delete(*args, **kwargs)

    def make_findable(self):
        """Makes the DOI findable in DataCite.

        After this, the DOI cannot be modified or deleted anymore.
        """
        assert self.state == 'draft', 'can only make drafts findable'
        publish_doi(self.doi)

        self.state = 'findable'
        self.save()


class EmbeddedDOI(EmbeddedDocument):
    meta = {'strict': False}

    id = StringField(primary_key=True)
