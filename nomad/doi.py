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

'''
This module contains all functions necessary to manage DOI via datacite.org and its
MDS API (https://support.datacite.org/docs/mds-api-guide).
'''
import xml.etree.ElementTree as ET
import datetime
import requests
from requests.auth import HTTPBasicAuth
from mongoengine import Document, StringField, DateTimeField, BinaryField
from mongoengine.errors import NotUniqueError

from nomad.datamodel import User
from nomad import config, utils


def edit_url(doi: str, url: str = None):
    ''' Changes the URL of an already findable DOI. '''
    if url is None:
        url = 'https://nomad-lab.eu/prod/rae/gui/datasets/doi/%s' % doi

    metadata_url = '%s/doi/%s' % (config.datacite.mds_host, doi)
    response = requests.put(
        metadata_url,
        headers={'Content-Type': 'text/plain;charset=UTF-8'},
        data='doi=%s\nurl=%s' % (doi, url), **_requests_args())

    if response.status_code >= 300:
        raise Exception('Unexpected datacite response (status code %d): %s' % (
            response.status_code, response.text))


def _xml(parent, element: str, value: str = None):
    path = element.split('/')
    el = parent
    for segment in path:
        el = ET.SubElement(el, segment)

    if value is not None:
        el.text = value

    return el


def _requests_args():
    return dict(auth=HTTPBasicAuth(config.datacite.user, config.datacite.password))


class DOI(Document):
    doi = StringField(primary_key=True)
    url = StringField()
    metadata_url = StringField()
    doi_url = StringField()
    state = StringField()
    create_time = DateTimeField()
    matadata_xml = BinaryField()

    @staticmethod
    def create(title: str, user: User) -> 'DOI':
        ''' Creates a unique DOI with the NOMAD DOI prefix. '''
        # TODO We use a collection of all DOIs in mongo to ensure uniqueness. We attempt
        # to create new DOIs based on a counter per day until we find a non existing DOI.
        # This might be bad if many DOIs per day are to be expected.
        counter = 1
        create_time = datetime.datetime.now()

        while True:
            doi_str = '%s/NOMAD/%s-%d' % (
                config.datacite.prefix, create_time.strftime('%Y.%m.%d'), counter)

            try:
                doi = DOI(doi=doi_str)
                doi.save(force_insert=True)
                break
            except NotUniqueError:
                counter += 1

        doi.metadata_url = '%s/metadata/%s' % (config.datacite.mds_host, doi_str)
        doi.doi_url = '%s/doi/%s' % (config.datacite.mds_host, doi_str)
        doi.state = 'created'
        doi.create_time = create_time
        doi.url = '%s/dataset/doi/%s' % (config.gui_url(), doi_str)

        affiliation = ''
        if user.affiliation is not None:
            affiliation += user.affiliation.strip()
        if user.affiliation_address is not None:
            affiliation += '; ' + user.affiliation_address.strip()

        if title is None or title.strip() == '':
            title = 'NOMAD Repository Dataset'

        mds_resource = ET.Element("resource")
        mds_resource.attrib['xsi:schemaLocation'] = 'http://datacite.org/schema/kernel-3 http://schema.datacite.org/meta/kernel-3.1/metadata.xsd'
        mds_resource.attrib['xmlns'] = 'http://datacite.org/schema/kernel-3'
        mds_resource.attrib['xmlns:xsi'] = 'http://www.w3.org/2001/XMLSchema-instance'

        mds_identifier = ET.SubElement(mds_resource, 'identifier')
        mds_identifier.text = doi_str
        mds_identifier.attrib['identifierType'] = 'DOI'

        mds_creator = _xml(mds_resource, 'creators/creator')
        _xml(mds_creator, 'creatorName', user.name)
        if affiliation != '':
            _xml(mds_creator, 'affiliation', affiliation)
        _xml(mds_resource, 'titles/title', title.strip())
        _xml(mds_resource, 'publisher', 'NOMAD Repository')
        _xml(mds_resource, 'publicationYear', str(datetime.datetime.now().year))

        doi.metadata_xml = ET.tostring(mds_resource, encoding='UTF-8', method='xml')
        doi.save()

        return doi

    def __handle_datacite_errors(self, response, msg: str):
        if response.status_code >= 300:
            utils.get_logger(__name__).error(
                'could not %s' % msg,
                status_code=response.status_code, body=response.content,
                doi=self.doi)

            return False
        else:
            return True

    def create_draft(self):
        if config.datacite.enabled:
            assert self.state == 'created', 'can only create a draft for created DOIs'
            response = requests.post(
                self.metadata_url,
                headers={'Content-Type': 'application/xml;charset=UTF-8'},
                data=self.metadata_xml, **_requests_args())

            if self.__handle_datacite_errors(response, 'create draft DOI'):
                self.state = 'draft'
                self.save()

    def delete(self, *args, **kwargs):
        if config.datacite.enabled:
            assert self.state == 'draft', 'can only delete drafts'
            response = requests.delete(self.metadata_url, **_requests_args())

            self.__handle_datacite_errors(response, 'delete draft DOI')

        super().delete(*args, **kwargs)

    def make_findable(self):
        if config.datacite.enabled:
            assert self.state == 'draft', 'can only make drafts findable'
            body = ('doi=%s\nurl=%s' % (self.doi, self.url)).encode('utf-8')
            response = requests.put(
                self.doi_url, **_requests_args(),
                headers={'Content-Type': 'text/plain;charset=UTF-8'}, data=body)

            if self.__handle_datacite_errors(response, 'make DOI findable'):
                self.state = 'findable'
                self.save()
