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

'''
This module contains all functions necessary to manage DOI via datacite.org and its
MDS API (https://support.datacite.org/docs/mds-api-guide).
'''
import xml.etree.ElementTree as ET
import datetime
import requests
from requests.auth import HTTPBasicAuth
from mongoengine import Document, StringField, DateTimeField
from mongoengine.errors import NotUniqueError

from nomad.datamodel import User
from nomad import config, utils
from fastapi import HTTPException


class DOIException(Exception):
    ''' Datacite-requests related errors. '''
    pass


def _create_dataset_url(doi: str) -> str:
    '''
    Returns:
        The url that set in the DOI record and is used to resolve the DOI. The URL
        points to the dataset page in the NOMAD GUI.
    '''
    return f'{config.gui_url()}/dataset/doi/{doi}'


def edit_doi_url(doi: str, url: str = None):
    ''' Changes the URL of an already findable DOI. '''
    if url is None:
        url = _create_dataset_url(doi)

    doi_url = '%s/doi/%s' % (config.datacite.mds_host, doi)
    headers = {'Content-Type': 'text/plain;charset=UTF-8'}
    data = f'doi={doi}\nurl={url}'
    response = requests.put(doi_url, headers=headers, data=data, **_requests_args())

    # There seems to be a bug in datacite. Some old records might have somehow invalid
    # xml stored at datacite and for those the DOI update might fail. We try to
    # get the xml string, parse and re-serialize and put it again. After this the
    # url update might work.
    if response.status_code == 422 and 'No matching global declaration available' in response.text:
        metadata_url = f'{config.datacite.mds_host}/metadata/{doi}'
        response = requests.get(metadata_url, **_requests_args())
        original_xml = response.text
        tree = ET.fromstring(original_xml)
        repaired_xml = ET.tostring(tree, encoding='UTF-8', method='xml').decode('utf-8')
        response = requests.put(metadata_url, **_requests_args())
        requests.put(
            metadata_url,
            headers={'Content-Type': 'application/xml;charset=UTF-8'},
            data=repaired_xml.encode('utf-8'), **_requests_args())
        response = requests.put(doi_url, headers=headers, data=data, **_requests_args())
        if response.status_code >= 300:
            raise Exception(f'Encountered known xml problems for {doi}. But could not fix.')

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
    metadata_xml = StringField()

    @staticmethod
    def create(title: str, user: User) -> 'DOI':
        ''' Creates a unique DOI with the NOMAD DOI prefix. '''
        # TODO We use a collection of all DOIs in mongo to ensure uniqueness. We attempt
        # to create new DOIs based on a counter per day until we find a non existing DOI.
        # This might be bad if many DOIs per day are to be expected.
        counter = 1
        create_time = datetime.datetime.utcnow()

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
        doi.url = _create_dataset_url(doi_str)

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

        doi.metadata_xml = ET.tostring(mds_resource, encoding='UTF-8', method='xml').decode('utf-8')
        doi.save()

        return doi

    def __handle_datacite_errors(self, response, msg: str):
        if response is None or response.status_code >= 300:
            utils.get_logger(__name__).error(
                'could not %s' % msg,
                status_code=response.status_code, body=response.content,
                doi=self.doi)

            raise DOIException()

        return True

    def create_draft(self):
        if config.datacite.enabled:
            assert self.state == 'created', 'can only create a draft for created DOIs'
            response = None
            try:
                response = requests.post(
                    self.metadata_url,
                    headers={'Content-Type': 'application/xml;charset=UTF-8'},
                    data=self.metadata_xml.encode('utf-8'), **_requests_args())
            except HTTPException:
                pass

            if self.__handle_datacite_errors(response, 'create draft DOI'):
                self.state = 'draft'
                self.save()

    def delete(self, *args, **kwargs):
        if config.datacite.enabled:
            assert self.state == 'draft', 'can only delete drafts'
            response = None

            try:
                response = requests.delete(self.metadata_url, **_requests_args())
            except HTTPException:
                pass

            self.__handle_datacite_errors(response, 'delete draft DOI')

        super().delete(*args, **kwargs)

    def make_findable(self):
        if config.datacite.enabled:
            assert self.state == 'draft', 'can only make drafts findable'
            body = ('doi=%s\nurl=%s' % (self.doi, self.url)).encode('utf-8')
            response = None

            try:
                response = requests.put(
                    self.doi_url, **_requests_args(),
                    headers={'Content-Type': 'text/plain;charset=UTF-8'}, data=body)
            except HTTPException:
                pass

            if self.__handle_datacite_errors(response, 'make DOI findable'):
                self.state = 'findable'
                self.save()
