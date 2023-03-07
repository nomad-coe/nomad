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
import json

import yaml
import requests
import re
from urllib.parse import urlparse, parse_qs
from lxml.html.clean import clean_html  # pylint: disable=no-name-in-module

from nomad.metainfo import MSection, Section, Quantity, SubSection, Package, Datetime, MEnum, JSON, Reference
from nomad.datamodel.data import EntryData, ElnIntegrationCategory
from nomad.metainfo.metainfo import SectionProxy

m_package = Package(name='labfolder')


class LabfolderDataElementDataContent(MSection):
    '''The content of a labfolder data grid.'''

    title = Quantity(type=str, description='the title of the table')
    type = Quantity(type=str, description='type')
    value = Quantity(type=str, description='value')
    unit = Quantity(type=str, description='unit')
    physical_quantity_id = Quantity(type=str, description='physical_quantity_id')
    description = Quantity(type=str, description='physical_quantity_id')
    children = SubSection(sub_section=SectionProxy('LabfolderDataElementGrid'), repeats=True)


class LabfolderDataElementGrid(MSection):
    '''A labfolder grid containing data elements.'''

    title = Quantity(type=str, description='the title of the table')
    type = Quantity(type=str, description='the title of the table')
    children = SubSection(sub_section=LabfolderDataElementDataContent, repeats=True)


class LabfolderImportError(Exception):
    pass


class LabfolderElement(MSection):
    m_def = Section(label_quantity='element_type')

    id = Quantity(type=str, description='the stable pointer to the element')
    entry_id = Quantity(type=str, description='the id of the stable pointer to the entry')
    version_id = Quantity(type=str, description='the unique id of the element')
    version_date = Quantity(type=Datetime, description='the creation date of the entry element version (same with the creation date on the first version)')
    creation_date = Quantity(type=Datetime, description='the creation date of the entry element (first version)')
    owner_id = Quantity(type=str, description='the id of the original author')
    element_type = Quantity(type=MEnum('TEXT', 'DATA', 'FILE'), description='Denotes that this is a file element. The value is always `FILE`')

    def download_files(self, labfolder_api_method, archive, logger):
        pass

    def post_process(self, labfolder_api_method, archive, logger, res_data={}):
        pass


class LabfolderTextElement(LabfolderElement):
    content = Quantity(
        type=str, description='The text based content of this element',
        a_browser=dict(value_component='HtmlValue'))

    def post_process(self, *args, **kwargs):
        self.content = clean_html(self.content)


class LabfolderFileElement(LabfolderElement):
    file_name = Quantity(type=str, description='The name of the file')
    file_size = Quantity(type=int, description='The size of the file in bytes')
    content_type = Quantity(type=str, description='The type of the binary content which is sent on header parameter `Content-Type`')

    file = Quantity(type=str, a_browser=dict(adaptor='RawFileAdaptor'))

    def download_files(self, labfolder_api_method, archive, logger):
        response = labfolder_api_method(requests.get, f'/elements/file/{self.id}/download')
        try:
            with archive.m_context.raw_file(self.file_name, 'wb') as f:
                f.write(response.content)
        except Exception as e:
            logger.error('could not download file', exc_info=e, data=dict(file_name=self.file_name))

        self.file = self.file_name

    def post_process(self, labfolder_api_method, archive, logger, res_data={}):
        self.download_files(labfolder_api_method, archive, logger)


class LabfolderImageElement(LabfolderElement):
    title = Quantity(type=str, description='the title of the image element')
    file_size = Quantity(type=int, description='the size of the image file in bytes')
    preview_height = Quantity(type=int, description='height of the downscaled image version, in px')
    preview_width = Quantity(type=int, description='width of the downscaled image version, in px')
    preview_zoom = Quantity(type=float, description='image zoom in the ELN UI, in percentage')
    original_file_content_type = Quantity(type=str, description='the content type of the original uploaded image file')
    annotation_layer_svg = Quantity(type=str, description='The vector graphic used for the image annotation layer, defined in SVG format')

    original_image_file = Quantity(type=str, a_browser=dict(adaptor='RawFileAdaptor'))
    preview_image_file = Quantity(type=str, a_browser=dict(adaptor='RawFileAdaptor'))

    def download_files(self, labfolder_api_method, archive, logger):
        def download(path, file_quantity):
            response = labfolder_api_method(requests.get, f'/elements/image/{self.id}/{path}')

            content_disposition = response.headers.get('Content-Disposition', '')
            match = re.match(r'^attachment; filename="(.+)"$', content_disposition)
            if match:
                file_name = match.group(1)
            else:
                file_name = self.id
                logger.warn('there is no filename for an image', data=dict(element_id=self.id))
            try:
                with archive.m_context.raw_file(file_name, 'wb') as f:
                    f.write(response.content)
            except Exception as e:
                logger.error('could not download file', exc_info=e, data=dict(file_name=self.file_name))

            self.m_set(file_quantity, file_name)

        download('original-data', LabfolderImageElement.original_image_file)
        download('preview-data', LabfolderImageElement.preview_image_file)

    def post_process(self, labfolder_api_method, archive, logger, res_data={}):
        self.download_files(labfolder_api_method, archive, logger)


class LabfolderTableElement(LabfolderElement):
    title = Quantity(type=str, description='the title of the table')
    content = Quantity(
        type=JSON, description='The JSON content of the table element',
        a_browser=dict(value_component='JsonValue'))


class LabfolderDataElement(LabfolderElement):

    data_elements = SubSection(section=LabfolderDataElementGrid, repeats=True)
    labfolder_data = Quantity(
        type=JSON, description='The JSON content of the table element',
        a_browser=dict(value_component='JsonValue'))

    nomad_data = SubSection(sub_section=LabfolderElement)

    nomad_data_schema = Quantity(
        type=Reference(Section),
        a_eln=dict(component='ReferenceEditQuantity')
    )

    def parse_data(self, data_from_response, data_converted):
        for item in data_from_response:
            children = item.get('children', None)
            if children is not None:
                item_dict_name = item.get('title', None)
                data_converted[item_dict_name] = {}
                self.parse_data(children, data_converted[item_dict_name])
            else:
                child_dict = {}
                child_dict_name = item.get('title', None)
                child_dict[child_dict_name] = {}
                child_dict[child_dict_name].update({'value': item.get('value', None)})
                child_dict[child_dict_name].update({'unit': item.get('unit', None)})
                child_dict[child_dict_name].update({'description': item.get('description', None)})
                data_converted.update(child_dict)

    def post_process(self, labfolder_api_method, archive, logger, res_data={}):
        data_from_response = res_data.get('data_elements', None)
        if data_from_response is not None:
            data_converted = {}
            self.parse_data(data_from_response, data_converted)
            self.labfolder_data = data_converted
        else:
            logger.warning('the labfolder api returned no data')


class LabfolderWellPlateElement(LabfolderElement):
    title = Quantity(type=str, description='The title of the well plate template')
    content = Quantity(
        type=JSON, description='The title of the well plate template',
        a_browser=dict(value_component='JsonValue'))
    meta_data = Quantity(
        type=JSON, description='JSON meta data for visualization processing, used to store information about layer colors and well identifiers',
        a_browser=dict(value_component='JsonValue'))


_element_type_path_mapping = {
    'TEXT': 'text',
    'FILE': 'file',
    'IMAGE': 'image',
    'DATA': 'data',
    'TABLE': 'table',
    'WELL_PLATE': 'well-plate'
}

_element_type_section_mapping = {
    'TEXT': LabfolderTextElement,
    'FILE': LabfolderFileElement,
    'IMAGE': LabfolderImageElement,
    'DATA': LabfolderDataElement,
    'TABLE': LabfolderTableElement,
    'WELL_PLATE': LabfolderWellPlateElement
}


class LabfolderProject(EntryData):
    m_def = Section(label='Labfolder Project Import', categories=[ElnIntegrationCategory])

    def __init__(self, *args, **kwargs):
        super(LabfolderProject, self).__init__(*args, **kwargs)

        self.__headers = None
        self.logger = None

    project_url = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity'))
    labfolder_email = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity'))
    password = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity', props=dict(type='password')))
    resync_labfolder_repository = Quantity(
        type=bool,
        a_eln=dict(component='BoolEditQuantity'))

    id = Quantity(type=str)
    version_id = Quantity(type=str)
    author_id = Quantity(type=str)
    project_id = Quantity(type=str)
    version_date = Quantity(type=Datetime)
    creation_date = Quantity(type=Datetime)
    custom_dates = Quantity(type=Datetime, shape=['*'])
    tags = Quantity(type=str, shape=['*'])
    title = Quantity(type=str)
    hidden = Quantity(type=bool)
    editable = Quantity(type=bool)

    elements = SubSection(sub_section=LabfolderElement, repeats=True)

    def _labfolder_api_method(self, method, url, msg='cannot do labfolder api request', **kwargs):
        response = method(
            f'{self._api_base_url}{url}',
            headers=self._headers, **kwargs)

        if response.status_code >= 400:
            self.logger.error(
                msg,
                data=dict(status_code=response.status_code, text=response.text))
            raise LabfolderImportError()

        return {} if url.endswith('/logout') else response

    @property
    def _api_base_url(self):
        match = re.match(r'^(.+)/eln/notebook.*$', self.project_url)
        if not match:
            self.logger.error('unexpected labfolder url format', data=dict(project_url=self.project_url))
            raise LabfolderImportError()

        return f'{match.group(1)}/api/v2'

    @property
    def _headers(self):
        if not self.__headers:
            response = requests.post(
                f'{self._api_base_url}/auth/login',
                json=dict(user=self.labfolder_email, password=self.password))

            if response.status_code != 200:
                self._clear_user_data()
                self.logger.error(
                    'cannot login',
                    data=dict(status_code=response.status_code, text=response.text))
                raise LabfolderImportError()

            self.__headers = dict(Authorization=f'Token {response.json()["token"]}')

        return self.__headers

    def _clear_user_data(self):
        self.labfolder_email = None
        self.password = None

        archive = self.m_root()
        with archive.m_context.raw_file(archive.metadata.mainfile, 'wt') as f:
            if archive.metadata.mainfile.endswith('json'):
                json.dump(dict(data=archive.data.m_to_dict()), f)
            else:
                yaml.dump(dict(data=archive.data.m_to_dict()), f)

    def normalize(self, archive, logger):
        super(LabfolderProject, self).normalize(archive, logger)
        self.logger = logger

        if not self.elements:
            self.resync_labfolder_repository = True

        if self.resync_labfolder_repository:
            if not self.project_url or not self.labfolder_email or not self.password:
                logger.error('missing information, cannot import project')
                raise LabfolderImportError()

            try:
                project_ids = parse_qs(urlparse(self.project_url).fragment[1:])['projectIds']
            except KeyError as e:
                logger.error('cannot parse project ids from url', exc_info=e)
                raise LabfolderImportError()

            data = self._labfolder_api_method(
                requests.get, f'/entries?project_ids={",".join(project_ids)}'
            ).json()

            elements = data[0]['elements']
            del data[0]['elements']

            try:
                self.m_update_from_dict(data[0])
            except Exception as e:
                logger.error('cannot update archive with labfolder data', exc_info=e)
                raise LabfolderImportError()

            # remove potential old content
            self.elements.clear()

            for element in elements:
                element_type = element['type']

                if element_type not in _element_type_path_mapping:
                    logger.warn('unknown element type', data=dict(element_type=element_type))
                    continue

                data = self._labfolder_api_method(
                    requests.get,
                    f'/elements/{_element_type_path_mapping[element_type]}/{element["id"]}/version/{element["version_id"]}'
                ).json()
                nomad_element = _element_type_section_mapping[element_type]()

                nomad_element.m_update_from_dict(data)
                nomad_element.post_process(self._labfolder_api_method, archive, logger, res_data=data)
                self.elements.append(nomad_element)

            # Resetting Token and Logging out: Invalidating all access tokens
            self.resync_labfolder_repository = False
            self._clear_user_data()
            self._labfolder_api_method(requests.post, '/auth/logout')
            logger.info('reached the end')

        elif not self.resync_labfolder_repository and len(self.elements) > 0:
            for element in self.elements:
                if isinstance(element, LabfolderDataElement) and element.labfolder_data and element.nomad_data_schema:
                    try:
                        element.nomad_data = element.nomad_data_schema.m_from_dict(element.labfolder_data)
                    except Exception as e:
                        logger.error('could not apply schema to labfolder data element', exc_info=e)


m_package.init_metainfo()
