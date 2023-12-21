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
import os

from nomad.metainfo import (
    MSection,
    Section,
    Quantity,
    SubSection,
    Package,
    Datetime,
    MEnum,
    JSON,
    Reference,
)
from nomad.datamodel.data import EntryData, ElnIntegrationCategory

m_package = Package(name='openbis')


class OpenbisImportError(Exception):
    pass


class OpenbisBaseSection(MSection):
    code = Quantity(type=str, description='Code for the current section.')
    description = Quantity(type=str, description='Description of the current section.')
    frozen = Quantity(type=bool, description='If current space is frozen.')
    permId = Quantity(type=str, description='Permanent id of the current section.')
    modificationDate = Quantity(type=Datetime)
    registration_date = Quantity(type=Datetime)
    custom_dates = Quantity(type=Datetime, shape=['*'])

    def download_files(self, openbis_element, archive, logger):
        pass

    def post_process(self, openbis_element, archive, logger):
        pass


class OpenbisAttachment(MSection):
    file = Quantity(type=str, a_browser=dict(adaptor='RawFileAdaptor'))


class OpenbisExperiment(OpenbisBaseSection):
    attachments = SubSection(sub_section=OpenbisAttachment, repeats=True)

    def download_files(self, experiment, archive, logger):
        try:
            datasets = experiment.get_datasets()
            for ds in datasets:
                filepath = os.path.join(
                    archive.m_context.upload_files.external_os_path,
                    'raw',
                    experiment.identifier[1:].lower(),
                )
                ds.download(
                    destination=filepath,
                    create_default_folders=False,
                    wait_until_finished=False,
                )
                new_attachment = OpenbisAttachment()
                file_name = '/'.join(
                    [experiment.identifier[1:].lower(), ds.file_list[0].split('/')[1]]
                )
                new_attachment.file = file_name
                self.attachments.append(new_attachment)
        except Exception as e:
            logger.error(
                f'Could not download the file {filepath}.',
                exc_info=e,
                data=dict(filepath=filepath),
            )


class OpenbisProject(OpenbisBaseSection):
    m_def = Section(label_quantity='code')

    identifier = Quantity(
        type=str, description='Path of the current item in the Openbis file system.'
    )

    experiments = SubSection(sub_section=OpenbisExperiment, repeats=True)


class OpenbisSpace(OpenbisBaseSection):
    m_def = Section(
        label_quantity='code', description='Name of the current Openbis Space.'
    )
    projects = SubSection(sub_section=OpenbisProject, repeats=True)


class OpenbisEntry(EntryData):
    m_def = Section(label='Openbis Project Import', categories=[ElnIntegrationCategory])

    def __init__(self, *args, **kwargs):
        super(OpenbisEntry, self).__init__(*args, **kwargs)
        self.logger = None

    project_url = Quantity(type=str, a_eln=dict(component='StringEditQuantity'))
    username = Quantity(type=str, a_eln=dict(component='StringEditQuantity'))
    password = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity', props=dict(type='password')),
    )

    spaces = SubSection(sub_section=OpenbisSpace, repeats=True)

    def _clear_login_info(self, archive):
        self.username = None
        self.password = None
        archive.data.username = None
        archive.data.password = None

    def normalize(self, archive, logger):
        super(OpenbisEntry, self).normalize(archive, logger)
        self.logger = logger

        if not self.project_url or not self.username or not self.password:
            logger.warning(
                'Please make sure all fields of `project_url`, `username` and `password` are filled.'
            )
        else:
            # Initializing pybis object
            try:
                from pybis import Openbis

                openbis = Openbis(self.project_url, verify_certificates=False)
            except Exception:
                self._clear_login_info(archive)
                raise OpenbisImportError(
                    'Failed to connect to the openbis server. Check the project url to be correct.'
                )

            try:
                openbis.login(self.username, self.password, save_token=True)
            except Exception:
                self._clear_login_info(archive)
                raise OpenbisImportError(
                    'Failed to login to the openbis server. Check the username and password to be correct.'
                )

            # remove potential old content
            self.spaces.clear()

            try:
                spaces = openbis.get_spaces()
            except Exception:
                self._clear_login_info(archive)
                raise OpenbisImportError(
                    'Failed to fetch spaces. The pybis package might have been changed.'
                )

            # parsing the content in the spaces
            for space in spaces:
                space_element = OpenbisSpace()
                space_element.m_update_from_dict(space.attrs.all())

                try:
                    projects = space.get_projects()
                except Exception:
                    self._clear_login_info(archive)
                    raise OpenbisImportError(
                        'Failed to fetch projects. The pybis package might have been changed.'
                    )
                for project in projects:
                    project_element = OpenbisProject()
                    project_element.m_update_from_dict(project.attrs.all())
                    space_element.projects.append(project_element)

                    try:
                        experiments = project.get_experiments()
                    except Exception:
                        self._clear_login_info(archive)
                        raise OpenbisImportError(
                            'Failed to fetch experiments. The pybis package might have been changed.'
                        )
                    for experiment in experiments:
                        experiment_element = OpenbisExperiment()
                        experiment_element.m_update_from_dict(experiment.attrs.all())
                        try:
                            experiment_element.download_files(
                                experiment, archive, logger
                            )
                        except Exception as e:
                            logger.error('Failed to download attachments.', exec_info=e)

                        project_element.experiments.append(experiment_element)

                self.spaces.append(space_element)

            openbis.logout()

            self._clear_login_info(archive)

            logger.info('Transferring openbis data is finished.')


m_package.init_metainfo()
