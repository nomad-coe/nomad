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
from typing import Union, Iterable
import json
import re
import copy

from nomad.datamodel import EntryArchive, EntryData, Results
from nomad.datamodel.data import ElnIntegrationCategory, ArchiveSection
from nomad.datamodel.metainfo.annotations import ELNAnnotation
from nomad.metainfo import (
    Package,
    Quantity,
    JSON,
    MEnum,
    MSection,
    Datetime,
    SubSection,
    Section,
)
from nomad import utils
from nomad.metainfo.util import camel_case_to_snake_case
from nomad.parsing.parser import MatchingParser

m_package = Package(name='elabftw')


def _remove_at_sign_from_keys(obj):
    obj = copy.deepcopy(obj)

    for k, v in list(obj.items()):
        if k.startswith('@'):
            obj[k.lstrip('@')] = v
            del obj[k]
            k = k.lstrip('@')
        if isinstance(v, dict):
            obj[k] = _remove_at_sign_from_keys(v)
        if isinstance(v, list):
            for i, item in enumerate(v):
                if isinstance(item, dict):
                    obj[k][i] = _remove_at_sign_from_keys(item)

    return obj


def _map_response_to_dict(data: list) -> dict:
    mapped_dict: dict = {}
    cleaned_data = _remove_at_sign_from_keys(data)
    for item in cleaned_data['graph']:
        id = item['id']
        if id not in mapped_dict:
            mapped_dict[id] = item
        else:
            num = int(id.split('__internal_')[1])
            id_internal = f'{id}__internal_{num + 1}'
            mapped_dict[id_internal] = item
    return mapped_dict


def _create_file_section(file, graph, parent_folder_raw_path, logger=None):
    try:
        section = _element_type_section_mapping[file['type']]()
    except Exception:
        logger.error(f"Could not find type fo the file {file['id']}")
        raise ELabFTWParserError(f"Could not find type fo the file {file['id']}")
    section.m_update_from_dict(graph[file['id']])
    try:
        file_name = file['id'].split('./')[1]
        full_path = os.path.join(parent_folder_raw_path, file_name)
        section.post_process(file_name=full_path)
    except Exception:
        logger.error(f"Could not set the file path for file {file['id']}")
    return section


class ELabFTWRef(MSection):
    """Represents a referenced item in ELabFTW entry."""

    row_refs = Quantity(
        type=ArchiveSection,
        a_eln=ELNAnnotation(
            component='ReferenceEditQuantity',
            label='ELabFTW Reference',
        ),
        description='References that connect to each ELabFTW ref. Each item is stored in it individual entry.',
    )


class ELabFTWParserError(Exception):
    """eln parser related errors."""

    pass


class ELabFTWExperimentLink(MSection):
    """
    This class contains information from other experiments that are linked to this specific experiment.
    the external link can be accessed using the query parameter #id:
    https://demo.elabftw.net/experiments.php?mode=view&id={itemid}
    """

    itemid = Quantity(
        type=str, description='id of the external experiment linked to this experiment'
    )
    title = Quantity(type=str, description='title of the external experiment')
    elabid = Quantity(type=str, description='hashed id')
    category = Quantity(
        type=str, description='Category/status of the external experiment link'
    )


class ELabFTWItemLink(ELabFTWExperimentLink):
    """
    This class holds information of the items related to this specific experiment.
    The external link can be accessed via setting the query parameter #related:
    https://demo.elabftw.net/database.php?mode=show&related={itemid}
    """

    bookable = Quantity(type=bool)


class ELabFTWSteps(MSection):
    """
    Steps recorded for the current experiment
    """

    id = Quantity(type=str, description='id of the current step')
    item_id = Quantity(type=str, description='item_id of the current experiment')
    body = Quantity(type=str, description='title of the step')
    ordering = Quantity(
        type=str, description='location of the current step in the overall order'
    )
    finished = Quantity(
        type=bool, description='a boolean if the step is taken/finished'
    )
    finished_time = Quantity(
        type=Datetime, description='time at which the step is finished'
    )
    deadline = Quantity(type=Datetime, description='deadline time')


class ELabFTWExperimentData(MSection):
    """
    Detailed information of the given ELabFTW experiment, such as links to external resources and extra fields, are
    stored here.
    """

    body = Quantity(
        type=str,
        description='an html-tagged string containing the information of this experiment',
        a_browser=dict(render_value='HtmlValue'),
    )
    created_at = Quantity(
        type=Datetime,
        description='Date and time of when this experiment is created at.',
    )
    sharelink = Quantity(
        type=str,
        a_eln=dict(component='URLEditQuantity'),
        description='URL link to this experiment in the ELabFTW repository',
    )
    extra_fields = Quantity(
        type=JSON,
        description='data in the extra_fields field',
        a_browser=dict(value_component='JsonValue'),
    )
    firstname = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity'),
        description="Author's first name",
    )
    fullname = Quantity(
        type=str,
        a_eln=dict(component='StringEditQuantity'),
        description="Author's full name",
    )

    items_links = SubSection(sub_section=ELabFTWItemLink, repeats=True)
    experiments_links = SubSection(sub_section=ELabFTWExperimentLink, repeats=True)
    steps = SubSection(sub_section=ELabFTWSteps, repeats=True)
    references = SubSection(sub_section=ELabFTWRef, repeats=True)

    def normalize(self, archive, logger) -> None:
        exp_ids = [('experiments', exp.itemid) for exp in self.experiments_links]
        res_ids = [('database', exp.itemid) for exp in self.items_links]

        for item in exp_ids + res_ids:
            from nomad.search import search, MetadataPagination

            query = {'external_id': item[1]}
            search_result = search(
                owner='all',
                query=query,
                pagination=MetadataPagination(page_size=1),
                user_id=archive.metadata.main_author.user_id,
            )
            if search_result.pagination.total > 0:
                entry_id = search_result.data[0]['entry_id']
                upload_id = search_result.data[0]['upload_id']
                ref = ELabFTWRef()
                ref.row_refs = f'../uploads/{upload_id}/archive/{entry_id}#data'
                self.references.append(ref)
                if search_result.pagination.total > 1:
                    logger.warn(
                        f'Found {search_result.pagination.total} entries with external id: '
                        f'"{item[1]}". Will use the first one found.'
                    )
            else:
                logger.warn(f'Found no entries with metadata.external_id: "{item[1]}".')


class ELabFTWComment(MSection):
    """
    A section containing comments made on the experiment. It contains a user object that refers to the id of the
    comment creator
    """

    date_created = Quantity(type=Datetime, description='Creation date of the comment')
    text = Quantity(type=str, description="Comment's content")
    author = Quantity(
        type=JSON,
        description='author information',
        a_browser=dict(value_component='JsonValue'),
    )


class ELabFTWBaseSection(MSection):
    """
    General information on the exported files/experiment of the .eln file
    """

    m_def = Section(label_quantity='type')

    id = Quantity(type=str, description='id of the current data-type')
    type = Quantity(type=str, description='type of the data')

    def post_process(self, **kwargs):
        pass


class ELabFTWFile(ELabFTWBaseSection):
    """
    Information of the exported files
    """

    description = Quantity(type=str, description='Description of the file')
    name = Quantity(type=str, description='Name of the file')
    content_size = Quantity(type=str, description='Size of the file')
    content_type = Quantity(type=str, description='Type of this file')
    file = Quantity(type=str, a_browser=dict(adaptor='RawFileAdaptor'))

    def post_process(self, **kwargs):
        file_name = kwargs.get('file_name', None)
        self.file = file_name


class ElabFTWDataset(ELabFTWBaseSection):
    """
    Information of the dataset type. The author information goes here.
    """

    # TODO: the author section here needs to be integrated into NOMAD user authentication (keycloak) system.
    author = Quantity(
        type=JSON,
        description='author information',
        a_browser=dict(value_component='JsonValue'),
    )
    name = Quantity(type=str, description='Name of the Dataset')
    text = Quantity(type=str, description='Body content of the dataset')
    url = Quantity(
        type=str,
        a_eln=dict(component='URLEditQuantity'),
        description='Link to this dataset in ELabFTW repository',
    )
    date_created = Quantity(type=Datetime, description='Creation date')
    date_modified = Quantity(type=Datetime, description='Last modification date')
    keywords = Quantity(
        type=str,
        shape=['*'],
        description='keywords associated with the current experiment',
    )

    comment = SubSection(sub_section=ELabFTWComment, repeats=True)


class ELabFTW(EntryData):
    """
    Each exported .eln formatted file contains ro-crate-metadata.json file which is parsed into this class.
    Important Quantities are:
        id: id of the file which holds metadata info
        date_created: date of when the file is exported

    title is used as an identifier for the GUI to differentiate between the parsed entries and the original file.
    """

    m_def = Section(label='ELabFTW Project Import', categories=[ElnIntegrationCategory])

    id = Quantity(
        type=str,
        description='id of the file containing the metadata information. It should always be ro-crate-metadata.json',
    )
    title = Quantity(type=str, description='Title of the entry')

    date_created = Quantity(type=Datetime, description='Creation date of the .eln')
    sd_publisher = Quantity(
        type=JSON,
        description='Publisher information',
        a_browser=dict(value_component='JsonValue'),
    )

    author = Quantity(type=str, description="Full name of the experiment's author")
    project_id = Quantity(type=str, description='Project ID')
    status = Quantity(
        type=MEnum(
            'Not set', 'Running', 'Waiting', 'Success', 'Need to be redone', 'Fail'
        ),
        description='Status of the Experiment',
    )
    keywords = Quantity(
        type=str, shape=['*'], description='Keywords associated with the Experiment'
    )

    experiment_data = SubSection(sub_section=ELabFTWExperimentData)
    experiment_files = SubSection(sub_section=ELabFTWBaseSection, repeats=True)

    def post_process(self, **kwargs):
        full_name = kwargs.get('full_name')
        self.author = full_name


_element_type_section_mapping = {'File': ELabFTWFile, 'Dataset': ElabFTWDataset}


class ELabFTWParser(MatchingParser):
    creates_children = True

    def __init__(self) -> None:
        super().__init__(
            name='parsers/elabftw_parser',
            code_name='ElabFTW',
            domain=None,
            mainfile_mime_re=r'text/plain|application/json|text/html',
            mainfile_name_re=r'.*ro-crate-metadata.json$',
        )

    def is_mainfile(
        self,
        filename: str,
        mime: str,
        buffer: bytes,
        decoded_buffer: str,
        compression: str = None,
    ) -> Union[bool, Iterable[str]]:
        is_ro_crate = super().is_mainfile(
            filename, mime, buffer, decoded_buffer, compression
        )
        if not is_ro_crate:
            return False
        try:
            with open(filename, 'rt') as f:
                data = json.load(f)
        except Exception:
            return False

        try:
            if any(
                item.get('@type') == 'SoftwareApplication' for item in data['@graph']
            ):
                root_experiment = next(
                    (item for item in data['@graph'] if item.get('@id') == './')
                )
                no_of_experiments = len(root_experiment['hasPart'])
            else:
                no_of_experiments = len(data['@graph'][-2]['hasPart'])
        except (KeyError, IndexError, TypeError):
            return False

        return [str(item) for item in range(0, no_of_experiments)]

    def parse(
        self, mainfile: str, archive: EntryArchive, logger=None, child_archives=None
    ):
        if logger is None:
            logger = utils.get_logger(__name__)

        title_pattern = re.compile(r'^\d{4}-\d{2}-\d{2} - ([a-zA-Z0-9\-]+) - .*$')

        lab_ids: list[tuple[str, str]] = []
        with open(mainfile, 'rt') as f:
            data = json.load(f)

        snake_case_data = camel_case_to_snake_case(data)
        clean_data = _remove_at_sign_from_keys(snake_case_data)
        graph = {item['id']: item for item in clean_data['graph']}
        experiments = graph['./']

        for index, experiment in enumerate(experiments['has_part']):
            exp_id = experiment['id']
            raw_experiment, exp_archive = graph[exp_id], child_archives[str(index)]

            # hook for matching the older .eln files from Elabftw exported files
            if not any(
                item.get('type') == 'SoftwareApplication'
                for item in clean_data['graph']
            ):
                elabftw_experiment = _parse_legacy(
                    graph,
                    raw_experiment,
                    exp_archive,
                    data,
                    mainfile,
                    exp_id,
                    title_pattern,
                    lab_ids,
                    archive,
                    logger,
                )

                if not archive.results:
                    archive.results = Results()
                    archive.results.eln = Results.eln.sub_section.section_cls()
                    archive.results.eln.lab_ids = [str(lab_id[1]) for lab_id in lab_ids]
                    archive.results.eln.tags = [lab_id[0] for lab_id in lab_ids]

            else:
                elabftw_experiment = _parse_latest(
                    graph,
                    raw_experiment,
                    exp_archive,
                    mainfile,
                    exp_id,
                    archive,
                    logger,
                )

            exp_archive.data = elabftw_experiment

        logger.info('eln parsed successfully')


def _parse_legacy(
    graph,
    raw_experiment,
    exp_archive,
    data,
    mainfile,
    exp_id,
    title_pattern,
    lab_ids,
    archive,
    logger,
) -> ELabFTW:
    elabftw_experiment = ELabFTW()
    try:
        del graph['ro-crate-metadata.json']['type']
    except Exception:
        pass
    elabftw_experiment.m_update_from_dict(graph['ro-crate-metadata.json'])
    elabftw_entity_type = _set_experiment_metadata(
        raw_experiment, exp_archive, elabftw_experiment, logger
    )

    try:
        author_full_name = ' '.join(
            [
                data['graph'][-1]['given_name'],
                data['graph'][-1]['family_name'],
            ]
        )
        elabftw_experiment.post_process(full_name=author_full_name)
    except Exception:
        logger.error('Could not extract the author name')

    mainfile_raw_path = os.path.dirname(mainfile)
    parent_folder_raw_path = mainfile.split('/')[-2]

    extracted_title = _set_child_entry_name(exp_id, exp_archive, archive, logger)

    matched = title_pattern.findall(extracted_title)
    if matched:
        title = matched[0]
    else:
        title = extracted_title
    elabftw_experiment.title = title

    path_to_export_json = os.path.join(mainfile_raw_path, exp_id, 'export-elabftw.json')
    try:
        with open(path_to_export_json, 'rt') as f:
            export_data = json.load(f)
    except FileNotFoundError:
        raise ELabFTWParserError(f"Couldn't find export-elabftw.json file.")

    def clean_nones(value):
        if isinstance(value, list):
            return [
                clean_nones(x) for x in value
            ]  # ??? what to do if this list has None values
        if isinstance(value, dict):
            return {k: clean_nones(v) for k, v in value.items() if v is not None}

        return value

    experiment_data = ELabFTWExperimentData()
    try:
        experiment_data.m_update_from_dict(clean_nones(export_data[0]))
    except (IndexError, KeyError, TypeError):
        logger.warning(
            f"Couldn't read and parse the data from export-elabftw.json file"
        )
    try:
        experiment_data.extra_fields = export_data[0]['metadata']['extra_fields']
    except Exception:
        pass
    elabftw_experiment.experiment_data = experiment_data

    try:
        exp_archive.metadata.comment = elabftw_entity_type

        lab_ids.extend(
            ('experiment_link', experiment_link['itemid'])
            for experiment_link in export_data[0]['experiments_links']
        )
        lab_ids.extend(
            ('item_link', experiment_link['itemid'])
            for experiment_link in export_data[0]['items_links']
        )
    except Exception:
        pass
    for file_id in raw_experiment['has_part']:
        file_section = _create_file_section(
            graph[file_id['id']], graph, parent_folder_raw_path, logger
        )
        elabftw_experiment.experiment_files.append(file_section)

    return elabftw_experiment


def _set_child_entry_name(exp_id, child_archive, archive, logger):
    matched_title = exp_id.split('/')
    if len(matched_title) > 1:
        extracted_title = matched_title[1]
        archive.metadata.m_update_from_dict(dict(entry_name='ELabFTW Schema'))
        child_archive.metadata.m_update_from_dict(dict(entry_name=extracted_title))
    else:
        logger.warning(f"Couldn't extract the title from {exp_id}")
        extracted_title = None
    return extracted_title


def _set_experiment_metadata(raw_experiment, exp_archive, elab_instance, logger):
    try:
        exp_external_id = raw_experiment['url'].split('&id=')[1]
        if match := re.search(r'.*/([^/]+)\.php', raw_experiment['url']):
            elabftw_entity_type = match.group(1)
        exp_archive.metadata.external_id = str(exp_external_id)
        elab_instance.project_id = exp_external_id
    except Exception:
        logger.error('Could not set the the external_id from the experiment url')
        elabftw_entity_type = None
    return elabftw_entity_type


def _parse_latest(
    graph,
    raw_experiment,
    exp_archive,
    mainfile,
    exp_id,
    archive,
    logger,
) -> ELabFTW:
    latest_elab_instance = ELabFTW(
        author=raw_experiment.get('author').get('id'),
        title=raw_experiment.get('name', None),
        keywords=raw_experiment.get('keywords', '').split(','),
        id=raw_experiment.get('id', ''),
        status=raw_experiment.get('creative_work_status', 'Not set'),
    )

    _ = _set_experiment_metadata(
        raw_experiment, exp_archive, latest_elab_instance, logger
    )
    data_section = ELabFTWExperimentData(
        body=raw_experiment.get('text', None),
        created_at=raw_experiment.get('date_created', None),
        extra_fields={
            i: value
            for i, value in enumerate(raw_experiment.get('variable_measured', []))
        },
    )
    data_section.steps.extend(
        [
            ELabFTWSteps(
                ordering=step.get('position', None),
                deadline=step.get('expires', None),
                finished=step.get('creative_work_status', None) == 'finished',
                body=step.get('item_list_element', [{'text': None}])[0]['text'],
            )
            for step in raw_experiment.get('step', [])
        ]
    )

    data_section.experiments_links.extend(
        [
            ELabFTWExperimentLink(title=exp_link.get('id', None))
            for exp_link in raw_experiment.get('mentions', [])
        ]
    )

    parent_folder_raw_path = mainfile.split('/')[-2]
    for file_id in raw_experiment.get('has_part', []):
        file_section = _create_file_section(
            graph[file_id['id']], graph, parent_folder_raw_path, logger
        )
        latest_elab_instance.experiment_files.append(file_section)

    latest_elab_instance.m_add_sub_section(
        latest_elab_instance.m_def.all_sub_sections['experiment_data'], data_section
    )

    _ = _set_child_entry_name(exp_id, exp_archive, archive, logger)

    return latest_elab_instance


m_package.__init_metainfo__()
