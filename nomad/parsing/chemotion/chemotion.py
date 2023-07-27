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

import numpy as np
from nomad.datamodel import EntryArchive, EntryData
from nomad.datamodel.data import ElnIntegrationCategory
from nomad.metainfo import Package, Quantity, JSON, MSection, Datetime, Section, SubSection
from nomad import utils
from nomad.parsing.parser import MatchingParser

m_package = Package(name='chemotion')


class ChemotionGeneralMetainfo(MSection):
    user_id = Quantity(type=str)
    created_at = Quantity(type=Datetime)
    updated_at = Quantity(type=Datetime)
    deleted_at = Quantity(type=Datetime)
    molecule_id = Quantity(type=str)
    sample_id = Quantity(type=str)
    collection_id = Quantity(type=str)


class ChemotionReactionSample(ChemotionGeneralMetainfo):
    reaction_id = Quantity(type=str)
    reference = Quantity(type=bool)
    position = Quantity(type=int)
    waste = Quantity(type=bool)
    coefficient = Quantity(type=np.float16)


class ChemotionCollection(ChemotionGeneralMetainfo):
    ancestry = Quantity(type=str)
    label = Quantity(type=str)
    shared_by_id = Quantity(type=str)
    is_shared = Quantity(type=bool)
    permission_level = Quantity(type=int)
    sample_detail_level = Quantity(type=int)
    reaction_detail_level = Quantity(type=int)
    wellplate_detail_level = Quantity(type=int)
    position = Quantity(type=int)
    screen_detail_level = Quantity(type=int)
    is_locked = Quantity(type=bool)
    is_synchronized = Quantity(type=bool)
    researchplan_detail_level = Quantity(type=int)


class ChemotionSample(ChemotionGeneralMetainfo):
    name = Quantity(type=str)
    target_amount_value = Quantity(type=np.float16)
    target_amount_unit = Quantity(type=str)
    description = Quantity(type=str)
    molfile = Quantity(type=str)
    purity = Quantity(type=str)
    solvent = Quantity(type=str)
    impurities = Quantity(type=str)
    location = Quantity(type=str)
    is_top_secret = Quantity(type=bool)
    ancestry = Quantity(type=str)
    external_label = Quantity(type=str)
    created_by = Quantity(type=str)
    short_label = Quantity(type=str)
    real_amount_value = Quantity(type=str)
    real_amount_unit = Quantity(type=str)
    imported_readout = Quantity(type=str)
    deleted_at = Quantity(type=str)
    sample_svg_file = Quantity(type=str)
    identifier = Quantity(type=str)
    density = Quantity(type=np.float16)
    melting_point = Quantity(type=np.float16)
    boiling_point = Quantity(type=np.float16)
    fingerprint_id = Quantity(type=str)
    xref = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    molarity_value = Quantity(type=np.float16)
    molarity_unit = Quantity(type=str)
    molecule_name_id = Quantity(type=str)
    molfile_version = Quantity(type=str)
    stereo = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    file = Quantity(type=str, a_browser=dict(adaptor='RawFileAdaptor'))

    def post_process(self, **kwargs):
        full_path = os.path.join('images', 'samples', self.sample_svg_file)
        self.file = full_path


class ChemotionCollectionsSample(ChemotionGeneralMetainfo):
    collection_id = Quantity(type=str)


class ChemotionFingerprint(ChemotionGeneralMetainfo):
    fp0 = Quantity(type=str)
    fp1 = Quantity(type=str)
    fp2 = Quantity(type=str)
    fp3 = Quantity(type=str)
    fp4 = Quantity(type=str)
    fp5 = Quantity(type=str)
    fp6 = Quantity(type=str)
    fp7 = Quantity(type=str)
    fp8 = Quantity(type=str)
    fp9 = Quantity(type=str)
    fp10 = Quantity(type=str)
    fp11 = Quantity(type=str)
    fp12 = Quantity(type=str)
    fp13 = Quantity(type=str)
    fp14 = Quantity(type=str)
    fp15 = Quantity(type=str)
    num_set_bits = Quantity(type=int)


class ChemotionMolecule(ChemotionGeneralMetainfo):
    inchikey = Quantity(type=str)
    inchistring = Quantity(type=str)
    density = Quantity(type=np.float16)
    molecular_weight = Quantity(type=np.float16)
    molfile = Quantity(type=str)
    melting_point = Quantity(type=str)
    boiling_point = Quantity(type=str)
    sum_formular = Quantity(type=str)
    names = Quantity(type=str)
    iupac_name = Quantity(type=str)
    molecule_svg_file = Quantity(type=str)
    is_partial = Quantity(type=bool)
    exact_molecular_weight = Quantity(type=np.float16)
    cano_smiles = Quantity(type=str)
    cas = Quantity(type=str)
    molfile_version = Quantity(type=str)
    file = Quantity(type=str, a_browser=dict(adaptor='RawFileAdaptor'))

    def post_process(self, **kwargs):
        full_path = os.path.join('images', 'molecules', self.molecule_svg_file)
        self.file = full_path


class ChemotionMoleculeName(ChemotionGeneralMetainfo):
    description = Quantity(type=str)
    name = Quantity(type=str)


class ChemotionResidue(ChemotionGeneralMetainfo):
    residue_type = Quantity(type=str)
    custom_info = Quantity(type=str)


class ChemotionReaction(ChemotionGeneralMetainfo):
    name = Quantity(type=str)
    description = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    timestamp_start = Quantity(type=str)
    timestamp_stop = Quantity(type=str)
    observation = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    purification = Quantity(type=JSON)
    dangerous_products = Quantity(type=JSON)
    tlc_solvents = Quantity(type=str)
    tlc_description = Quantity(type=str)
    rf_value = Quantity(type=str)
    temperature = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    status = Quantity(type=str)
    reaction_svg_file = Quantity(type=str)
    solvent = Quantity(type=str)
    deleted_at = Quantity(type=str)
    short_label = Quantity(type=str)
    created_by = Quantity(type=str)
    role = Quantity(type=str)
    origin = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    rinchi_string = Quantity(type=str)
    rinchi_long_key = Quantity(type=str)
    rinchi_short_key = Quantity(type=str)
    rinchi_web_key = Quantity(type=str)
    duration = Quantity(type=str)
    file = Quantity(type=str, a_browser=dict(adaptor='RawFileAdaptor'))

    def post_process(self, **kwargs):
        full_path = os.path.join('images', 'reactions', self.reaction_svg_file)
        self.file = full_path


class ChemotionCollectionsReaction(ChemotionGeneralMetainfo):

    reaction_id = Quantity(type=str)


class ChemotionReactionsStartingMaterialSample(ChemotionReactionSample):
    pass


class ChemotionReactionsSolventSample(ChemotionReactionSample):
    pass


class ChemotionReactionsPurificationSolventSample(ChemotionReactionSample):
    pass


class ChemotionReactionsReactantSample(ChemotionReactionSample):
    pass


class ChemotionReactionsProductSample(ChemotionReactionSample):
    pass


class ChemotionWellplate(ChemotionGeneralMetainfo):
    name = Quantity(type=str)
    size = Quantity(type=int)
    description = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))


class ChemotionCollectionsWellplate(ChemotionGeneralMetainfo):
    wellplate_id = Quantity(type=str)


class ChemotionWell(ChemotionGeneralMetainfo):
    wellplate_id = Quantity(type=str)
    position_x = Quantity(type=int)
    position_y = Quantity(type=int)
    readout = Quantity(type=Datetime)
    additive = Quantity(type=Datetime)


class ChemotionScreen(ChemotionGeneralMetainfo):
    description = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    name = Quantity(type=str)
    result = Quantity(type=str)
    collaborator = Quantity(type=str)
    conditions = Quantity(type=str)
    requirements = Quantity(type=str)


class ChemotionCollectionsScreen(ChemotionGeneralMetainfo):
    screen_id = Quantity(type=str)


class ChemotionScreensWellplate(ChemotionGeneralMetainfo):
    screen_id = Quantity(type=str)
    wellplate_id = Quantity(type=str)


class ChemotionResearchPlan(ChemotionGeneralMetainfo):
    name = Quantity(type=str)
    body = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    sdf_file = Quantity(type=str)
    svg_file = Quantity(type=str)
    created_by = Quantity(type=str)

    def post_process(self):
        pass


class ChemotionCollectionsResearchPlan(ChemotionGeneralMetainfo):
    research_plan_id = Quantity(type=str)


class ChemotionContainer(MSection):
    ancestry = Quantity(type=str)
    containable_id = Quantity(type=str)
    containable_type = Quantity(type=str)
    name = Quantity(type=str)
    container_type = Quantity(type=str)
    description = Quantity(type=str)
    extended_metadata = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    created_at = Quantity(type=Datetime)
    updated_at = Quantity(type=Datetime)
    parent_id = Quantity(type=str)


class ChemotionAttachment(ChemotionGeneralMetainfo):
    attachable_id = Quantity(type=str)
    filename = Quantity(type=str)
    identifier = Quantity(type=str)
    checksum = Quantity(type=str)
    storage = Quantity(type=str)
    created_by = Quantity(type=str)
    created_for = Quantity(type=str)
    version = Quantity(type=int)
    content_type = Quantity(type=str)
    bucket = Quantity(type=str)
    key = Quantity(type=str)
    thumb = Quantity(type=bool)
    folder = Quantity(type=bool)
    attachable_type = Quantity(type=str)
    aasm_state = Quantity(type=str)
    file = Quantity(type=str, a_browser=dict(adaptor='RawFileAdaptor'))

    def post_process(self, **kwargs):
        full_path = os.path.join('attachments', self.key)
        self.file = full_path


class ChemotionLiteral(ChemotionGeneralMetainfo):
    literature_id = Quantity(type=str)
    element_id = Quantity(type=str)
    element_type = Quantity(type=str)
    category = Quantity(type=str)


class ChemotionLiterature(ChemotionGeneralMetainfo):
    title = Quantity(type=str)
    url = Quantity(type=str)
    refs = Quantity(type=JSON, a_browser=dict(value_component='JsonValue'))
    doi = Quantity(type=str)


class Chemotion(EntryData):
    '''
    Each exported .eln formatted file contains ro-crate-metadata.json file which is parsed into this class.
    Important Quantities are:
        id: id of the file which holds metadata info
        date_created: date of when the file is exported

    title is used as an identifier for the GUI to differentiate between the parsed entries and the original file.
    '''
    m_def = Section(label='Chemotion Project Import', categories=[ElnIntegrationCategory])

    name = Quantity(type=str, description='Name of the project', a_eln=dict(component='StringEditQuantity'))
    title = Quantity(type=str, description='Title of the entry', a_eln=dict(component='StringEditQuantity'))

    date_created = Quantity(type=Datetime, description='Creation date of the .eln',
                            a_eln=dict(component='DateTimeEditQuantity'))
    author = Quantity(type=str, description='Full name of the experiment\'s author',
                      a_eln=dict(component='StringEditQuantity'))

    Collection = SubSection(sub_section=ChemotionCollection, repeats=True)
    Sample = SubSection(sub_section=ChemotionSample, repeats=True)
    CollectionsSample = SubSection(sub_section=ChemotionCollectionsSample, repeats=True)
    Fingerprint = SubSection(sub_section=ChemotionFingerprint, repeats=True)
    Molecule = SubSection(sub_section=ChemotionMolecule, repeats=True)
    MoleculeName = SubSection(sub_section=ChemotionMoleculeName, repeats=True)
    Container = SubSection(sub_section=ChemotionContainer, repeats=True)
    Attachment = SubSection(sub_section=ChemotionAttachment, repeats=True)
    Reactions = SubSection(sub_section=ChemotionReaction, repeats=True)
    CollectionsReaction = SubSection(sub_section=ChemotionCollectionsReaction, repeats=True)
    ReactionsStartingMaterialSample = SubSection(sub_section=ChemotionReactionsStartingMaterialSample, repeats=True)
    ReactionsSolventSample = SubSection(sub_section=ChemotionReactionsSolventSample, repeats=True)
    ReactionsProductSample = SubSection(sub_section=ChemotionReactionsProductSample, repeats=True)
    ResearchPlan = SubSection(sub_section=ChemotionResearchPlan, repeats=True)
    CollectionsResearchPlan = SubSection(sub_section=ChemotionCollectionsResearchPlan, repeats=True)


_element_type_section_mapping = {
    'Collection': ChemotionCollection,
    'Sample': ChemotionSample,
    'CollectionsSample': ChemotionCollectionsSample,
    'Fingerprint': ChemotionFingerprint,
    'Molecule': ChemotionMolecule,
    'MoleculeName': ChemotionMoleculeName,
    'Container': ChemotionContainer,
    'Attachment': ChemotionAttachment,
    'Reaction': ChemotionReaction,
    'CollectionsReaction': ChemotionCollectionsReaction,
    'ReactionsStartingMaterialSample': ChemotionReactionsStartingMaterialSample,
    'ReactionsSolventSample': ChemotionReactionsSolventSample,
    'ReactionsProductSample': ChemotionReactionsProductSample,
    'ResearchPlan': ChemotionResearchPlan,
    'CollectionsResearchPlan': ChemotionCollectionsResearchPlan
}


def _set_inf_to_nan_if_string(dct, key):
    if key in dct and isinstance(dct[key], str):
        dct[key] = np.NaN


class ChemotionParser(MatchingParser):
    creates_children = True

    def __init__(self) -> None:
        super().__init__(
            name='parsers/chemotion', code_name='Chemotion',
            domain=None,
            mainfile_mime_re=r'application/json|text/plain',
            mainfile_name_re=r'^.*export.json$')

    def is_mainfile(
        self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
        compression: str = None
    ) -> Union[bool, Iterable[str]]:
        is_export_json = super().is_mainfile(filename, mime, buffer, decoded_buffer, compression)
        if not is_export_json:
            return False

        return [str(0)]

    def parse(self, mainfile: str, archive: EntryArchive, logger=None, child_archives=None):

        if logger is None:
            logger = utils.get_logger(__name__)

        chemotion = Chemotion()

        with open(mainfile, 'rt') as f:
            data = json.load(f)

        for item_name, item_content in data.items():
            for sub_item in item_content.values():
                try:
                    chemotion_subsection = _element_type_section_mapping[item_name]()
                    _set_inf_to_nan_if_string(sub_item, 'melting_point')
                    _set_inf_to_nan_if_string(sub_item, 'boiling_point')

                    chemotion_subsection.m_update_from_dict(sub_item)

                    if item_name in ['Sample', 'Molecule', 'Reaction', 'ResearchPlan']:
                        chemotion_subsection.post_process()
                    item_name = 'Reactions' if item_name == 'Reaction' else item_name
                    sub_section_def = getattr(chemotion.m_def.section_cls, item_name)
                    chemotion.m_add_sub_section(sub_section_def, chemotion_subsection, -1)
                except Exception as e:
                    logger.error('No dot (.) is allowed in the column name.', details=dict(column=item_name), exc_info=e)
        for child_archive in child_archives.values():
            child_archive.data = chemotion
        logger.info('eln parsed successfully')


m_package.__init_metainfo__()
