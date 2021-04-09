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

from typing import Any, Tuple, Dict, Union, List

from nomad import infrastructure, config
from nomad.metainfo import MSection, Definition, Quantity, Reference, SubSection, Section
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.common import FastAccess


def create_partial_archive(archive: EntryArchive) -> Dict:
    '''
    Creates a partial archive JSON serializable dict that can be stored directly.
    The given archive is filtered based on the metainfo category ``FastAccess``.
    Selected sections and other data that they reference (recursively) comprise the
    resulting partial archive.

    TODO at the moment is hard coded and NOT informed by the metainfo. We simply
    add sections EntryMetadata and Workflow.

    Arguments:
        archive: The archive as an :class:`EntryArchive` instance.

    Returns: the partial archive in JSON serializable dictionary form.
    '''
    # A list with all referenced sections that might not yet been ensured to be in the
    # resulting partial archive
    referenceds: List[MSection] = []
    # contents keeps track of all sections in the partial archive by keeping their
    # JSON serializable form and placeholder status in a dict
    contents: Dict[MSection, Tuple[dict, bool]] = dict()

    def partial(definition: Definition, section: MSection) -> bool:
        '''
        ``m_to_dict`` partial function that selects what goes into the partial archive
        and what not. It also collects references as a side-effect.
        '''
        if section.m_def == EntryArchive.m_def:
            if definition.m_def == Quantity:
                return True
            return FastAccess.m_def in definition.categories

        if isinstance(definition, Quantity) and isinstance(definition.type, Reference) \
                and FastAccess.m_def in definition.categories:
            # Reference list in partial archives are not supported
            if definition.is_scalar:
                referenced = getattr(section, definition.name)
                if referenced is not None and referenced not in contents:
                    referenceds.append(referenced)

        if isinstance(definition, SubSection):
            return FastAccess.m_def in definition.categories

        return True

    # add the main content
    partial_contents = archive.m_to_dict(partial=partial)

    # add the referenced data
    def add(section, placeholder=False) -> dict:
        '''
        Adds the given section to partial_contents at the right place. If not a placeholder,
        the section's serialization is added (or replacing an existing placeholder).
        Otherwise, an empty dict is added as a placeholder for referenced children.
        '''
        result: Dict[str, Any] = None
        content, content_is_placeholder = contents.get(section, (None, True))
        if content is not None:
            if content_is_placeholder and not placeholder:
                # the placeholder gets replaced later
                pass
            else:
                return content

        if section.m_parent is None:
            contents[section] = partial_contents, False
            return partial_contents

        parent_dict = add(section.m_parent, placeholder=True)
        if placeholder:
            result = {}
        else:
            result = section.m_to_dict(partial=partial)

        sub_section = section.m_parent_sub_section
        if sub_section.repeats:
            sections = parent_dict.setdefault(sub_section.name, [])
            while len(sections) < section.m_parent_index + 1:
                sections.append(None)
            sections[section.m_parent_index] = result
        else:
            parent_dict[sub_section.name] = result

        contents[section] = result, placeholder
        return result

    # we add referenced objects as long as they are added by subsequent serialization
    # of referenced sections to implement the recursive nature of further references in
    # already referenced sections.
    while len(referenceds) > 0:
        referenced = referenceds.pop()
        add(referenced)

    return partial_contents


def write_partial_archive_to_mongo(archive: EntryArchive):
    ''' Partially writes the given archive to mongodb. '''
    mongo_db = infrastructure.mongo_client[config.mongo.db_name]
    mongo_collection = mongo_db['archive']
    mongo_id = archive.section_metadata.calc_id

    partial_archive_dict = create_partial_archive(archive)
    partial_archive_dict['_id'] = mongo_id
    mongo_collection.replace_one(dict(_id=mongo_id), partial_archive_dict, upsert=True)


def read_partial_archive_from_mongo(entry_id: str, as_dict=False) -> Union[EntryArchive, Dict]:
    '''
    Reads the partial archive for the given id from mongodb.

    Arguments:
        entry_id: The entry id for the entry.
        as_dict: Return the JSON serializable dictionary form of the archive not the
            :class:`EntryArchive` form.
    '''
    mongo_db = infrastructure.mongo_client[config.mongo.db_name]
    mongo_collection = mongo_db['archive']
    archive_dict = mongo_collection.find_one(dict(_id=entry_id))

    if as_dict:
        return archive_dict

    return EntryArchive.m_from_dict(archive_dict)


def delete_partial_archives_from_mongo(entry_ids: List[str]):
    mongo_db = infrastructure.mongo_client[config.mongo.db_name]
    mongo_collection = mongo_db['archive']
    mongo_collection.delete_many(dict(_id={'$in': entry_ids}))


def read_partial_archives_from_mongo(entry_ids: List[str], as_dict=False) -> Dict[str, Union[EntryArchive, Dict]]:
    '''
    Reads the partial archives for a set of entries.

    Arguments:
        entry_ids: A list of entry ids.
        as_dict: Return the JSON serializable dictionary form of the archive not the
            :class:`EntryArchive` form.

    Returns:
        A dictionary with entry_ids as keys.
    '''
    mongo_db = infrastructure.mongo_client[config.mongo.db_name]
    mongo_collection = mongo_db['archive']
    archive_dicts = mongo_collection.find(dict(_id={'$in': entry_ids}))

    if as_dict:
        return {archive_dict.pop('_id'): archive_dict for archive_dict in archive_dicts}

    return {
        archive_dict.pop('_id'): EntryArchive.m_from_dict(archive_dict)
        for archive_dict in archive_dicts}


__all_parent_sections: Dict[Section, Tuple[str, Section]] = {}


def _all_parent_sections():
    if len(__all_parent_sections) == 0:
        def add(section):
            for sub_section in section.all_sub_sections.values():
                sub_section_section = sub_section.sub_section.m_resolved()
                __all_parent_sections.setdefault(sub_section_section, []).append((sub_section.qualified_name(), section, ))
                add(sub_section_section)

        add(EntryArchive.m_def)

    return __all_parent_sections


class _Incomplete(Exception): pass


def compute_required_with_referenced(required):
    '''
    Updates the given required dictionary to ensure that references to non required
    sections within a partial fast access archive are included. Only references that
    are directly contained in required are added. References from wildcarded sub sections
    are ignored.

    Returns: A new required dict or None. None is returned if it is unclear if the required
    is only accessing information of fast access partial archives.
    '''
    # TODO this function should be based on the metainfo

    if not isinstance(required, dict):
        return None

    if any(key.startswith('section_run') for key in required):
        return None

    required = dict(**required)

    def add_parent_section(section, child_required):
        parent_sections = _all_parent_sections().get(section, [])
        if len(parent_sections) == 0:
            return [required]

        result = []
        for name, parent_section in parent_sections:
            child_key = name.split('.')[-1]
            for parent_required in add_parent_section(parent_section, None):
                result.append(parent_required.setdefault(child_key, child_required if child_required else {}))

        return result

    def traverse(
            current: Union[dict, str],
            parent: Section = EntryArchive.m_def):

        if isinstance(current, str):
            return

        current_updates = {}
        for key, value in current.items():
            prop = key.split('[')[0]
            prop_definition = parent.all_properties[prop]
            if isinstance(prop_definition, SubSection):
                if FastAccess.m_def not in prop_definition.categories:
                    raise _Incomplete()

                traverse(value, prop_definition.sub_section)
            if isinstance(prop_definition, Quantity) and isinstance(prop_definition.type, Reference):
                current_updates[prop] = '*'
                if FastAccess.m_def not in prop_definition.categories:
                    continue

                target_section_def = prop_definition.type.target_section_def.m_resolved()
                # TODO is this a bug, should the result of the call not be used to traverse?
                add_parent_section(target_section_def, value)
                traverse(value, target_section_def)
        current.update(**current_updates)

    try:
        traverse(dict(**required))
    except _Incomplete:
        # We realized that something is required that is not in the partial archive
        return None

    return required
