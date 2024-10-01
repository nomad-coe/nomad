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
import math
import re
from typing import Callable, Any

import numpy as np

from nomad import utils
from nomad.metainfo import (
    Section,
    AnnotationModel,
    MSection,
    SubSection,
    Property,
    MetainfoError,
)
from nomad.units import ureg

# ../entries/<entry_id>/archive#<path>
# /entries/<entry_id>/archive#<path>
_regex_form_a = re.compile(r'^(?:\.\.)?/entries/([^?]+)/(archive|raw)#([^?]+?)$')

# ../upload/<upload_id>/archive/<entry_id>#<path>
# /uploads/<upload_id>/archive/<entry_id>#<path>
# <installation>/uploads/<upload_id>/archive/<entry_id>#<path>
_regex_form_b = re.compile(
    r'^([^?]+?)?/uploads?/([\w=-]*)/?(archive|raw)/([^?]+?)#([^?]+?)$'
)


def parse_path(url: str, upload_id: str = None):
    """
    Parse a reference path.

    The upload_id of current upload is taken as the input to account for that the relative reference has no
    information about the upload_id, and it may also contain no entry_id. Has to know the upload_id when only
    path to mainfile is given.

    Returns:
        (installation, upload_id, entry_id_or_mainfile, kind, path): successfully parsed path
        None: fail to parse path
    """

    url_match = _regex_form_b.match(url)
    if not url_match:
        # try another form
        url_match = _regex_form_a.match(url)
        if not url_match:
            # not valid
            return None

        entry_id_or_mainfile = url_match.group(1)
        kind = url_match.group(2)  # archive or raw
        path = url_match.group(3)
        if not path.startswith('/'):
            path = '/' + path

        return None, upload_id, entry_id_or_mainfile, kind, path

    installation = url_match.group(1)
    if installation == '':
        installation = None
    elif installation == '..':
        installation = None

    # if empty, it is a local reference to the same upload, use the current upload_id
    other_upload_id = upload_id if url_match.group(2) == '' else url_match.group(2)

    kind = url_match.group(3)  # archive or raw
    entry_id_or_mainfile = url_match.group(4)
    path = url_match.group(5)
    if not path.startswith('/'):
        path = '/' + path

    if kind == 'archive':
        if entry_id_or_mainfile.startswith('mainfile/'):
            entry_id_or_mainfile = utils.generate_entry_id(
                other_upload_id, entry_id_or_mainfile.replace('mainfile/', '')
            )
        elif '/' in entry_id_or_mainfile:  # should not contain '/' in entry_id
            return None

    return installation, other_upload_id, entry_id_or_mainfile, kind, path


def create_custom_mapping(
    section_def: Section,
    annotation_type: AnnotationModel,
    annotation_name: str,
    annotation_attr: str,
) -> list[tuple[str, Callable[[MSection, Any], MSection]]]:
    mapping: list[tuple[str, Callable[[MSection, Any], MSection]]] = []
    base_sections: set[MSection] = set()

    def add_section_def(
        section_def: Section, path: list[tuple[SubSection, Section]], base_sections
    ):
        base_sections.add(section_def)
        properties: set[Property] = set()
        for quantity in section_def.all_quantities.values():
            if quantity in properties:
                continue
            properties.add(quantity)

            annotation = quantity.m_get_annotations(annotation_name)
            annotation = annotation[0] if isinstance(annotation, list) else annotation
            annotation = annotation_type.parse_obj(annotation) if annotation else None
            if annotation and getattr(annotation, annotation_attr, None):
                prop = getattr(annotation, annotation_attr)

                def set_value(
                    section: MSection,
                    value,
                    path=path,
                    quantity=quantity,
                    annotation=annotation,
                ):
                    for sub_section, path_section_def in path:
                        next_section = None
                        try:
                            next_section = section.m_get_sub_section(sub_section, -1)
                        except (KeyError, IndexError):
                            pass
                        if not next_section:
                            next_section = path_section_def.section_cls()
                            section.m_add_sub_section(sub_section, next_section)
                        section = next_section

                    if annotation and annotation.unit:
                        value *= ureg(annotation.unit)

                    # NaN values are not supported in the metainfo. Set as None
                    # which means that they are not stored.
                    if isinstance(value, float) and math.isnan(value):
                        value = None

                    if isinstance(value, (int, float, str)):
                        value = np.array([value])

                    if value is not None:
                        if len(value.shape) == 1 and len(quantity.shape) == 0:
                            if len(value) == 1:
                                value = value[0]
                            elif len(value) == 0:
                                value = None
                            else:
                                raise MetainfoError(
                                    f'The shape of {quantity.name} does not match the given data.'
                                )
                        elif len(value.shape) != len(quantity.shape):
                            raise MetainfoError(
                                f'The shape of {quantity.name} does not match the given data.'
                            )

                    section.m_set(quantity, value)

                mapping.append((prop, set_value))

        for sub_section in section_def.all_sub_sections.values():
            if sub_section in properties:
                continue
            next_base_section = sub_section.sub_section
            properties.add(sub_section)
            for sub_section_section in next_base_section.all_inheriting_sections + [
                next_base_section if next_base_section not in base_sections else None
            ]:
                if sub_section_section:
                    add_section_def(
                        sub_section_section,
                        path
                        + [
                            (
                                sub_section,
                                sub_section_section,
                            )
                        ],
                        base_sections,
                    )

    add_section_def(section_def, [], base_sections)
    return mapping
