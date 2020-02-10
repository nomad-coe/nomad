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
Adds elastic search support to the metainfo.
"""

from . import Section, MSection


def elastic_mapping(section: Section, base_cls: type) -> type:
    """ Creates an elasticsearch_dsl document class from a section definition. """

    dct = {
        name: quantity.m_annotations['elastic']['type']()
        for name, quantity in section.all_quantities.items()
        if 'elastic' in quantity.m_annotations}

    return type(section.name, (base_cls,), dct)


def elastic_obj(source: MSection, target_cls: type):
    if source is None:
        return None

    assert isinstance(source, MSection), '%s must be an MSection decendant' % source.__class__.__name__

    target = target_cls()

    for name, quantity in source.m_def.all_quantities.items():
        elastic_annotation = quantity.m_annotations.get('elastic')
        if elastic_annotation is None:
            continue

        if 'mapping' in elastic_annotation:
            value = elastic_annotation['mapping'](source)
        else:
            value = getattr(source, name)

        setattr(target, name, value)

    return target
