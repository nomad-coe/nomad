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
This module contains classes that allow to represent the core
nomad data entities (entries/calculations, users, datasets) on a high level of abstraction
independent from their representation in the different modules
:py:mod:`nomad.processing`, :py:mod:`nomad.coe_repo`, :py:mod:`nomad.parsing`,
:py:mod:`nomad.search`, :py:mod:`nomad.app`, :py:mod:`nomad.migration`.

It is not about representing every detail, but those parts that are directly involved in
api, processing, mirroring, or other 'infrastructure' operations.

The class :class:`User` is used to represent users and their attributes.

.. autoclass:: nomad.datamodel.User
    :members:

The class :class:`Dataset` is used to represent datasets and their attributes.

.. autoclass:: nomad.datamodel.Dataset
    :members:

The class :class:`UserMetadata` is used to represent user determined entry metadata.

.. autoclass:: nomad.datamodel.UserMetadata
    :members:

The class :class:`EntryMetadata` is used to represent all metadata about an entry.

.. autoclass:: nomad.datamodel.EntryMetadata
    :members:
'''

from .dft import DFTMetadata
from .ems import EMSMetadata
from .metainfo import Dataset, User, EditableUserMetadata, UserMetadata, EntryMetadata

domains = {
    'dft': {
        'metadata': DFTMetadata,
        'metainfo_all_package': 'all.nomadmetainfo.json',
        'root_section': 'section_run'
    },
    'ems': {
        'metadata': EMSMetadata,
        'metainfo_all_package': 'all.experimental.nomadmetainfo.json',
        'root_section': 'section_experiment'
    }
}

root_sections = [domain['root_section'] for domain in domains.values()] + ['section_entry_info']
