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

import os.path

from nomad.metainfo import MSection


class ArchiveSection(MSection):
    '''
    Base class for sections in a NOMAD archive. Provides a framework for custom
    section normalization via the `normalize` function.
    '''
    def normalize(self, archive, logger):
        '''
        Is called during entry normalization. If you overwrite this with custom
        normalization code, make sure to call `super(YourClass, self).normalize(archive, logger)`.
        Otherwise, not all normalize functions might be called for section definitions
        with multiple base-classes.

        Arguments:
            archive: The whole archive that is normalized.
            logger: The structlog logger used during normalization.
        '''
        pass


class EntryData(ArchiveSection):
    '''
    An empty base section definition. This can be used to add new top-level sections
    to an entry.
    '''

    def normalize(self, archive, logger):
        super(EntryData, self).normalize(archive, logger)

        from nomad.datamodel.results import Results

        archive.metadata.entry_type = self.m_def.name
        if archive.metadata.entry_name is None and archive.metadata.mainfile:
            archive.metadata.entry_name = os.path.basename(archive.metadata.mainfile)

        if not archive.results:
            archive.results = Results()
