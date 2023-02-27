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

from nomad.datamodel.data import EntryData
from . import Normalizer


class MetainfoNormalizer(Normalizer):
    domain = None

    def normalize_section(self, section, logger):
        normalize = None
        try:
            normalize = getattr(section, 'normalize')
        except Exception:
            pass

        if normalize:
            try:
                normalize(self.entry_archive, logger)
            except Exception as e:
                logger.error(
                    'could not normalize section',
                    section=section.m_def.name,
                    exc_info=e)

    def normalize(self, logger=None) -> None:
        if logger is None:
            from nomad import utils
            logger = utils.get_logger(__name__)

        logger = logger.bind(normalizer=self.__class__.__name__)
        self.logger = logger

        for sub_section in self.entry_archive.m_contents():
            if isinstance(sub_section, EntryData):
                for section in sub_section.m_all_contents():
                    self.normalize_section(section, logger)
            self.normalize_section(sub_section, logger)

        self.normalize_section(self.entry_archive, logger)
