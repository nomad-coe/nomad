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

from nomad.datamodel import EntryArchive
from nomad.datamodel.data import ArchiveSection
from nomad.datamodel import EntryArchive
from typing import Optional
from . import Normalizer


class MetainfoNormalizer(Normalizer):
    domain: Optional[str] = None

    def normalize_section(self, archive: EntryArchive, section, logger):
        normalize = None
        try:
            normalize = getattr(section, 'normalize')
        except Exception:
            pass

        if normalize:
            try:
                normalize(archive, logger)
            except Exception as e:
                logger.error(
                    'could not normalize section',
                    section=section.m_def.name,
                    exc_info=e,
                )

    def normalize(self, archive: EntryArchive, logger=None) -> None:
        if logger is None:
            from nomad import utils

            logger = utils.get_logger(__name__)

        logger = logger.bind(normalizer=self.__class__.__name__)
        self.logger = logger

        def _normalize(section):
            sub_sections = [sub_section for sub_section in section.m_contents()]
            # TODO eln test failing because non-ArchiveSection may be normalization
            sub_sections.sort(
                key=lambda x: x.normalizer_level if isinstance(x, ArchiveSection) else 0
            )
            for sub_section in sub_sections:
                _normalize(sub_section)
            self.normalize_section(archive, section, logger)

        _normalize(archive)
