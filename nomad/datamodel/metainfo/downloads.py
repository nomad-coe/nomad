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

from nomad.datamodel.data import ArchiveSection
from nomad.metainfo import MSection, Package, Quantity, SubSection

m_package = Package(name='downloads')


class Download(MSection):
    url = Quantity(
        type=str,
        description="""
        Deprecated compatibility field for a downloadable URL.""",
    )

    output = Quantity(
        type=str,
        default='./',
        description="""
        Deprecated compatibility field for the download target path.""",
    )

    extract = Quantity(
        type=bool,
        description="""
        Deprecated compatibility field for extracted downloads.""",
    )


class Downloads(ArchiveSection):
    """
    Deprecated compatibility section.

    Existing archives can still deserialize this section, but processing ignores
    it and will not perform any downloads or trigger follow-up processing.
    """

    description = Quantity(
        type=str,
        description="""Provides some additional description for these downloads.""",
    )

    mainfiles = Quantity(
        type=str,
        shape=['*'],
        description="""
        Deprecated compatibility field for mainfiles that used to be triggered
        after downloads.""",
    )

    skip_download = Quantity(
        type=bool,
        description="""
        Compatibility flag retained for old archives. Processing always ignores
        this section and sets the flag to true.""",
    )

    downloads = SubSection(
        section=Download,
        repeats=True,
        description="""
        Defines URLs and how to download them.""",
    )

    def normalize(self, archive, logger):
        super().normalize(archive, logger)

        if archive.data == self and archive.metadata:
            archive.metadata.entry_type = 'Downloads'
            archive.metadata.entry_name = archive.metadata.mainfile

        logger.warning(
            'Downloads archive sections are deprecated and ignored during processing.'
        )
        self.skip_download = True


m_package.__init_metainfo__()
