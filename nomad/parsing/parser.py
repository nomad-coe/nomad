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

from typing import List
from abc import ABCMeta, abstractmethod
import re
import os
import os.path
from functools import lru_cache

from nomad import config
from nomad.datamodel import EntryArchive, UserProvidableMetadata, EntryMetadata


class Parser(metaclass=ABCMeta):
    '''
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.
    '''
    name = "parsers/parser"

    def __init__(self):
        self.domain = 'dft'

    @abstractmethod
    def is_mainfile(
            self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
            compression: str = None) -> bool:
        '''
        Checks if a file is a mainfile for the parsers.

        Arguments:
            filename: The filesystem path to the mainfile
            mime: The mimetype of the mainfile guessed with libmagic
            buffer: The first 2k of the mainfile contents
            compression: The compression of the mainfile ``[None, 'gz', 'bz2']``
        '''
        pass

    @abstractmethod
    def parse(self, mainfile: str, archive: EntryArchive, logger=None) -> None:
        '''
        Runs the parser on the given mainfile and populates the result in the given
        archive root_section. It allows to be run repeatedly for different mainfiles.

        Args:
            mainfile: A path to a mainfile that this parser can parse.
            archive: An instance of the section :class:`EntryArchive`. It might contain
                a ``section_metadata`` with information about the entry.
            logger: A optional logger
        '''
        pass

    @classmethod
    def main(cls, mainfile):
        archive = EntryArchive()
        archive.m_create(EntryMetadata)
        cls().parse(mainfile, archive)  # pylint: disable=no-value-for-parameter
        return archive


class BrokenParser(Parser):
    '''
    A parser implementation that just fails and is used to match mainfiles with known
    patterns of corruption.
    '''
    name = 'parsers/broken'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.code_name = 'currupted mainfile'
        self._patterns = [
            re.compile(r'^pid=[0-9]+'),  # some 'mainfile' contain list of log-kinda information with pids
            re.compile(r'^Can\'t open .* library:.*')  # probably bad code runs
        ]

    def is_mainfile(
            self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
            compression: str = None) -> bool:

        if decoded_buffer is not None:
            for pattern in self._patterns:
                if pattern.search(decoded_buffer) is not None:
                    return True

        return False

    def parse(self, mainfile: str, archive, logger=None):
        raise Exception('Failed on purpose.')


class MatchingParser(Parser):
    '''
    A parser implementation that uses regular expressions to match mainfiles.

    Arguments:
        code_name: The name of the code or input format
        code_homepage: The homepage of the code or input format
        mainfile_mime_re: A regexp that is used to match against a files mime type
        mainfile_contents_re: A regexp that is used to match the first 1024 bytes of a
            potential mainfile.
        mainfile_name_re: A regexp that is used to match the paths of potential mainfiles
        mainfile_alternative: If True files are mainfile if no mainfile_name_re matching file
            is present in the same directory.
        domain: The domain that this parser should be used for. Default is 'dft'.
        supported_compressions: A list of [gz, bz2], if the parser supports compressed files
    '''
    def __init__(
            self, name: str, code_name: str, code_homepage: str = None,
            mainfile_contents_re: str = None,
            mainfile_binary_header: bytes = None,
            mainfile_binary_header_re: bytes = None,
            mainfile_mime_re: str = r'text/.*',
            mainfile_name_re: str = r'.*',
            mainfile_alternative: bool = False,
            domain='dft',
            supported_compressions: List[str] = []) -> None:

        super().__init__()
        self.name = name
        self.code_name = code_name
        self.code_homepage = code_homepage
        self.domain = domain
        self._mainfile_binary_header = mainfile_binary_header
        self._mainfile_mime_re = re.compile(mainfile_mime_re)
        self._mainfile_name_re = re.compile(mainfile_name_re)
        self._mainfile_alternative = mainfile_alternative
        # Assign private variable this way to avoid static check issue.
        if mainfile_contents_re is not None:
            self._mainfile_contents_re = re.compile(mainfile_contents_re)
        else:
            self._mainfile_contents_re = None
        if mainfile_binary_header_re is not None:
            self._mainfile_binary_header_re = re.compile(mainfile_binary_header_re)
        else:
            self._mainfile_binary_header_re = None
        self._supported_compressions = supported_compressions

        self._ls = lru_cache(maxsize=16)(lambda directory: os.listdir(directory))

    def is_mainfile(
            self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
            compression: str = None) -> bool:

        if self._mainfile_binary_header is not None:
            if self._mainfile_binary_header not in buffer:
                return False

        if self._mainfile_binary_header_re is not None:
            if buffer is not None:
                if self._mainfile_binary_header_re.search(buffer) is None:
                    return False
            else:
                return False

        if self._mainfile_contents_re is not None:
            if decoded_buffer is not None:
                if self._mainfile_contents_re.search(decoded_buffer) is None:
                    return False
            else:
                return False

        if self._mainfile_mime_re.match(mime) is None:
            return False

        if compression is not None and compression not in self._supported_compressions:
            return False

        if self._mainfile_name_re.fullmatch(filename) is None:
            if not self._mainfile_alternative:
                return False

            directory = os.path.dirname(filename)
            for sibling in self._ls(directory):
                sibling = os.path.join(directory, sibling)
                sibling_is_mainfile = sibling != filename and \
                    self._mainfile_name_re.fullmatch(sibling) is not None and \
                    os.path.isfile(sibling)
                if sibling_is_mainfile:
                    return False

        return True

    def parse(self, mainfile: str, archive: EntryArchive, logger=None) -> None:
        raise NotImplementedError()

    def __repr__(self):
        return self.name


# For backward compatibility
FairdiParser = MatchingParser


class ArchiveParser(MatchingParser):
    def __init__(self):
        super().__init__(
            name='parsers/archive',
            code_name=config.services.unavailable_value,
            domain=None,
            mainfile_mime_re=r'application/json',
            mainfile_name_re=r'.*(archive|metainfo)\.json',
            mainfile_contents_re=r'section_metadata|results')

    def parse(self, mainfile: str, archive: EntryArchive, logger=None):
        import json
        with open(mainfile, 'rt') as f:
            archive_data = json.load(f)

        metadata_data = archive_data.get(EntryArchive.section_metadata.name, None)

        if metadata_data is not None:
            metadata = archive.section_metadata
            for key, value in metadata_data.items():
                if UserProvidableMetadata.m_def not in getattr(EntryMetadata, key).categories:
                    continue

                metadata.m_update_from_dict({key: value})
            del(archive_data[EntryArchive.section_metadata.name])

        archive.m_update_from_dict(archive_data)


class MissingParser(MatchingParser):
    '''
    A parser implementation that just fails and is used to match mainfiles with known
    patterns of corruption.
    '''
    name = "parsers/missing"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def parse(self, mainfile: str, archive: EntryArchive, logger=None):
        raise Exception('The code %s is not yet supported.' % self.code_name)
