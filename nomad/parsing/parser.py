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

from typing import List, Set, Dict, Union
from abc import ABCMeta, abstractmethod
import re
import os
import os.path
from functools import lru_cache
import importlib

from nomad import config, utils
from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.metainfo import Package


class Parser(metaclass=ABCMeta):
    '''
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.
    '''
    name = "parsers/parser"
    level = 0
    creates_children = False
    '''
    Level 0 parsers are run first, then level 1, and so on. Normally the value should be 0,
    use higher values only when a parser depends on other parsers.
    '''

    def __init__(self):
        self.domain = 'dft'

    @abstractmethod
    def is_mainfile(
            self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
            compression: str = None) -> Union[bool, Set[str]]:
        '''
        Checks if a file is a mainfile for the parser. Should return True or a set of
        *keys* (non-empty strings) if it is a mainfile, otherwise a falsey value.

        The option to return a set of keys should only be used by parsers that have
        `creates_children == True`. These create multiple entries for one mainfile, namely
        a *main* entry and some number of *child* entries. Most parsers, however, have
        `creates_children == False` and thus only generate a main entry, no child entries,
        and these should thus just return a boolean value.

        If the return value is a set of keys, a main entry will be created when parsing,
        plus one child entry for each key in the returned set. The key value will be stored
        in the field `mainfile_key` of the corresponding child entry. Main entries have
        `mainfile_key == None`.

        The combination (`upload_id`, `mainfile`, `mainfile_key`) uniquely identifies an entry
        (regardless of it's a main entry or child entry).

        Arguments:
            filename: The filesystem path to the mainfile
            mime: The mimetype of the mainfile guessed with libmagic
            buffer: The first 2k of the mainfile contents
            compression: The compression of the mainfile ``[None, 'gz', 'bz2']``
        '''
        pass

    @abstractmethod
    def parse(
            self,
            mainfile: str,
            archive: EntryArchive,
            logger=None,
            child_archives: Dict[str, EntryArchive] = None) -> None:
        '''
        Runs the parser on the given mainfile and populates the result in the given
        archive root_section. It allows to be run repeatedly for different mainfiles.

        Args:
            mainfile: A path to a mainfile that this parser can parse.
            archive: An instance of the section :class:`EntryArchive`. It might contain
                a section ``metadata`` with information about the entry.
            logger: A optional logger
            child_archives: a dictionary with {mainfile_key : EntryArchive} for each child,
                for the parse function to populate with data.
        '''
        pass

    def after_normalization(self, archive: EntryArchive, logger=None) -> None:
        '''
        This is called after the archive produced by `parsed` has been normalized. This
        allows to apply additional code-specific processing steps based on the normalized data.

        Args:
            archive: An instance of the section :class:`EntryArchive`. It might contain
                a section ``metadata`` with information about the entry.
            logger: A optional logger
        '''
        pass

    @classmethod
    def main(cls, mainfile, mainfile_keys: List[str] = None):
        archive = EntryArchive()
        archive.m_create(EntryMetadata)
        if mainfile_keys:
            child_archives = {}
            for mainfile_key in mainfile_keys:
                child_archive = EntryArchive()
                child_archive.m_create(EntryMetadata)
                child_archives[mainfile_key] = child_archive
            kwargs = dict(child_archives=child_archives)
        else:
            kwargs = {}

        cls().parse(mainfile, archive, **kwargs)  # pylint: disable=no-value-for-parameter
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

    def parse(self, mainfile: str, archive, logger=None, child_archives=None):
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

    def parse(self, mainfile: str, archive: EntryArchive, logger=None, child_archives=None) -> None:
        raise NotImplementedError()

    def __repr__(self):
        return self.name


class MatchingParserInterface(MatchingParser):
    '''
    An interface to the NOMAD parsers.

    Arguments:
        parser_class_name: concatenation of module path and parser class name
    '''
    def __init__(self, parser_class_name: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parser_class_name = parser_class_name
        self._mainfile_parser = None

    @property
    def mainfile_parser(self):
        if self._mainfile_parser is None:
            try:
                module_path, parser_class = self._parser_class_name.rsplit('.', 1)
                module = importlib.import_module(module_path)
                self._mainfile_parser = getattr(module, parser_class)()
            except Exception as e:
                logger = utils.get_logger(__name__)
                logger.error('Error importing parser.', exc_info=e)
                raise e
        return self._mainfile_parser

    def parse(self, mainfile: str, archive: EntryArchive, logger=None, child_archives=None):
        self.mainfile_parser.parse(mainfile, archive, logger)


class ArchiveParser(MatchingParser):
    def __init__(self):
        super().__init__(
            name='parsers/archive',
            code_name=config.services.unavailable_value,
            domain=None,
            mainfile_mime_re='.*',
            mainfile_name_re=r'.*(archive|metainfo)\.(json|yaml|yml)$')

    def parse(self, mainfile: str, archive: EntryArchive, logger=None, child_archives=None):
        if mainfile.endswith('.json'):
            import json
            with open(mainfile, 'rt') as f:
                archive_data = json.load(f)
        else:
            import yaml
            with open(mainfile, 'rt') as f:
                archive_data = yaml.load(f, Loader=getattr(yaml, 'FullLoader'))

        metadata_data = archive_data.get(EntryArchive.metadata.name, None)

        if metadata_data is not None:
            # Setting metadata in this way is not supported (any more)
            del(archive_data[EntryArchive.metadata.name])

        # ensure that definitions are parsed first to make them available for the
        # parsing itself.
        if 'definitions' in archive_data:
            archive.definitions = Package.m_from_dict(
                archive_data['definitions'], m_context=archive.m_context)
            del archive_data['definitions']

        archive.m_update_from_dict(archive_data)


class MissingParser(MatchingParser):
    '''
    A parser implementation that just fails and is used to match mainfiles with known
    patterns of corruption.
    '''
    name = "parsers/missing"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def parse(self, mainfile: str, archive: EntryArchive, logger=None, child_archives=None):
        raise Exception('The code %s is not yet supported.' % self.code_name)
