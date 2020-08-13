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

from typing import List
from abc import ABCMeta, abstractmethod
import re
import importlib

from nomad.metainfo import Environment
from nomad.datamodel import EntryArchive


class Parser(metaclass=ABCMeta):
    '''
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.

    TODO: There are currently two "run" functions. :func:`run` and :func:`parse`.
    Because we are in the middle of transitioning out of the backend dependence we currently
    have both, where 'run' creates a backend and 'parse' simply gets an archive that the
    parser is supposed to populate. Eventually, we will only have the 'parse' function.
    '''
    name = "parsers/parser"

    def __init__(self):
        self.domain = 'dft'
        self._metainfo_env: Environment = None

    @property
    def metainfo_env(self):
        return self._metainfo_env

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
    def run(self, mainfile: str, logger=None):
        '''
        Runs the parser on the given mainfile. It uses :class:`Backend` as
        a backend. The meta-info access is handled by the underlying NOMAD-coe parser.

        Args:
            mainfile: A path to a mainfile that this parser can parse.
            logger: A optional logger

        Returns:
            The used :class:`Backend` with status information and result data.
        '''

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

    def run(self, mainfile: str, logger=None):
        raise Exception('Failed on purpose.')


class MatchingParser(Parser):
    '''
    A parser implementation that used regular experessions to match mainfiles.

    Arguments:
        code_name: The name of the code or input format
        code_homepage: The homepage of the code or input format
        mainfile_mime_re: A regexp that is used to match against a files mime type
        mainfile_contents_re: A regexp that is used to match the first 1024 bytes of a
            potential mainfile.
        mainfile_name_re: A regexp that is used to match the paths of potential mainfiles
        domain: The domain that this parser should be used for. Default is 'dft'.
        supported_compressions: A list of [gz, bz2], if the parser supports compressed files
    '''
    def __init__(
            self, name: str, code_name: str, code_homepage: str = None,
            mainfile_contents_re: str = None,
            mainfile_binary_header: bytes = None,
            mainfile_mime_re: str = r'text/.*',
            mainfile_name_re: str = r'.*',
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
        # Assign private variable this way to avoid static check issue.
        if mainfile_contents_re is not None:
            self._mainfile_contents_re = re.compile(mainfile_contents_re)
        else:
            self._mainfile_contents_re = None
        self._supported_compressions = supported_compressions

    def is_mainfile(
            self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
            compression: str = None) -> bool:

        if self._mainfile_binary_header is not None:
            if self._mainfile_binary_header not in buffer:
                return False
        if self._mainfile_contents_re is not None:
            if decoded_buffer is not None:
                if self._mainfile_contents_re.search(decoded_buffer) is None:
                    return False
            else:
                return False
        return self._mainfile_mime_re.match(mime) is not None and \
            self._mainfile_name_re.fullmatch(filename) is not None and \
            (compression is None or compression in self._supported_compressions)

    def __repr__(self):
        return self.name


class FairdiParser(MatchingParser):

    def run(self, mainfile: str, logger=None):
        from .legacy import Backend
        python_module = importlib.import_module(self.__module__ + '.metainfo')
        metainfo = getattr(python_module, 'm_env')
        backend = Backend(metainfo, domain=self.domain, logger=logger)
        self.parse(mainfile, backend.entry_archive, logger=logger)
        return backend

    def parse(self, mainfile: str, archive: EntryArchive, logger=None):
        raise NotImplementedError()

    @classmethod
    def main(cls, mainfile):
        archive = EntryArchive()
        cls().parse(mainfile, archive)  # pylint: disable=no-value-for-parameter
        return archive


class MissingParser(MatchingParser):
    '''
    A parser implementation that just fails and is used to match mainfiles with known
    patterns of corruption.
    '''
    name = "parsers/missing"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def run(self, mainfile: str, logger=None):
        raise Exception('The code %s is not yet supported.' % self.code_name)
