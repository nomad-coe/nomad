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
import sys
import re
import importlib
import inspect
from unittest.mock import patch
import logging
import os.path
import os
import glob

from nomad import utils, config
from nomad.parsing.backend import LocalBackend


class Parser(metaclass=ABCMeta):
    """
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.
    """

    def __init__(self):
        self.domain = 'DFT'

    @abstractmethod
    def is_mainfile(self, filename: str, mime: str, buffer: bytes, compression: str = None) -> bool:
        """
        Checks if a file is a mainfile for the parsers.

        Arguments:
            filename: The filesystem path to the mainfile
            mime: The mimetype of the mainfile guessed with libmagic
            buffer: The first 2k of the mainfile contents
            compression: The compression of the mainfile ``[None, 'gz', 'bz2']``
        """
        pass

    @abstractmethod
    def run(self, mainfile: str, logger=None) -> LocalBackend:
        """
        Runs the parser on the given mainfile. It uses :class:`LocalBackend` as
        a backend. The meta-info access is handled by the underlying NOMAD-coe parser.

        Args:
            mainfile: A path to a mainfile that this parser can parse.
            logger: A optional logger

        Returns:
            The used :class:`LocalBackend` with status information and result data.
        """


class BrokenParser(Parser):
    """
    A parser implementation that just fails and is used to match mainfiles with known
    patterns of corruption.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'parser/broken'
        self.code_name = 'currupted mainfile'
        self._patterns = [
            re.compile(r'^pid=[0-9]+'),  # some 'mainfile' contain list of log-kinda information with pids
            re.compile(r'^Can\'t open .* library:.*')  # probably bad code runs
        ]

    def is_mainfile(self, filename: str, mime: str, buffer: bytes, compression: str = None) -> bool:

        try:  # Try to open the file as a string for regex matching.
            decoded_buffer = buffer.decode('utf-8')
        except UnicodeDecodeError:
            # This file is binary, and should not be binary
            pass
        else:
            for pattern in self._patterns:
                if pattern.search(decoded_buffer) is not None:
                    return True

        return False

    def run(self, mainfile: str, logger=None) -> LocalBackend:
        raise Exception('Failed on purpose.')


class MatchingParser(Parser):
    """
    A parser implementation that used regular experessions to match mainfiles.

    Arguments:
        mainfile_mime_re: A regexp that is used to match against a files mime type
        mainfile_contents_re: A regexp that is used to match the first 1024 bytes of a
            potential mainfile.
        mainfile_name_re: A regexp that is used to match the paths of potential mainfiles
        domain: The domain that this parser should be used for. Default is 'DFT'.
        supported_compressions: A list of [gz, bz2], if the parser supports compressed files
    """
    def __init__(
            self, name: str, code_name: str,
            mainfile_contents_re: str = None,
            mainfile_binary_header: bytes = None,
            mainfile_mime_re: str = r'text/.*',
            mainfile_name_re: str = r'.*',
            domain='DFT',
            supported_compressions: List[str] = []) -> None:

        super().__init__()
        self.name = name
        self.code_name = code_name
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

    def is_mainfile(self, filename: str, mime: str, buffer: bytes, compression: str = None) -> bool:
        if self._mainfile_binary_header is not None:
            if self._mainfile_binary_header not in buffer:
                return False
        if self._mainfile_contents_re is not None:
            try:  # Try to open the file as a string for regex matching.
                decoded_buffer = buffer.decode('utf-8')
                print('DECODED',decoded_buffer)
            except UnicodeDecodeError as e:
                return False  # We're looking for a string match in a file that can't be converted to string.
            if self._mainfile_contents_re.search(decoded_buffer) is None:
                return False
        return self._mainfile_mime_re.match(mime) is not None and \
            self._mainfile_name_re.fullmatch(filename) is not None and \
            (compression is None or compression in self._supported_compressions)

    def __repr__(self):
        return self.name


class MissingParser(MatchingParser):
    """
    A parser implementation that just fails and is used to match mainfiles with known
    patterns of corruption.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def run(self, mainfile: str, logger=None) -> LocalBackend:
        raise Exception('The code %s is not yet supported.' % self.code_name)


class LegacyParser(MatchingParser):
    """
    A parser implementation for legacy NOMAD-coe parsers. It assumes that parsers
    are installed to the python environment.

    Arguments:
        parser_class_name: the main parser class that implements NOMAD-coe's
        backend_factory: a callable that returns a backend, takes meta_info and logger as argument
    """
    def __init__(self, parser_class_name: str, *args, backend_factory=None, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        self.parser_class_name = parser_class_name
        self.backend_factory = backend_factory

    def run(self, mainfile: str, logger=None) -> LocalBackend:
        # TODO we need a homogeneous interface to parsers, but we dont have it right now.
        # There are some hacks to distringuish between ParserInterface parser and simple_parser
        # using hasattr, kwargs, etc.
        def create_backend(meta_info):
            if self.backend_factory is not None:
                return self.backend_factory(meta_info, logger=logger)

            return LocalBackend(meta_info, debug=False, logger=logger)

        module_name = self.parser_class_name.split('.')[:-1]
        parser_class = self.parser_class_name.split('.')[-1]
        module = importlib.import_module('.'.join(module_name))
        Parser = getattr(module, parser_class)

        init_signature = inspect.getargspec(Parser.__init__)
        kwargs = dict(backend=create_backend, log_level=logging.DEBUG, debug=True)
        kwargs = {key: value for key, value in kwargs.items() if key in init_signature.args}

        with utils.legacy_logger(logger):
            self.parser = Parser(**kwargs)

            with patch.object(sys, 'argv', []):
                backend = self.parser.parse(mainfile)
                os.chdir(config.fs.working_directory)

            if backend is None or not hasattr(backend, 'status'):
                backend = self.parser.parser_context.super_backend

        return backend


class VaspOutcarParser(LegacyParser):
    """
    LegacyParser that only matches mailfiles, if there is no .xml in the
    same directory, i.e. to use the VASP OUTCAR parser in absence of .xml
    output file.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = 'parsers/vaspoutcar'

    def is_mainfile(self, filename: str, *args, **kwargs) -> bool:
        is_mainfile = super().is_mainfile(filename, *args, **kwargs)

        if is_mainfile:
            directory = os.path.dirname(filename)
            if len(glob.glob('%s/*.xml*' % directory)) > 0:
                return False

        return is_mainfile
