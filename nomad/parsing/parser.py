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

from abc import ABCMeta, abstractmethod
import sys
import re
import importlib
import inspect
from unittest.mock import patch
import logging
import os.path
import glob

from nomad import utils
from nomad.parsing.backend import LocalBackend


class Parser(metaclass=ABCMeta):
    """
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.

    Arguments:
        name: The name of the parser
        parser_class_name: Full qualified name of the main parser class. We assume it have one
                           parameter for the backend.
        main_file_re: A regexp that matches main file paths that this parser can handle.
        main_contents_re: A regexp that matches main file headers that this parser can parse.
    """

    @abstractmethod
    def is_mainfile(self, filename: str, mime: str, buffer: str) -> bool:
        """ Checks if a file is a mainfile for the parsers. """
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


class LegacyParser(Parser):
    """
    A parser implementation for legacy NOMAD-coe parsers. It assumes that parsers
    are installed to the python environment. It
    uses regular expessions to match parsers to mainfiles.

    Arguments:
        python_get: the git repository and commit that contains the legacy parser
        parser_class_name: the main parser class that implements NOMAD-coe's
            python-common *ParserInterface*. Instances of this class are currently not reused.
        mainfile_mime_re: A regexp that is used to match against a files mime type
        mainfile_contents_re: A regexp that is used to match the first 1024 bytes of a
            potential mainfile.
        mainfile_name_re: A regexp that is used to match the paths of potential mainfiles
    """
    def __init__(
            self, name: str, parser_class_name: str,
            mainfile_contents_re: str,
            mainfile_mime_re: str = r'text/.*',
            mainfile_name_re: str = r'.*') -> None:

        self.name = name
        self.parser_class_name = parser_class_name
        self._mainfile_mime_re = re.compile(mainfile_mime_re)
        self._mainfile_name_re = re.compile(mainfile_name_re)
        self._mainfile_contents_re = re.compile(mainfile_contents_re)

    def is_mainfile(self, filename: str, mime: str, buffer: str) -> bool:
        return self._mainfile_name_re.match(filename) is not None and \
            self._mainfile_mime_re.match(mime) is not None and \
            self._mainfile_contents_re.search(buffer) is not None

    def run(self, mainfile: str, logger=None) -> LocalBackend:
        # TODO we need a homogeneous interface to parsers, but we dont have it right now.
        # There are some hacks to distringuish between ParserInterface parser and simple_parser
        # using hasattr, kwargs, etc.
        def create_backend(meta_info, logger=None):
            return LocalBackend(meta_info, debug=False, logger=logger)

        module_name = self.parser_class_name.split('.')[:-1]
        parser_class = self.parser_class_name.split('.')[1]
        module = importlib.import_module('.'.join(module_name))
        Parser = getattr(module, parser_class)

        init_signature = inspect.getargspec(Parser.__init__)
        kwargs = dict(
            backend=lambda meta_info: create_backend(meta_info, logger=logger),
            log_level=logging.DEBUG, debug=True)
        kwargs = {key: value for key, value in kwargs.items() if key in init_signature.args}

        with utils.legacy_logger(logger):
            self.parser = Parser(**kwargs)

            with patch.object(sys, 'argv', []):
                backend = self.parser.parse(mainfile)

            if backend is None or not hasattr(backend, 'status'):
                backend = self.parser.parser_context.super_backend

        return backend

    def __repr__(self):
        return self.name


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
            os.path.dirname(filename)
            if len(glob.glob('%s/*.xml*' % filename)) > 0:
                return False

        return is_mainfile
