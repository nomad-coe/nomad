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

"""
The *parsing* modules is currenlty an abstraction for the existin NOMAD-coe parsers.
The parser code is used via :mod:`nomad.dependencies`. This module redefines
some of the old NOMAD-coe python-common functionality to create a more coherent
interface to the parsers.

Assumption about parsers
------------------------
For now, we make a few assumption about parsers
- they always work on the same *meta-info* version
- they have no conflicting python requirments
- they can be loaded at the same time and can be used within the same python process
- they are uniquely identified by a GIT URL and publicly accessible
- their version is uniquly identified by a GIT commit SHA

Each parser is defined via an instance of :class:`Parser`.

.. autoclass:: nomad.parsing.Parser
    :members:

The parser definitions are available via the following two variables.

.. autodata:: nomad.parsing.parsers
.. autodata:: nomad.parsing.parser_dict

Parsers in NOMAD-coe use a *backend* to create output.

.. autoclass:: nomad.parsing.AbstractParserBackend
    :members:
.. autoclass:: nomad.parsing.LocalBackend
    :members:
"""

from typing import TextIO, Tuple, List, Any, Callable, IO
from abc import ABCMeta, abstractmethod
from io import StringIO
import sys
import json
import re
import importlib
import inspect
from unittest.mock import patch
import io

from nomadcore.local_backend import LocalBackend as LegacyLocalBackend
from nomadcore.local_backend import Section, Results

from nomad.dependencies import dependencies_dict as dependencies, PythonGit
from nomad.utils import get_logger

logger = get_logger(__name__)

ParserStatus = Tuple[str, List[str]]


class DelegatingMeta(ABCMeta):
    def __new__(meta, name, bases, dct):
        abstract_method_names = frozenset.union(*(base.__abstractmethods__ for base in bases))
        for name in abstract_method_names:
            if name not in dct:
                dct[name] = DelegatingMeta._make_delegator_method(name)

        return super(DelegatingMeta, meta).__new__(meta, name, bases, dct)

    @staticmethod
    def _make_delegator_method(name):
        def delegator(self, *args, **kwargs):
            return getattr(self._delegate, name)(*args, **kwargs)
        return delegator


class BadContextURI(Exception):
    pass


class WrongContextState(Exception):
    pass


class AbstractParserBackend(metaclass=ABCMeta):
    """
    This ABS provides the parser backend interface used by the NOMAD-coe parsers
    and normalizers.
    """
    @abstractmethod
    def metaInfoEnv(self):
        """ Returns the meta info used by this backend. """
        pass

    @abstractmethod
    def startedParsingSession(
            self, mainFileUri, parserInfo, parserStatus=None, parserErrors=None):
        """
        Should be called when the parsing starts.
        ParserInfo should be a valid json dictionary.
        """
        pass

    @abstractmethod
    def finishedParsingSession(
            self, parserStatus, parserErrors, mainFileUri=None, parserInfo=None,
            parsingStats=None):
        """ Called when the parsing finishes. """
        pass

    @abstractmethod
    def openContext(self, contextUri: str):
        """ Open existing archive data to introduce new data into an existing section. """
        pass

    @abstractmethod
    def closeContext(self, contextUri: str):
        """ Close priorly opened existing archive data again. """
        pass

    @abstractmethod
    def openSection(self, metaName):
        """ Opens a new section and returns its new unique gIndex. """
        pass

    @abstractmethod
    def closeSection(self, metaName, gIndex):
        """
        Closes the section with the given meta name and index. After this, no more
        value can be added to this section.
        """
        pass

    @abstractmethod
    def openNonOverlappingSection(self, metaName):
        """ Opens a new non overlapping section. """
        pass

    @abstractmethod
    def closeNonOverlappingSection(self, metaName):
        """
        Closes the current non overlapping section for the given meta name. After
        this, no more value can be added to this section.
        """
        pass

    @abstractmethod
    def openSections(self):
        """ Returns the sections that are still open as metaName, gIndex tuples. """
        pass

    @abstractmethod
    def addValue(self, metaName, value, gIndex=-1):
        """
        Adds a json value for the given metaName. The gIndex is used to identify
        the right parent section.
        """
        pass

    @abstractmethod
    def addRealValue(self, metaName, value, gIndex=-1):
        """
        Adds a float value for the given metaName. The gIndex is used to identify
        the right parent section.
        """
        pass

    @abstractmethod
    def addArray(self, metaName, shape, gIndex=-1):
        """
        Adds an unannitialized array of the given shape for the given metaName.
        The gIndex is used to identify the right parent section.
        This is neccessary before array values can be set with :func:`setArrayValues`.
        """

    @abstractmethod
    def setArrayValues(self, metaName, values, offset=None, gIndex=-1):
        """
        Adds values of the given numpy array to the last array added for the given
        metaName and parent gIndex.
        """
        pass

    @abstractmethod
    def addArrayValues(self, metaName, values, gIndex=-1):
        """
        Adds an array with the given numpy array values for the given metaName and
        parent section gIndex.
        """
        pass

    @abstractmethod
    def pwarn(self, msg):
        """ Used to catch parser warnings. """
        pass

    # The following are extensions to the origin NOMAD-coe parser backend. And allow
    # access to existing data

    @property
    @abstractmethod
    def data(self) -> Results:
        pass

    @abstractmethod
    def get_sections(self, meta_name: str) -> List[int]:
        """ Return all gIndices for existing sections of the given meta_name. """
        pass

    @abstractmethod
    def get_value(self, metaName: str, g_index=-1) -> Any:
        """
        Return the value set to the given meta_name in its parent section of the given index.
        An index of -1 (default) is only allowed if there is exactly one parent section.
        """
        pass


class JSONStreamWriter():
    START = 0
    OBJECT = 1
    ARRAY = 2
    KEY_VALUE = 3

    """
    A generator that allows to output JSON based on calling 'event' functions.
    Its pure python and could be replaced by some faster implementation, e.g. yajl-py.
    It uses standard json decode to write values. This allows to mix streaming with
    normal encoding. Furthermore, this writer understands numpy values and encodes
    them properly.

    Arguments:
        file: A file like to write to.
        pretty: True to indent and use separators.

    Raises:
        AssertionError: If methods were called in a non JSON fashion. Call :func:`close`
        to make sure everything was closed properly.
    """
    def __init__(self, file, pretty=False):
        self._fp = file
        self._pretty = pretty

        self._indent = ''  # the current indent
        self._separators = ['']  # a stack of the next necessary separator
        self._states = [JSONStreamWriter.START]  # a stack of what is currenty open

    def _write(self, str):
        self._fp.write(str)

    def _write_seperator(self):
        self._write(self._separators.pop())

    def _seperator_with_newline(self, base=None):
        pretty_ext = ('\n%s' % self._indent) if self._pretty else ''
        if base is None:
            return pretty_ext
        else:
            return '%s%s' % (base, pretty_ext)

    def _open(self, open_char):
        self._write_seperator()
        self._write(open_char)
        self._indent = '%s  ' % self._indent
        self._separators.append(self._seperator_with_newline())

    def _close(self, close_char):
        self._separators.pop()
        self._indent = self._indent[:-2]
        self._write(self._seperator_with_newline())
        self._write(close_char)
        self._separators.append(self._seperator_with_newline(','))

    def open_object(self):
        assert self._states[-1] != JSONStreamWriter.OBJECT, "Cannot open object in object."
        if self._states[-1] == JSONStreamWriter.KEY_VALUE:
            self._states.pop()
        self._open('{')
        self._states.append(JSONStreamWriter.OBJECT)

    def close_object(self):
        assert self._states.pop() == JSONStreamWriter.OBJECT, "Can only close object in object."
        self._close('}')

    def open_array(self):
        assert self._states[-1] != JSONStreamWriter.OBJECT, "Cannot open array in object."
        if self._states[-1] == JSONStreamWriter.KEY_VALUE:
            self._states.pop()
        self._open('[')
        self._states.append(JSONStreamWriter.ARRAY)

    def close_array(self):
        assert self._states.pop() == JSONStreamWriter.ARRAY, "Can only close array in array."
        self._close(']')

    def key_value(self, key, value):
        self.key(key)
        self.value(value)

    def key(self, key: str):
        assert self._states[-1] == JSONStreamWriter.OBJECT, "Key can only be in objects."
        self._write_seperator()
        json.dump(key, self._fp)
        self._separators.append(': ' if self._pretty else ':')
        self._states.append(JSONStreamWriter.KEY_VALUE)

    @staticmethod
    def _json_serializable_value(value):
        if hasattr(value, 'tolist'):
            # run tolist of pentential numpy array types
            return value.tolist()

        else:
            return value

    def write(self, str):
        if self._pretty:
            str = str.replace('\n', '\n%s' % self._indent)
        self._fp.write(str)

    def value(self, value):
        assert self._states[-1] != JSONStreamWriter.OBJECT, "Values can not be in objects."
        if self._states[-1] == JSONStreamWriter.KEY_VALUE:
            self._states.pop()

        self._write_seperator()
        json.dump(
            JSONStreamWriter._json_serializable_value(value), self,
            indent=2 if self._pretty else 1,
            separators=(', ', ': ') if self._pretty else (',', ':'))
        self._separators.append(self._seperator_with_newline(','))

    def close(self):
        assert self._states[-1] == JSONStreamWriter.START, "Something was not closed."


class LegacyParserBackend(AbstractParserBackend, metaclass=DelegatingMeta):
    """
    Simple implementation of :class:`AbstractParserBackend` that delegates all calls to
    another parser object that not necessarely need to decend from the abstract base class.
    """
    def __init__(self, legacy_backend):
        self._delegate = legacy_backend


class LocalBackend(LegacyParserBackend):
    """
    This implementation of :class:`AbstractParserBackend` is a extended version of
    NOMAD-coe's ``LocalBackend`` that allows to write the results in an *archive*-style .json.
    It can be used like the original thing, but also allows to output archive JSON
    after parsing via :func:`write_json`.
    """
    def __init__(self, *args, **kwargs):
        delegate = LegacyLocalBackend(*args, **kwargs)
        super().__init__(delegate)
        self._status = 'none'
        self._errors = None

        self._open_context: Tuple[str, int] = None
        self._context_section = None

        # things that have no real purpos, but are required by some legacy code
        self._unknown_attributes = {}
        self.fileOut = io.StringIO()

    def __getattr__(self, name):
        """ Support for unimplemented and unexpected methods. """
        if self._unknown_attributes.get(name) is None:
            logger.debug('Access of unexpected backend attribute/method', attribute=name)
            self._unknown_attributes[name] = name
        return lambda *args, **kwargs: None

    def finishedParsingSession(self, parserStatus, parserErrors, *args, **kwargs):
        self._delegate.finishedParsingSession(parserStatus, parserErrors, *args, **kwargs)
        self._status = parserStatus
        self._errors = parserErrors

    def pwarn(self, msg):
        logger.debug('Warning in parser', parse_msg=msg)

    def _parse_context_uri(self, context_uri: str) -> Tuple[str, int]:
        """
        Returns the last segment of the given context uri, i.e. the section that
        constitutes the context.
        """
        path_str = re.sub(r'^(nmd://[^/]+/[^/]+)?/', '', context_uri, count=1)
        path = path_str.split('/')[::-1]  # reversed path via extended slice syntax

        if len(path) == 0:
            raise BadContextURI('Uri %s has not path.' % context_uri)

        while len(path) > 0:
            meta_name = path.pop()
            potential_index = path[-1] if len(path) > 0 else 'none'
            try:
                index = int(potential_index)
                path.pop()
            except ValueError:
                index = 0

        return meta_name, index

    def openSection(self, metaName: str) -> int:
        if self._open_context is None:
            return super().openSection(metaName)
        else:
            assert self._context_section is not None

            child_sections = list()

            def find_child_sections(section):
                for subsections in section.subsections.values():
                    for subsection in subsections:
                        if subsection.name == metaName:
                            child_sections.append(subsection)
                        find_child_sections(subsection)

            find_child_sections(self._context_section)

            if len(child_sections) == 0:
                return super().openSection(metaName)
            elif len(child_sections) == 1:
                index = child_sections[0].gIndex  # TODO  this also needs to be reversed, on closing sections
                self._delegate.sectionManagers[metaName].lastSectionGIndex = index
                return index
            else:
                raise WrongContextState(
                    'You cannot re-open %s with multiple instances in the context.' % metaName)

    def openContext(self, contextUri: str):
        if self._open_context is not None:
            raise WrongContextState('There is already an open context on this backend.')

        meta_name, index = self._parse_context_uri(contextUri)
        try:
            section_manager = self._delegate.sectionManagers[meta_name]
        except KeyError:
            raise BadContextURI('The section %s does not exist.' % meta_name)

        if section_manager.lastSectionGIndex < index:
            raise BadContextURI(
                'Last index of section %s is %d, cannot open %d.' %
                (meta_name, section_manager.lastSectionGIndex, index))

        self._context_section = section_manager.openSections[index]
        self._open_context = meta_name, section_manager.lastSectionGIndex
        section_manager.lastSectionGIndex = index

    def closeContext(self, contextUri):
        if self._open_context is None:
            raise WrongContextState('There is no context to close on this backend.')

        meta_name, old_index = self._open_context
        context_meta_name, _ = self._parse_context_uri(contextUri)
        if context_meta_name != meta_name:
            raise BadContextURI(
                '%d is not the URI that his context was opened with.' % contextUri)

        self._delegate.sectionManagers[context_meta_name].lastSectionGIndex = old_index
        self._open_context = None
        self._context_section = None

    @property
    def data(self) -> Results:
        return self._delegate.results

    def get_value(self, meta_name, g_index=-1):
        datamanager = self._delegate.results._datamanagers.get(meta_name)
        if datamanager is not None:
            sectionmanager = datamanager.superSectionManager
            sections = sectionmanager.openSections
            if g_index != -1:
                sections = [section for section in sections if section.gIndex == g_index]

            assert len(sections) == 1
            section = sections[0]

            return section[meta_name]

    def get_sections(self, meta_name):
        sections = self._delegate.results[meta_name]
        return [section.gIndex for section in sections]

    @staticmethod
    def _write(
            json_writer: JSONStreamWriter,
            value: Any,
            filter: Callable[[str, Any], Any]=None):

        if isinstance(value, list):
            json_writer.open_array()
            for item in value:
                LocalBackend._write(json_writer, item, filter=filter)
            json_writer.close_array()

        elif isinstance(value, Section):
            section = value
            json_writer.open_object()
            json_writer.key_value('_name', section.name)
            json_writer.key_value('_gIndex', section.gIndex)

            for name, value in section.items():
                if filter is not None:
                    value = filter(name, value)

                if value is not None:
                    json_writer.key(name)
                    LocalBackend._write(json_writer, value, filter=filter)

            json_writer.close_object()

        else:
            json_writer.value(value)

    @property
    def status(self) -> ParserStatus:
        """ Returns status and potential errors. """
        return (self._status, self._errors)

    def write_json(self, out: TextIO, pretty=True, filter: Callable[[str, Any], Any]=None):
        """
        Writes the results stored in the backend after parsing in an 'archive'.json
        style format.

        Arguments:
            out: The file-like that is used to write the json to.
            pretty: Format the json or not.
            filter: Optional filter that takes metaname, value pairs and returns a new value.
        """
        json_writer = JSONStreamWriter(out, pretty=pretty)
        json_writer.open_object()

        json_writer.key_value('parser_status', self._status)
        if self._errors is not None and len(self._errors) > 0:
            json_writer.key('parser_errors')
            json_writer.open_array
            for error in self._errors:
                json_writer.value(error)
            json_writer.close_array

        json_writer.key('section_run')
        json_writer.open_array()
        for run in self._delegate.results['section_run']:
            LocalBackend._write(json_writer, run, filter=filter)
        json_writer.close_array()

        json_writer.close_object()
        json_writer.close()

    def __repr__(self):
        def filter(name, value):
            if name.startswith('section_'):
                return value

            if name.startswith('x_'):
                return None

            if getattr(value, 'tolist', None) or isinstance(value, list):
                return '<some array>'
            else:
                return value

        out = StringIO()
        self.write_json(JSONStreamWriter(out), filter=filter)
        return out.getvalue()


class Parser():
    """
    Instances specify a parser. It allows to find *main files* from  given uploaded
    and extracted files. Further, allows to run the parser on those 'main files'.

    Arguments:
        python_git: The :class:`PythonGit` that describes the parser code.
        parser_class_name: Full qualified name of the main parser class. We assume it have one
                           parameter for the backend.
        main_file_re: A regexp that matches main file paths that this parser can handle.
        main_contents_re: A regexp that matches main file headers that this parser can parse.
    """
    def __init__(
            self, python_git: PythonGit, parser_class_name: str, main_file_re: str,
            main_contents_re: str) -> None:

        self.name = python_git.name
        self.python_git = python_git
        self.parser_class_name = parser_class_name
        self._main_file_re = re.compile(main_file_re)
        self._main_contents_re = re.compile(main_contents_re)

    def is_mainfile(self, filename: str, open: Callable[[str], IO[Any]]) -> bool:
        """ Checks if a file is a mainfile via the parsers ``main_contents_re``. """
        if self._main_file_re.match(filename):
            file = None
            try:
                file = open(filename)
                contents = file.read(500)
                return self._main_contents_re.match(contents) is not None
            finally:
                if file:
                    file.close()

        return False

    def run(self, mainfile: str) -> LocalBackend:
        """
        Runs the parser on the given mainfile. It uses :class:`LocalBackend` as
        a backend. The meta-info access is handled by the underlying NOMAD-coe parser.

        Args:
            mainfile: A path to a mainfile that this parser can parse.

        Returns:
            The used :class:`LocalBackend` with status information and result data.
        """
        def create_backend(meta_info):
            return LocalBackend(meta_info, debug=False)

        module_name = self.parser_class_name.split('.')[:-1]
        parser_class = self.parser_class_name.split('.')[1]
        module = importlib.import_module('.'.join(module_name))
        Parser = getattr(module, parser_class)

        init_signature = inspect.getargspec(Parser.__init__)
        kwargs = dict(
            backend=create_backend,
            mainfile=mainfile, main_file=mainfile,
            debug=True)
        kwargs = {key: value for key, value in kwargs.items() if key in init_signature.args}
        parser = Parser(**kwargs)

        with patch.object(sys, 'argv', []):
            backend = parser.parse(mainfile)

        # TODO we need a homogeneous interface to parsers, but we dont have it right now
        # thats a hack to distringuish between ParserInterface parser and simple_parser
        if backend is None or not hasattr(backend, 'status'):
            backend = parser.parser_context.super_backend

        return backend

    def __repr__(self):
        return self.python_git.__repr__()


parsers = [
    Parser(
        python_git=dependencies['parsers/vasp'],
        parser_class_name='vaspparser.VASPParser',
        main_file_re=r'^.*\.xml(\.[^\.]*)?$',
        main_contents_re=(
            r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
            r'?\s*<modeling>'
            r'?\s*<generator>'
            r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
            r'?')
    ),
    Parser(
        python_git=dependencies['parsers/exciting'],
        parser_class_name='parser_exciting.ExcitingParser',
        main_file_re=r'^.*/INFO\.OUT?',
        main_contents_re=(
            r'^\s*=================================================+\s*'
            r'\s*\|\s*EXCITING\s+\S+\s+started\s*='
            r'\s*\|\s*version hash id:\s*\S*\s*=')
    ),
    Parser(
        python_git=dependencies['parsers/fhi-aims'],
        parser_class_name='fhiaimsparser.FHIaimsParser',
        main_file_re=r'^.*\.out$',
        main_contents_re=(
            r'^(.*\n)*'
            r'?\s*Invoking FHI-aims \.\.\.'
            r'?\s*Version')
    )
]
""" Instanciation and constructor based config of all parsers. """

parser_dict = {parser.name: parser for parser in parsers}
""" A dict to access parsers by name. Usually 'parsers/<...>', e.g. 'parsers/vasp'. """
