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

'''
This module contains functionality to use old 'legacy' NOMAD CoE parsers with the
new nomad@fairdi infrastructure. This covers aspects like the new metainfo, a unifying
wrapper for parsers, parser logging, and a parser backend.
'''

from typing import List, Union, Any, Tuple, Type, cast
from abc import ABCMeta, abstractmethod
import importlib
import os.path
import inspect
from unittest.mock import patch
import logging
import glob
import sys

from nomad import utils, datamodel, config
from nomad.metainfo import (
    SubSection, Quantity, Section, SectionProxy, Reference, MResource, MSection, MSectionBound, Property)
from nomad.metainfo.legacy import (
    LegacyMetainfoEnvironment, python_package_mapping, normalize_name)

from .parser import MatchingParser


class BackendError(Exception):
    pass


class AbstractParserBackend(metaclass=ABCMeta):
    '''
    This ABS provides the parser backend interface used by the NOMAD-coe parsers.
    '''

    @abstractmethod
    def metaInfoEnv(self):
        ''' Returns the meta info used by this backend. '''
        pass

    @abstractmethod
    def startedParsingSession(
            self, mainFileUri, parserInfo, parserStatus=None, parserErrors=None):
        '''
        Should be called when the parsing starts.
        ParserInfo should be a valid json dictionary.
        '''
        pass

    @abstractmethod
    def finishedParsingSession(
            self, parserStatus, parserErrors, mainFileUri=None, parserInfo=None,
            parsingStats=None):
        ''' Called when the parsing finishes. '''
        pass

    @abstractmethod
    def openSection(self, metaName, parent_index=-1):
        ''' Opens a new section and returns its new unique gIndex. '''
        pass

    @abstractmethod
    def closeSection(self, metaName, gIndex):
        '''
        Closes the section with the given meta name and index. After this, no more
        value can be added to this section.
        '''
        pass

    @abstractmethod
    def openNonOverlappingSection(self, metaName):
        ''' Opens a new non overlapping section. '''
        pass

    @abstractmethod
    def setSectionInfo(self, metaName, gIndex, references):
        '''
        Sets info values of an open section references should be a dictionary with the
        gIndexes of the root sections this section refers to.
        '''
        pass

    @abstractmethod
    def closeNonOverlappingSection(self, metaName):
        '''
        Closes the current non overlapping section for the given meta name. After
        this, no more value can be added to this section.
        '''
        pass

    @abstractmethod
    def openSections(self):
        ''' Returns the sections that are still open as metaName, gIndex tuples. '''
        pass

    @abstractmethod
    def addValue(self, metaName, value, gIndex=-1):
        '''
        Adds a json value for the given metaName. The gIndex is used to identify
        the right parent section.
        '''
        pass

    @abstractmethod
    def addRealValue(self, metaName, value, gIndex=-1):
        '''
        Adds a float value for the given metaName. The gIndex is used to identify
        the right parent section.
        '''
        pass

    @abstractmethod
    def addArray(self, metaName, shape, gIndex=-1):
        '''
        Adds an unannitialized array of the given shape for the given metaName.
        The gIndex is used to identify the right parent section.
        This is neccessary before array values can be set with :func:`setArrayValues`.
        '''

    @abstractmethod
    def setArrayValues(self, metaName, values, offset=None, gIndex=-1):
        '''
        Adds values of the given numpy array to the last array added for the given
        metaName and parent gIndex.
        '''
        pass

    @abstractmethod
    def addArrayValues(self, metaName, values, gIndex=-1, override: bool = False):
        '''
        Adds an array with the given numpy array values for the given metaName and
        parent section gIndex. Override determines whether to rewrite exisiting values
        in the backend.
        '''
        pass

    @abstractmethod
    def pwarn(self, msg):
        ''' Used to catch parser warnings. '''
        pass

    @abstractmethod
    def get_sections(self, meta_name: str, g_index: int = -1) -> List[int]:
        ''' Return all gIndices for existing sections of the given meta_name and parent section index. '''
        pass

    @abstractmethod
    def get_value(self, metaName: str, g_index=-1) -> Any:
        '''
        Return the value set to the given meta_name in its parent section of the given index.
        An index of -1 (default) is only allowed if there is exactly one parent section.
        '''
        pass

    @abstractmethod
    def __getitem__(self, key):
        pass


class Backend(AbstractParserBackend):
    '''
    A backend that uses the new metainfo to store all data.

    Arguments:
        metainfo: The main legacy metainfo package name or a legacy metainfo environment
            instance.
        logger: An optional logger.
        domain: An optional domain name.

    Attributes:
        domain: The domain that this backend contains data for.
        env: The metainfo environment (all available definitions).
        resource: The metainfo resource that contains all data.
        entry_archive: The root section of the archive behind this backend.
        logger: A logger that can be used to log metainfo and backend operation related
            warnings and errors.
    '''

    @staticmethod
    def filter_legacy_defs(definition):
        return definition.m_get_annotations('legacy') is not None

    def __init__(
            self, metainfo: Union[str, LegacyMetainfoEnvironment], domain: str = None,
            logger=None):
        assert metainfo is not None

        if logger is None:
            logger = utils.get_logger(__name__)
        self.logger = logger
        self.domain = domain if domain is not None else 'dft'  # TODO

        if isinstance(metainfo, str):
            python_package_name, _ = python_package_mapping(metainfo)
            python_package_name = '.'.join(python_package_name.split('.')[:-1])
            python_module = importlib.import_module(python_package_name)
            metainfo = getattr(python_module, 'm_env')

        self.env: LegacyMetainfoEnvironment = cast(LegacyMetainfoEnvironment, metainfo)
        self.__legacy_env = None
        self.resource = MResource(logger=logger)
        self.entry_archive = datamodel.EntryArchive()
        self.resource.add(self.entry_archive)

        self.strict = False  # TODO

        self.reset_status()
        # things that have no real purpose, but are required by some legacy code
        # self._unknown_attributes = {}
        # self._known_attributes = ['results']
        # self.fileOut = io.StringIO()

    def __open_section(self, section_def: Union[Section, SectionProxy], index: int = -1) -> MSection:
        ''' Returns the last opened section with the given parent index. '''
        sections = self.resource.all(section_def.section_cls)

        if sections is None:
            raise KeyError(
                'section %s with given parent index %s does not exist' % (section_def.name, index))

        if index == -1:
            return sections[-1]

        for section in reversed(sections):
            if section.m_parent_index == index:
                return section

        raise KeyError(
            'section %s with given parent index %s does not exist (B)' % (section_def.name, index))

    def __getitem__(self, key):
        property_def = self.resolve_definition(key, Property)
        section_def = property_def.m_parent
        if section_def.extends_base_section:
            section_def = section_def.base_sections[0]

        section = self.__open_section(section_def)
        if isinstance(property_def, Quantity):
            return section.m_get(property_def)

        elif isinstance(property_def, SubSection):
            return section.m_get_sub_sections(property_def)

    def metaInfoEnv(self):
        if self.__legacy_env is None:
            self.__legacy_env = self.env.legacy_info_env()
        return self.__legacy_env

    def resolve_definition(self, name, section_cls: Type[MSectionBound]) -> MSectionBound:
        definition = self.env.from_legacy_name(name, section_cls)
        if definition:
            return definition

        return self.env.resolve_definition(
            normalize_name(name), section_cls, Backend.filter_legacy_defs)

    def openSection(self, name, parent_index: int = -1, return_section=False):
        '''
        It will assume that there is a sub-section def with the given name.
        It will use the latest opened section of the sub-sections parent as the parent
        for the new section.
        An Exception will be known root sections, e.g. 'section_run'.

        Args:
            name: The section name
            parent_index: Index of the parent section in which the section will
                be opened in.
            return_section: If True, returns the section instead of the section index.
        '''
        section_def = self.resolve_definition(name, Section)

        if name in datamodel.root_sections:
            section = self.entry_archive.m_create(section_def.section_cls)

        else:
            sub_section_def = self.resolve_definition(name, SubSection)
            parent_section_def = sub_section_def.m_parent_as(Section)
            if parent_section_def.extends_base_section:
                parent_section_def = parent_section_def.base_sections[0]

            parent = self.__open_section(parent_section_def, parent_index)
            section = parent.m_create(section_def.section_cls, sub_section_def)

        if return_section:
            return section
        return section.m_parent_index

    def get_open_section_for_quantity(self, name, g_index):
        ''' Returns the open section that contains the quantity of the given name. '''

        quantity_def = self.resolve_definition(name, Quantity)
        section_def = quantity_def.m_parent_as(Section)
        if section_def.extends_base_section:
            section_def = section_def.base_sections[0]

        section = self.__open_section(section_def, g_index)

        return section, quantity_def

    def closeSection(self, name, g_index):
        pass

    def openNonOverlappingSection(self, metaName):
        return self.openSection(metaName)

    def setSectionInfo(self, metaName, gIndex, references):
        '''
        Sets info values of an open section references should be a dictionary with the
        gIndexes of the root sections this section refers to.
        '''
        # TODO might be necessary to make references work?
        pass

    def closeNonOverlappingSection(self, name):
        return self.closeSection(name, -1)

    def openSections(self):
        ''' Returns the sections that are still open as metaName, gIndex tuples. '''
        raise NotImplementedError()
        # for section_def, sub_sections in self.open_sections_by_def:
        #     for sub_section in sub_sections:
        #         yield section_def.name, sub_section.m_parent_index

    def addValue(self, name, value, g_index=-1):
        section, quantity_def = self.get_open_section_for_quantity(name, g_index)
        if isinstance(quantity_def.type, Reference):
            # quantity is a reference
            try:
                referenced_target = self.__open_section(
                    quantity_def.type.target_section_def, value)
            except KeyError as e:
                raise BackendError(str(e))

            value = referenced_target

        setattr(section, name, value)

    def addRealValue(self, name, value, g_index=-1):
        self.addValue(name, value, g_index)

    def addArray(self, name, shape, g_index=-1):
        '''
        Adds an uninitialized array of the given shape for the given metaName.
        The gIndex is used to identify the right parent section.
        This is neccessary before array values can be set with :func:`setArrayValues`.
        '''
        raise NotImplementedError()

    def setArrayValues(self, metaName, values, offset=None, gIndex=-1):
        '''
        Adds values of the given numpy array to the last array added for the given
        metaName and parent gIndex.
        '''
        raise NotImplementedError()

    def addArrayValues(self, name, values, gIndex=-1, override: bool = False):
        '''
        Adds an array with the given numpy array values for the given metaName and
        parent section gIndex. Override determines whether to rewrite exisiting values
        in the backend.
        '''
        section, quantity_def = self.get_open_section_for_quantity(name, gIndex)
        if isinstance(quantity_def.type, Reference):
            # quantity is a reference
            resolved_values = []
            for value in values:
                # quantity is a reference
                try:
                    referenced_target = self.__open_section(quantity_def.type.target_section_def, value)
                except KeyError:
                    raise BackendError('There is not section for the given reference index')

                resolved_values.append(referenced_target)
            values = resolved_values

        if self.strict and not override:
            quantity_def = section.m_def.all_quantities[name]
            assert not section.m_is_set(quantity_def)

        setattr(section, name, values)

    # The following are extensions to the origin NOMAD-coe parser backend. And allow
    # access to existing data

    def get_sections(self, meta_name: str, g_index: int = -1) -> List[int]:
        ''' Return all gIndices for existing sections of the given meta_name and parent index. '''
        section_def = self.resolve_definition(meta_name, Section)
        return [
            section.m_parent_index for section in self.resource.all(section_def.section_cls)
            if g_index == -1 or section.m_parent.m_parent_index == g_index]

    def get_value(self, meta_name: str, g_index=-1) -> Any:
        '''
        Return the value set to the given meta_name in its parent section of the given index.
        An index of -1 (default) is only allowed if there is exactly one parent section.
        '''
        section, quantity_def = self.get_open_section_for_quantity(meta_name, g_index)
        value = section.m_get(quantity_def)

        if value is None:
            raise KeyError(quantity_def.name)

        if isinstance(quantity_def.type, Reference):
            return value.m_parent_index

        return value

    def startedParsingSession(
            self, mainFileUri, parserInfo, parserStatus=None, parserErrors=None):
        self.reset_status()

    def finishedParsingSession(
            self, parserStatus, parserErrors, mainFileUri=None, parserInfo=None,
            parsingStats=None):
        self._status = parserStatus
        self._errors = parserErrors

    def addMatchTelemetry(self, match_telemetry, gIndex=-1):
        # TODO
        pass

    def pwarn(self, msg):
        self.logger.warn(msg)
        if len(self._warnings) < 10:
            self._warnings.append(msg)
        elif len(self._warnings) == 10:
            self._warnings.append('There are more warnings, check the processing logs.')

    @property
    def status(self) -> Tuple[str, List[str]]:
        ''' Returns status and potential errors. '''
        return (self._status, self._errors)

    def reset_status(self) -> None:
        self._status = 'ParseSuccess'
        self._errors = None
        self._warnings: List[str] = []


class LegacyParser(MatchingParser):
    '''
    A parser implementation for legacy NOMAD-coe parsers. It assumes that parsers
    are installed to the python environment.

    Arguments:
        parser_class_name: the main parser class that implements NOMAD-coe's
        backend_factory: a callable that returns a backend, takes meta_info and logger as argument
    '''
    def __init__(self, parser_class_name: str, *args, backend_factory=None, **kwargs) -> None:
        super().__init__(*args, **kwargs)

        self.parser_class_name = parser_class_name
        self.backend_factory = backend_factory

        module_name = self.parser_class_name.split('.')[:-1]
        parser_class_name = self.parser_class_name.split('.')[-1]
        self.__parser_impl = module_name, parser_class_name
        self.__parser_class = None

    @property
    def metainfo_env(self):
        if self._metainfo_env is None:
            module_name, _ = self.__parser_impl
            module = importlib.import_module('.'.join(module_name + ['metainfo']))
            self._metainfo_env = getattr(module, 'm_env')

        return super().metainfo_env

    @property
    def parser_class(self):
        if self.__parser_class is None:
            module_name, parser_class_name = self.__parser_impl
            module = importlib.import_module('.'.join(module_name))
            self.__parser_class = getattr(module, parser_class_name)

        return self.__parser_class

    def run(self, mainfile: str, logger=None) -> Backend:
        # TODO we need a homogeneous interface to parsers, but we dont have it right now.
        # There are some hacks to distinguish between ParserInterface parser and simple_parser
        # using hasattr, kwargs, etc.
        if logger is None:
            logger = utils.get_logger('__name__')

        if issubclass(self.parser_class, CoEParser):
            # TODO reuse parser
            parser = self.parser_class()  # pylint: disable=not-callable
            return parser.run(mainfile, logger)

        def create_backend(meta_info):
            if self.backend_factory is not None:
                return self.backend_factory(meta_info, logger=logger)

            return Backend(meta_info, logger=logger, domain=self.domain)

        init_signature = inspect.getargspec(self.parser_class.__init__)
        kwargs = dict(backend=create_backend, log_level=logging.DEBUG, debug=True)
        kwargs = {key: value for key, value in kwargs.items() if key in init_signature.args}

        with utils.legacy_logger(logger):
            self.parser = self.parser_class(**kwargs)  # pylint: disable=not-callable

            with patch.object(sys, 'argv', []):
                backend = self.parser.parse(mainfile)
                os.chdir(config.fs.working_directory)

            if backend is None or not hasattr(backend, 'status'):
                backend = self.parser.parser_context.super_backend

        return backend


class CoEParser(metaclass=ABCMeta):

    @abstractmethod
    def run(self, mainfile, logger) -> Backend:
        pass


class CoEInterfaceParser(CoEParser):

    def __init__(self, interface_class):
        super().__init__()
        self.interface_class = interface_class
        self.__interface = None

    def run(self, mainfile, logger):
        if self.__interface is None:
            self.__interface = self.interface_class()

        self.__interface.setup_logger(logger)
        self.__interface.parse(mainfile)
        return self.__interface.parser_context.super_backend


class CoESimpleMatcherParser(CoEParser):

    def __init__(self):
        super().__init__()
        self.parser_description = self.create_parser_description()
        self.simple_matcher = self.create_simple_matcher()

        self._metainfo_env = self.metainfo_env()
        self.__legacy_metainfo_env = None

        self.caching_levels = self.create_caching_levels()

    @abstractmethod
    def metainfo_env(self) -> LegacyMetainfoEnvironment:
        pass

    @property
    def metaInfoEnv(self) -> LegacyMetainfoEnvironment:
        if self.__legacy_metainfo_env is None:
            self.__legacy_metainfo_env = self._metainfo_env.legacy_info_env()
        return self.__legacy_metainfo_env

    def create_caching_levels(self) -> dict:
        return dict()

    @abstractmethod
    def create_simple_matcher(self):
        pass

    @abstractmethod
    def create_parser_description(self) -> dict:
        raise NotImplementedError()

    @abstractmethod
    def create_super_context(self):
        pass

    def simple_parser(self, mainfile, logger) -> Backend:
        from nomadcore.simple_parser import mainFunction
        backend = Backend(self._metainfo_env, logger=logger)
        from unittest.mock import patch
        with patch.object(sys, 'argv', ['<exe>', '--uri', 'nmd://uri', mainfile]):
            mainFunction(
                mainFileDescription=self.simple_matcher,
                metaInfoEnv=self.metaInfoEnv,
                parserInfo=self.parser_description,
                cachingLevelForMetaName=self.caching_levels,
                superContext=self.create_super_context(),
                superBackend=backend)

        return backend

    def run(self, mainfile, logger) -> Backend:
        with utils.legacy_logger(logger):
            return self.simple_parser(mainfile, logger)


class VaspOutcarParser(LegacyParser):
    '''
    LegacyParser that only matches mailfiles, if there is no .xml in the
    same directory, i.e. to use the VASP OUTCAR parser in absence of .xml
    output file.
    '''
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
