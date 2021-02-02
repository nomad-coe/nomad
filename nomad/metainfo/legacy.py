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

'''
This module contains functionality to use old 'legacy' NOMAD CoE parsers with the
new nomad@fairdi infrastructure. This covers aspects like the new metainfo, a unifying
wrapper for parsers, parser logging, and a parser backend.
'''

from typing import cast, Dict, List, Any, Tuple, Type
import numpy as np
import os.path
import importlib


from nomadcore.local_meta_info import InfoKindEl, InfoKindEnv

from nomad import utils
from nomad.metainfo import (
    Definition, SubSection, Package, Quantity, Category, Section, Reference,
    Environment, MEnum, MSection, DefinitionAnnotation, MetainfoError, MSectionBound)

logger = utils.get_logger(__name__)


_ignored_packages = [
    'meta_types.nomadmetainfo.json',
    'repository.nomadmetainfo.json']


class LegacyDefinition(DefinitionAnnotation):

    def __init__(self, name: str):
        self.name = name


def def_name(definition):
    try:
        return definition.a_legacy.name
    except AttributeError:
        return definition.name


def normalize_name(name: str):
    return name.replace('.', '_').replace('-', '_')


def normalized_package_name(name: str):
    '''
    Transforms legacy metainfo '.nomadmetainfo.json' filenames into proper (python)
    identifier.
    '''
    name = name.replace('.nomadmetainfo.json', '')
    return normalize_name(name)


def python_package_mapping(metainfo_package_name: str) -> Tuple[str, str]:
    '''
    Compute the python package for the given metainfo package name. It returns
    a tuple containing a package name and a file path. The filepath denotes the file
    for this package within the nomad git project.
    '''
    prefix = metainfo_package_name.replace('.nomadmetainfo.json', '').split('.')[0]
    metainfo_package_name = normalized_package_name(metainfo_package_name)

    if prefix in ['common', 'general', 'public', 'dft', 'ems']:
        directory = 'nomad/datamodel/metainfo'
        python_package_name = 'nomad.datamodel.metainfo.%s' % metainfo_package_name

    else:
        parser_dir = prefix.replace('_', '-')
        prefix = prefix.replace('_', '')

        directory = 'dependencies/parsers/%s/%sparser/metainfo' % (parser_dir, prefix)
        python_package_name = '%sparser.metainfo.%s' % (prefix, metainfo_package_name)

    path = '%s/%s.py' % (directory, metainfo_package_name)

    return python_package_name, path


class LegacyMetainfoEnvironment(Environment):
    '''
    A metainfo environment with functions to create a legacy metainfo version of
    the environment.
    '''

    @staticmethod
    def from_legacy_package_path(path):
        metainfo_package_name = os.path.basename(path)
        package = metainfo_package_name
        if package.endswith('.nomadmetainfo.json'):
            package = package[:-19]
        if package.endswith('.json'):
            package = package[:-5]

        python_package_name, _ = python_package_mapping(package)
        python_package_name = '.'.join(python_package_name.split('.')[:-1])
        python_module = importlib.import_module(python_package_name)
        metainfo = getattr(python_module, 'm_env')

        return metainfo

    legacy_package_name = Quantity(type=str)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__section_to_sub_section_name = None
        self.__legacy_names = None

    def from_legacy_name(self, name: str, section_cls: Type[MSectionBound]) -> MSectionBound:
        ''' Returns the definition of the given globally unique legacy metainfo name. '''
        if self.__legacy_names is None:
            self.__legacy_names = dict()
            for definition in self.m_all_contents():
                try:
                    if isinstance(definition, Section):
                        if definition.extends_base_section:
                            continue
                    legacy = definition.a_legacy
                    key = (legacy.name, definition.m_def.section_cls)
                    if key in self.__legacy_names:
                        raise MetainfoError('Legacy name %s is not globally unique' % legacy.name)
                    self.__legacy_names[key] = definition
                except AttributeError:
                    pass

        return self.__legacy_names.get((name, section_cls))

    @property
    def section_to_sub_section_name(self) -> Dict[str, str]:
        if self.__section_to_sub_section_name is not None:
            return self.__section_to_sub_section_name

        self.__section_to_sub_section_name = dict()
        for definition in self.m_all_contents():
            if definition.m_def == SubSection.m_def:
                self.__section_to_sub_section_name[definition.sub_section.name] = definition.name

        return self.__section_to_sub_section_name

    def legacy_info(self, definition: Definition, *args, **kwargs) -> InfoKindEl:
        ''' Creates a legacy metainfo object for the given definition. '''
        super_names: List[str] = list()
        result: Dict[str, Any] = dict(
            name=def_name(definition),
            description=definition.description,
            superNames=super_names)

        for category in definition.categories:
            super_names.append(def_name(category))

        if isinstance(definition, Section):
            sub_section_name = self.section_to_sub_section_name.get(definition.name, definition.name)
            result['kindStr'] = 'type_section'
            result['repeats'] = any(
                sub_section.repeats
                for sub_section in self.resolve_definitions(sub_section_name, SubSection))

            for sub_section in self.resolve_definitions(sub_section_name, SubSection):
                super_names.append(def_name(sub_section.m_parent_as(Definition)))

        elif isinstance(definition, Quantity):
            result['kindStr'] = 'type_document_content'
            result['shape'] = definition.shape
            dtype_str = None
            if definition.type == int:
                dtype_str = 'i'
            elif definition.type == float:
                dtype_str = 'f'
            elif definition.type == bool:
                dtype_str = 'b'
            elif definition.type == str:
                dtype_str = 'C'
            elif isinstance(definition.type, Reference):
                dtype_str = 'r'
                result['referencedSections'] = [
                    def_name(definition.type.target_section_def.m_resolved())]
            elif isinstance(definition.type, MEnum):
                dtype_str = 'C'
            elif isinstance(definition.type, np.dtype):
                dtype_str = definition.type.name[0]
            elif definition.type == Any:
                dtype_str = 'D'
            else:
                dtype_str = str(definition.type)
                # raise TypeError(
                #     'Unsupported quantity type %s in %s.' % (definition.type, definition))
            result['dtypeStr'] = dtype_str
            if definition.unit is not None:
                result['units'] = str(definition.unit)
            super_names.append(def_name(definition.m_parent_as(Definition)))

        elif isinstance(definition, Category):
            result['kindStr'] = 'type_abstract_document_content'

        package = cast(MSection, definition)
        while not isinstance(package, Package):
            package = package.m_parent

        result['package'] = package.name

        return InfoKindEl(*args, **result, **kwargs)

    def legacy_info_env(self, packages: List[Package] = None, *args, **kwargs) -> InfoKindEnv:
        ''' Creates a legacy metainfo environment with all definitions from the given packages. '''
        if packages is None:
            packages = self.packages

        env = InfoKindEnv(*args, **kwargs)
        for package in packages:
            for definition in package.all_definitions.values():
                if not (isinstance(definition, Section) and definition.extends_base_section):
                    env.addInfoKindEl(self.legacy_info(definition))

                if isinstance(definition, Section):
                    for quantity in definition.quantities:
                        env.addInfoKindEl(self.legacy_info(quantity))

        return env

    def to_legacy_dict(
            self, packages: List[Package] = None, description: str = None,
            *args, **kwargs) -> Dict[str, Any]:
        '''
        Creates a dictionary that can be serialized to a legacy metainfo definition file
        (*.nomadmetainfo.json).

        Arguments:
            package: Will add all definitions of these packages as actual definitions,
                all other packages will be added by import.
            description: The description for the legacy file. If None the description of
                the firs package will be used.
        '''
        if packages is None:
            packages = []

        definitions = []
        dependencies = []
        for package in self.packages:
            if package in packages:
                if description is None:
                    description = package.description

                for definition in package.all_definitions.values():
                    if not (isinstance(definition, Section) and definition.extends_base_section):
                        definitions.append(self.legacy_info(definition).toDict())

                    if isinstance(definition, Section):
                        for quantity in definition.quantities:
                            definitions.append(self.legacy_info(quantity).toDict())
            else:
                dependencies.append(package)

        return {
            'type': 'nomad_meta_info_1_0',
            'description': description,
            'dependencies': [
                {'relativePath': def_name(dependency)}
                for dependency in dependencies],
            'metaInfos': definitions
        }
