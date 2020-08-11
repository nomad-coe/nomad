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

from typing import cast, Dict, List, Union, Any, Set, Iterable, Tuple, Type
import numpy as np
from pint.errors import UndefinedUnitError
import os.path
import importlib


from nomadcore.local_meta_info import loadJsonFile, InfoKindEl, InfoKindEnv

from nomad import utils
from nomad.units import ureg
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


class LegacyPackage(LegacyDefinition):
    def __init__(self, name, python_module, python_path):
        super().__init__(name)

        self.python_module = python_module
        self.python_path = python_path


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
            elif type(definition.type) == np.dtype:
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


class EnvironmentConversion:
    def __init__(self, legacy_env_or_path: Union[InfoKindEnv, str]):
        if isinstance(legacy_env_or_path, str):
            self.legacy_env, _ = loadJsonFile(filePath=legacy_env_or_path)

        else:
            self.legacy_env = cast(InfoKindEnv, legacy_env_or_path)

        self.__fix_legacy_super_names()

        self.package_conversions: Dict[str, PackageConversion] = {}

        for legacy_def in self.legacy_env.infoKindEls():
            if legacy_def.package in _ignored_packages:
                continue
            # legacy_def.package = normalized_package_name(legacy_def.package)
            package_conversion = self.package_conversions.get(legacy_def.package)
            if package_conversion is None:
                package_conversion = PackageConversion(self, legacy_def.package)
                self.package_conversions[legacy_def.package] = package_conversion

            package_conversion.legacy_defs.append(legacy_def)

        for package_conversion in self.package_conversions.values():
            package_conversion.create_definitions()

        for package_conversion in self.package_conversions.values():
            package_conversion.set_super_names()

        for package_conversion in self.package_conversions.values():
            package_conversion.init_definitions()

    def create_env(self) -> LegacyMetainfoEnvironment:
        env = LegacyMetainfoEnvironment()
        env.legacy_package_name = normalized_package_name(self.legacy_env.name)
        for package_conv in self.package_conversions.values():
            package = package_conv.package
            errors, warnings = package.m_all_validate()
            if len(errors) > 0:
                logger.error(
                    '%s. There are %d more errors in converted legacy package %s' %
                    (errors[0], len(errors) - 1, package))
            if len(warnings) > 0:
                logger.warn(
                    '%s. There are %d more warnings in converted legacy package %s' %
                    (warnings[0], len(warnings) - 1, package))
            env.m_add_sub_section(Environment.packages, package)
            package.init_metainfo()
        return env

    def __fix_legacy_super_names(self):

        def get_super_names(legacy_def: InfoKindEl, super_categories: Set[str] = None):
            super_section: str = None
            if super_categories is None:
                super_categories = set()

            for super_name in legacy_def.superNames:
                super_def = self.legacy_env.infoKindEl(super_name)

                if super_def.kindStr == 'type_section':
                    super_section = super_def.name

                elif super_def.kindStr == 'type_abstract_document_content':
                    super_categories.add(super_def.name)
                    super_super_section, _ = get_super_names(super_def, super_categories=super_categories)

                    if super_super_section is None:
                        pass

                    elif super_section is None:
                        super_section = super_super_section

                    elif super_section == super_super_section:
                        pass

                    else:
                        logger.error('conflicting parent sections %s, %s for %s' % (
                            super_section, super_def.name, legacy_def.name))

            return super_section, super_categories

        for legacy_def in self.legacy_env.infoKindEls():
            super_section, super_categories = get_super_names(legacy_def)

            if super_section is None:
                legacy_def.superNames = list(super_categories)

            else:
                legacy_def.superNames = [super_section] + list(super_categories)

    def resolve(self, name: str) -> Iterable[Definition]:
        for package_conversion in self.package_conversions.values():
            definition = package_conversion.package.all_definitions.get(name)
            if definition is not None:
                yield definition


class PackageConversion:

    def __init__(self, env_conversion: EnvironmentConversion, name: str):
        self.env_conversion = env_conversion
        self.legacy_defs: List[InfoKindEl] = []

        python_module, python_path = python_package_mapping(name)

        self.package = Package(
            name=normalize_name(name),
            a_legacy=LegacyPackage(name, python_module, python_path))

        self.quantities: Dict[str, Quantity] = {}

        self.logger = logger.bind(package=name)

    def create_definitions(self):
        for legacy_def in self.legacy_defs:
            name = normalize_name(legacy_def.name)

            if legacy_def.kindStr == 'type_abstract_document_content':
                self.package.m_create(
                    Category, name=name, a_legacy=LegacyDefinition(name=legacy_def.name))

            elif legacy_def.kindStr == 'type_section':
                self.package.m_create(
                    Section, name=name,
                    a_legacy=LegacyDefinition(name=legacy_def.name))

            elif legacy_def.kindStr in ['type_dimension', 'type_document_content']:
                definition = Quantity(
                    name=name,
                    a_legacy=LegacyDefinition(name=legacy_def.name))
                self.quantities[name] = (definition)

            else:
                logger.error('unknown kindStr %s for %s' % (legacy_def.kindStr, name))

    def __resolve(self, name: str, create_extends: bool = False):
        definition: Definition = self.package.all_definitions.get(name)
        if definition is None:
            definition = self.quantities.get(name)

        if definition is not None:
            if not (isinstance(definition, Section) and definition.extends_base_section) or create_extends:
                return definition

        for definition in self.env_conversion.resolve(name):
            if isinstance(definition, Section) and definition.extends_base_section:
                continue

            if create_extends and isinstance(definition, Section):
                extending_def = self.package.m_create(
                    Section, name=definition.name,
                    a_legacy=LegacyDefinition(name=definition.a_legacy.name))
                extending_def.base_sections = [definition]
                extending_def.extends_base_section = True
                return extending_def

            return definition

        assert False, 'definition %s must be created now' % name

    def set_super_names(self):
        for legacy_def in self.legacy_defs:
            name = normalize_name(legacy_def.name)
            definition = self.__resolve(name)
            assert definition is not None, 'definition %s must exist' % name

            if isinstance(definition, Section):
                parent_section: Section = None
                for super_name in legacy_def.superNames:
                    super_def = self.__resolve(normalize_name(super_name), create_extends=True)
                    if isinstance(super_def, Section):
                        parent_section = cast(Section, super_def)

                if parent_section is not None:
                    sub_section = parent_section.m_create(
                        SubSection, name=definition.name,
                        a_legacy=LegacyDefinition(name=legacy_def.name))
                    sub_section.sub_section = definition
                    sub_section.repeats = legacy_def.repeats is None or legacy_def.repeats

            if isinstance(definition, Quantity):
                parent_section: Section = None
                for super_name in legacy_def.superNames:
                    super_def = self.__resolve(normalize_name(super_name), create_extends=True)
                    if isinstance(super_def, Section):
                        parent_section = cast(Section, super_def)

                parent_section.m_add_sub_section(Section.quantities, definition)

    def init_definitions(self):
        for legacy_def in self.legacy_defs:
            name = normalize_name(legacy_def.name)
            definition = self.__resolve(name)
            assert definition is not None, 'definition %s must exist' % name
            logger = self.logger.bind(definition=definition.name)

            # common properties
            if legacy_def.description is not None and legacy_def.description.strip() != '':
                definition.description = legacy_def.description

            if isinstance(definition, Definition):
                # deal with categories
                categories: List[Category] = []
                for super_name in legacy_def.superNames:
                    super_def = self.__resolve(normalize_name(super_name))
                    if isinstance(super_def, Category):
                        categories.append(super_def)

                definition.categories = categories

            if isinstance(definition, Quantity):
                # type
                referenced_sections = legacy_def.extra_args.get('referencedSections', [])
                if len(referenced_sections) == 1:
                    referenced_section = self.__resolve(referenced_sections[0])
                    if referenced_section is None:
                        logger.error('could not find referencedSection %s of %s' % (
                            referenced_sections[0], definition.name))
                        definition.type = int
                    else:
                        definition.type = Reference(referenced_section)

                elif len(referenced_sections) > 1:
                    logger.error(
                        'higher dimensional references not yet supported: %s' % name)
                    definition.type = np.dtype(int)

                elif legacy_def.kindStr == 'type_dimension':
                    definition.type = int
                elif legacy_def.dtypeStr == 'D':
                    definition.type = Any
                elif legacy_def.dtypeStr == 'C':
                    definition.type = str
                elif legacy_def.dtypeStr == 'r':
                    logger.error('r typed quantity %s  doesn\'t have referencedSections' % name)
                    definition.type = int
                elif legacy_def.dtypeStr == 'b':
                    definition.type = bool
                elif legacy_def.dtypeStr == 'i64':
                    definition.type = np.dtype(np.int64)
                elif legacy_def.dtypeStr == 'f':
                    definition.type = np.dtype(np.float64)
                else:
                    definition.type = np.dtype(legacy_def.dtypeStr)

                # shapes
                legacy_shape = legacy_def.shape
                if legacy_shape is None:
                    legacy_shape = []

                definition.shape = legacy_shape
                if len(definition.shape) > 1 and definition.type == str:
                    # Usually only np types have higher shapes in old metainfo;
                    # str is one exception.
                    definition.type = np.dtype('U')

                # units
                if legacy_def.units is not None:
                    try:
                        definition.unit = ureg.parse_units(legacy_def.units)
                    except UndefinedUnitError:
                        logger.error('unknown unit %s' % legacy_def.units)
                    except ValueError as e:
                        logger.error('cannot parse unit %s' % legacy_def.units, exc_info=e)


def convert(metainfo_path: str) -> LegacyMetainfoEnvironment:
    return EnvironmentConversion(metainfo_path).create_env()


def generate_metainfo_code(metainfo_env: LegacyMetainfoEnvironment):
    '''
    Generates python code with metainfo definitions for all packages in the given
    environement

    Arguments:
        env: The metainfo environment.
        python_package_path: An optional directory path. The directory must exist. Default
            is the working directory. The path will be used to form the module prefix
            for generated Python modules.
    '''
    from jinja2 import Environment as JinjaEnvironment, PackageLoader, select_autoescape
    import textwrap

    def format_description(description, indent=0, width=90):
        paragraphs = [paragraph.strip() for paragraph in description.split('\n')]

        def format_paragraph(paragraph, first):
            lines = textwrap.wrap(text=paragraph, width=width - indent * 4)
            lines = [line.replace('\\', '\\\\') for line in lines]
            return textwrap.indent(
                '\n'.join(lines), ' ' * 4 * indent, lambda x: not (first and x.startswith(lines[0])))

        return '\n\n'.join([
            format_paragraph(p, i == 0)
            for i, p in enumerate(paragraphs) if p != ''])

    def format_type(pkg, mi_type):
        if type(mi_type) == np.dtype:
            if mi_type == np.dtype('U'):
                return 'np.dtype(\'U\')'

            return 'np.dtype(np.%s)' % mi_type

        if mi_type in [int, float, str, bool]:
            return mi_type.__name__

        if isinstance(mi_type, Reference):
            if pkg == mi_type.target_section_def.m_parent:
                return "Reference(SectionProxy('%s'))" % mi_type.target_section_def.name

            else:
                python_module = mi_type.target_section_def.m_parent.a_legacy.python_module
                return '%s.%s' % (python_module.split('.')[-1], mi_type.target_section_def.name)

        else:
            return str(mi_type)

    def format_unit(unit):
        return "'%s'" % unit

    def format_definition_refs(pkg, definitions):
        def format_definition_ref(definition: Definition):
            if pkg == definition.m_parent:
                return definition.name
            else:
                python_module = definition.m_parent.a_legacy.python_module
                return '%s.%s' % (python_module.split('.')[-1], definition.name)

        return ', '.join([format_definition_ref(definition) for definition in definitions])

    def fromat_package_import(pkg):
        python_module = pkg.a_legacy.python_module
        modules = python_module.split('.')
        return 'from %s import %s' % ('.'.join(modules[:-1]), modules[-1])

    def order_categories(categories):
        return sorted(categories, key=lambda c: len(c.categories))

    env = JinjaEnvironment(
        loader=PackageLoader('nomad.metainfo', 'templates'),
        autoescape=select_autoescape(['python']))
    env.globals.update(
        order_categories=order_categories,
        format_description=format_description,
        format_type=format_type,
        format_unit=format_unit,
        format_definition_refs=format_definition_refs,
        fromat_package_import=fromat_package_import)

    for package in metainfo_env.packages:
        path = package.a_legacy.python_path
        if not os.path.exists(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))

        with open(path, 'wt') as f:
            code = env.get_template('package.j2').render(pkg=package)
            code = '\n'.join([
                line.rstrip() if line.strip() != '' else ''
                for line in code.split('\n')])
            f.write(code)

    _, path = python_package_mapping(metainfo_env.legacy_package_name)
    with open(os.path.join(os.path.dirname(path), '__init__.py'), 'wt') as f:

        code = env.get_template('environment.j2').render(env=metainfo_env)
        code = '\n'.join([
            line.rstrip() if line.strip() != '' else ''
            for line in code.split('\n')])
        f.write(code)
