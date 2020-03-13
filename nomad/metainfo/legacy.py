
from typing import cast, Dict, List, Union, Any, Set, Iterable
import numpy as np
from pint.errors import UndefinedUnitError
import os.path
from jinja2 import Environment as JinjaEnvironment, PackageLoader, select_autoescape
import textwrap

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl, InfoKindEnv
import nomad_meta_info

from nomad import utils
from nomad.metainfo import (
    Definition, SubSection, Package, Quantity, Category, Section, Reference, units,
    Environment, MEnum)


logger = utils.get_logger(__name__)


_ignored_packages = [
    'meta_types.nomadmetainfo.json',
    'repository.nomadmetainfo.json']


class LegacyMetainfoEnvironment(Environment):
    '''
    A metainfo environment with functions to create a legacy metainfo version of
    the environment.
    '''
    def legacy_info(self, definition: Definition, *args, **kwargs) -> InfoKindEl:
        ''' Creates a legacy metainfo objects for the given definition. '''
        super_names: List[str] = list()
        result: Dict[str, Any] = dict(
            name=definition.name,
            description=definition.description,
            superNames=super_names)

        for category in definition.categories:
            super_names.append(category.name)

        if isinstance(definition, Section):
            result['kindStr'] = 'type_section'
            result['repeats'] = any(
                sub_section.repeats
                for sub_section in self.resolve_definitions(definition.name, SubSection))

            for sub_section in self.resolve_definitions(definition.name, SubSection):
                super_names.append(sub_section.m_parent_as(Definition).name)

        elif isinstance(definition, Quantity):
            result['kindStr'] = 'document_content'
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
                result['referencedSections'] = [definition.type.target_section_def.name]
            elif isinstance(definition.type, MEnum):
                dtype_str = 'C'
            elif type(definition.type) == np.dtype:
                dtype_str = definition.type.name[0]
            elif definition.type == Any:
                dtype_str = 'D'
            else:
                raise TypeError(
                    'Unsupported quantity type %s in %s.' % (definition.type, definition))
            result['dtypeStr'] = dtype_str
            if definition.unit is not None:
                result['units'] = str(definition.unit)
            super_names.append(definition.m_parent_as(Definition).name)

        elif isinstance(definition, Category):
            result['kindStr'] = 'abstract_document_content'

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


class EnvironmentConversion:
    def __init__(self, legacy_env_or_path: Union[InfoKindEnv, str]):
        if isinstance(legacy_env_or_path, str):
            legacy_env_or_path = os.path.normpath(os.path.join(
                os.path.dirname(nomad_meta_info.__file__), legacy_env_or_path))
            self.legacy_env, _ = loadJsonFile(filePath=legacy_env_or_path)

        else:
            self.legacy_env = cast(InfoKindEnv, legacy_env_or_path)

        self.__fix_legacy_super_names()

        self.package_conversions: Dict[str, PackageConversion] = {}

        for legacy_def in self.legacy_env.infoKindEls():
            if legacy_def.package in _ignored_packages:
                continue
            legacy_def.package = legacy_def.package.replace('.nomadmetainfo.json', '').replace('.', '_')
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

        self.package = Package(name=name)
        self.quantities: Dict[str, Quantity] = {}

        self.logger = logger.bind(package=name)

    def create_definitions(self):
        for legacy_def in self.legacy_defs:
            name = legacy_def.name

            if legacy_def.kindStr == 'type_abstract_document_content':
                self.package.m_create(Category, name=name)

            elif legacy_def.kindStr == 'type_section':
                self.package.m_create(Section, name=name)

            elif legacy_def.kindStr in ['type_dimension', 'type_document_content']:
                definition = Quantity(name=name)
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
                extending_def = self.package.m_create(Section, name=definition.name)
                extending_def.base_sections = [definition]
                extending_def.extends_base_section = True
                return extending_def

            return definition

        assert False, 'definition %s must be created now' % name

    def set_super_names(self):
        for legacy_def in self.legacy_defs:
            definition = self.__resolve(legacy_def.name)
            assert definition is not None, 'definition %s must exist' % legacy_def.name

            if isinstance(definition, Section):
                parent_section: Section = None
                for super_name in legacy_def.superNames:
                    super_def = self.__resolve(super_name, create_extends=True)
                    if isinstance(super_def, Section):
                        parent_section = cast(Section, super_def)

                if parent_section is not None:
                    sub_section = parent_section.m_create(SubSection, name=definition.name)
                    sub_section.sub_section = definition
                    sub_section.repeats = legacy_def.repeats is None or legacy_def.repeats

            if isinstance(definition, Quantity):
                parent_section: Section = None
                for super_name in legacy_def.superNames:
                    super_def = self.__resolve(super_name, create_extends=True)
                    if isinstance(super_def, Section):
                        parent_section = cast(Section, super_def)

                parent_section.m_add_sub_section(Section.quantities, definition)

    def init_definitions(self):
        for legacy_def in self.legacy_defs:
            definition = self.__resolve(legacy_def.name)
            assert definition is not None, 'definition %s must exist' % legacy_def.name
            logger = self.logger.bind(definition=definition.name)

            # common properties
            definition.description = legacy_def.description

            if isinstance(definition, Definition):
                # deal with categories
                categories: List[Category] = []
                for super_name in legacy_def.superNames:
                    super_def = self.__resolve(super_name)
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
                        'higher dimensional references not yet supported: %s' % legacy_def.name)
                    definition.type = np.dtype(int)

                elif legacy_def.kindStr == 'type_dimension':
                    definition.type = int
                elif legacy_def.dtypeStr == 'D':
                    definition.type = Any
                elif legacy_def.dtypeStr == 'C':
                    definition.type = str
                elif legacy_def.dtypeStr == 'r':
                    logger.error('r typed quantity %s  doesn\'t have referencedSections' % legacy_def.name)
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
                        definition.unit = units.parse_units(legacy_def.units)
                    except UndefinedUnitError:
                        logger.error('unknown unit %s' % legacy_def.units)
                    except ValueError as e:
                        logger.error('cannot parse unit %s' % legacy_def.units, exc_info=e)


def convert(metainfo_path: str) -> LegacyMetainfoEnvironment:
    return EnvironmentConversion(metainfo_path).create_env()


def generate_metainfo_code(metainfo_env: Environment, directory: str = None):
    '''
    Generates python code with metainfo definitions for all packages in the given
    environement

    Arguments:
        env: The metainfo environment.
        directory: An optional directory path. The directory must exist. Default
            is the working directory.
    '''

    if directory is None:
        directory = '.'

    def format_description(description, indent=0, width=90):
        paragraphs = [paragraph.strip() for paragraph in description.split('\n')]

        def format_paragraph(paragraph, first):
            lines = textwrap.wrap(text=paragraph, width=width - indent * 4)
            lines = [l.replace('\\', '\\\\') for l in lines]
            return textwrap.indent(
                '\n'.join(lines), ' ' * 4 * indent, lambda x: not (first and x.startswith(lines[0])))

        return '\n\n'.join([
            format_paragraph(p, i == 0)
            for i, p in enumerate(paragraphs) if p != ''])

    def format_type(mi_type):
        if type(mi_type) == np.dtype:
            return 'np.dtype(np.%s)' % mi_type
        if mi_type in [int, float, str, bool]:
            return mi_type.__name__
        if isinstance(mi_type, Reference):
            return "SectionProxy('%s')" % mi_type.target_section_def.name
        else:
            return str(mi_type)

    def format_unit(unit):
        return "'%s'" % unit

    def format_definition_refs(pkg, definitions):
        def format_definition_ref(definition: Definition):
            if pkg == definition.m_parent:
                return definition.name
            else:
                return definition.qualified_name()

        return ', '.join([format_definition_ref(definition) for definition in definitions])

    env = JinjaEnvironment(
        loader=PackageLoader('nomad.metainfo', 'templates'),
        autoescape=select_autoescape(['python']))
    env.globals.update(
        format_description=format_description,
        format_type=format_type,
        format_unit=format_unit,
        format_definition_refs=format_definition_refs)

    for package in metainfo_env.packages:
        file_name = package.name
        with open(os.path.join(directory, '%s.py' % file_name), 'wt') as f:
            code = env.get_template('package.j2').render(pkg=package)
            code = '\n'.join([
                line.rstrip() if line.strip() != '' else ''
                for line in code.split('\n')])
            f.write(code)


# if __name__ == '__main__':
#     output = 'output'

#     env = convert('vasp.nomadmetainfo.json')
#     assert env.resolve_definition('x_vasp_incar_EFIELD_PEAD', Quantity) is not None
#     assert 'x_vasp_incar_EFIELD_PEAD' in env.legacy_info_env()
#     generate_metainfo_code(env, output)

#     from output import public
#     import json

#     run = public.section_run()
#     system = run.m_create(public.section_system)
#     system.atom_labels = ['H', 'H', 'O']

#     print(json.dumps(run.m_to_dict(with_meta=True), indent=2))
