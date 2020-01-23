from typing import Tuple, Dict, List, Type, TypeVar, Any, Union, cast
import os.path
import numpy as np
from jinja2 import Environment as JinjaEnvironment, PackageLoader, select_autoescape
import textwrap

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl, InfoKindEnv
import nomad_meta_info

from nomad import utils
from nomad.metainfo import Definition, Package, Category, Section, Quantity, SubSection, \
    Environment, MEnum, Reference, MSection, units


T = TypeVar('T', bound=Definition)


def from_legacy_metainfo(meta_info_env, package_names: List[str] = None) \
        -> Tuple[Dict[str, InfoKindEl], Dict[str, List[InfoKindEl]]]:
    defs: Dict[str, InfoKindEl] = {}
    packages: Dict[str, List[InfoKindEl]] = {}
    for definition in meta_info_env.infoKindEls():
        if package_names is None or definition.package in package_names:
            defs[definition.name] = definition
            packages.setdefault(definition.package, []).append(definition)

    # resolve implicit super section
    # Some defintions in the old metainfo only have abstract document content as super
    # names, and thereby implicitely have a super section defined as super name of one
    # of their super abstract document contents
    def resolve_super_section(name) -> str:
        definition = defs.get(name)
        assert definition is not None, '%s could not be resolved' % name
        for super_name in definition.superNames:
            if super_name.startswith('section_'):
                return super_name

        for super_name in definition.superNames:
            super_section = resolve_super_section(super_name)
            if super_section is not None:
                definition.superNames.append(super_section)
                return super_section

        return None

    for definition in defs.values():
        resolve_super_section(definition.name)

    return defs, packages


class LegacyMetainfoEnvironment:
    """
    Args:
        env: The metainfo environment that is used to manage the definitions.
        orig_legacy_env: The old metainfo :class:`InfoKindEnv` environment with the
            original input.
        legacy_env: A old metainfo :class:`InfoKindEnv` environment generated from the
            converted metainfo environment.
        all_legacy_defs: A dict that stores the original :class:`InfoKindEl`s by name.
        all_defs: A dict that stroed the converted section and category definitions.
    """
    def __init__(self, metainfo=Union[InfoKindEnv, str], package_names: List[str] = None, logger=None):
        self.logger = utils.get_logger(__name__) if logger is None else logger
        self.env = Environment()

        if isinstance(metainfo, InfoKindEnv):
            self.orig_legacy_env = metainfo
        else:
            meta_info_path = os.path.normpath(os.path.join(
                os.path.dirname(nomad_meta_info.__file__), metainfo))
            self.orig_legacy_env, _ = loadJsonFile(filePath=meta_info_path)

        self.all_legacy_defs, self.legacy_packages = from_legacy_metainfo(
            self.orig_legacy_env, package_names=package_names)
        self.all_defs: Dict[str, Definition] = dict()

        self.__defs_in_flux: Dict[str, Definition] = dict()
        for legacy_package_name, legacy_package in self.legacy_packages.items():
            if legacy_package_name not in ['meta_types.nomadmetainfo.json']:
                self.convert_package(legacy_package, name=legacy_package_name)

        if len(self.__defs_in_flux) > 0:
            self.logger.error(
                '%s is referenced by not defined. There are %d more like this.' %
                (list(self.__defs_in_flux.keys())[0], len(self.__defs_in_flux) - 1))

    def __flux_box(self, name: str, section_cls: Type[T], is_new: bool = False) -> T:
        if name in self.all_defs and is_new:
            self.logger.error(
                'double definition in legacy metainfo',
                def_name=name, def_type='section')

        definition: Definition = None
        definition = self.all_defs.get(name)
        if definition is None:
            definition = self.__defs_in_flux.get(name)
        if definition is None:
            definition = section_cls(
                name=name,
                description='<this definition was referenced, but is not (yet) defined>')
            self.__defs_in_flux[name] = definition

        if is_new:
            if name in self.__defs_in_flux:
                del(self.__defs_in_flux[name])

        return cast(T, definition)

    def convert_package(
            self, legacy_definitions: List[InfoKindEl], **kwargs) -> Package:
        """ Converts a single legacy metainfo package, i.e. a list of :class:`InfoKindEl`
        into a metainfo package.
        """
        package = Package(**kwargs)

        definition: Definition = None
        for legacy_def in legacy_definitions:
            if legacy_def.kindStr == 'type_abstract_document_content':
                definition = self.__flux_box(legacy_def.name, Category, is_new=True)
                package.m_add_sub_section(Package.category_definitions, definition)

            elif legacy_def.kindStr == 'type_section':
                definition = self.__flux_box(legacy_def.name, Section, is_new=True)
                package.m_add_sub_section(Package.section_definitions, definition)

            elif legacy_def.kindStr in ['type_dimension', 'type_document_content']:
                definition = Quantity(name=legacy_def.name)
                referenced_sections = legacy_def.extra_args.get('referencedSections')
                if referenced_sections is not None and len(referenced_sections) > 0:
                    if len(referenced_sections) == 1:
                        definition.type = Reference(self.__flux_box(referenced_sections[0], Section))

                    else:
                        self.logger.error(
                            'Could not map non higher dimensional reference quantity %s.' %
                            definition.name)
                        definition.type = np.dtype(int)

                elif legacy_def.kindStr == 'type_dimension':
                    definition.type = int
                elif legacy_def.dtypeStr == 'D':
                    definition.type = Any
                elif legacy_def.dtypeStr == 'C':
                    definition.type = str
                elif legacy_def.dtypeStr == 'r':
                    definition.type = int
                elif legacy_def.dtypeStr == 'b':
                    definition.type = bool
                elif legacy_def.dtypeStr == 'i64':
                    definition.type = np.dtype(np.int64)
                else:
                    definition.type = np.dtype(legacy_def.dtypeStr)

                legacy_shape = legacy_def.shape
                if legacy_shape is None:
                    legacy_shape = []

                definition.shape = legacy_shape

                if legacy_def.units is not None:
                    definition.unit = units.parse_units(legacy_def.units)

            else:
                self.logger.error(
                    'unknown kindStr', def_name=legacy_def.name, kind_str=legacy_def.kindStr)
                continue

            self.all_defs[definition.name] = definition

            # common fields
            definition.description = legacy_def.description

            # superNames
            for legacy_super_name in legacy_def.superNames:
                legacy_super_def = self.all_legacy_defs.get(legacy_super_name)
                if legacy_super_def is None:
                    self.logger.error(
                        'super name does not exist', def_name=legacy_def.name,
                        super_name=legacy_super_name)

                if legacy_super_def.kindStr == 'type_section':
                    parent_def = self.__flux_box(legacy_super_name, Section)
                    if isinstance(definition, Section):
                        sub_section = parent_def.m_create(
                            SubSection, name=legacy_def.name, sub_section=definition)
                        sub_section.repeats = legacy_def.repeats is None or legacy_def.repeats

                    elif isinstance(definition, Quantity):
                        parent_def.m_add_sub_section(Section.quantities, definition)

                elif legacy_super_def.kindStr == 'type_abstract_document_content':
                    category = self.__flux_box(legacy_super_name, Category)
                    definition.categories += [category]

                else:
                    self.logger.error(
                        'super name is neither section nor category',
                        def_name=legacy_def.name, super_name=legacy_super_name)

        errors = package.m_all_validate()
        if len(errors) > 0:
            self.logger.error(
                '%s. There are %d more errors in converted legacy package %s' %
                (errors[0], len(errors) - 1, package))
        self.env.m_add_sub_section(Environment.packages, package)

        for section_def in package.section_definitions:  # pylint: disable=not-an-iterable
            attrs = dict(**section_def.all_properties)
            attrs.update(m_def=section_def, do_init=False)
            section_def.section_cls = type(section_def.name, (MSection,), attrs)
        return package

    def legacy_info(self, definition: Definition, *args, **kwargs) -> InfoKindEl:
        """ Creates a legacy metainfo objects for the given definition. """
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
                for sub_section in definition.parent_section_sub_section_defs)

            for sub_section in definition.parent_section_sub_section_defs:
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
        """ Creates a legacy metainfo environment with all definitions from the given packages. """
        if packages is None:
            packages = self.env.packages

        env = InfoKindEnv(*args, **kwargs)
        for package in packages:
            for definition in package.all_definitions.values():
                env.addInfoKindEl(self.legacy_info(definition))
                if isinstance(definition, Section):
                    for quantity in definition.quantities:
                        env.addInfoKindEl(self.legacy_info(quantity))

        return env

    def generate_metainfo_code(
            self, package: Package, directory: str = None, package_name: str = None):

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
                return "MProxy('%s')" % mi_type.target_section_def.name
            else:
                return str(mi_type)

        def format_unit(unit):
            return "'%s'" % unit

        env = JinjaEnvironment(
            loader=PackageLoader('nomad.metainfo', 'templates'),
            autoescape=select_autoescape(['python']))
        env.globals.update(
            format_description=format_description,
            format_type=format_type,
            format_unit=format_unit)

        with open(os.path.join(
                directory, '%s.py' % package_name
                if package_name is not None else package.name), 'wt') as f:
            code = env.get_template('package.j2').render(pkg=package)
            code = '\n'.join([
                line.rstrip() if line.strip() != '' else ''
                for line in code.split('\n')])
            f.write(code)


if __name__ == '__main__':
    """ Converts the old metainfo and code-generates definitions for the new metainfo """
    env = LegacyMetainfoEnvironment(
        metainfo='vasp.nomadmetainfo.json',
        package_names=['%s.nomadmetainfo.json' % pkg for pkg in ['common', 'public', 'vasp']])

    legacy_env = env.legacy_info_env()
    env.generate_metainfo_code(env.env.all_packages['public.nomadmetainfo.json'], package_name='public')
