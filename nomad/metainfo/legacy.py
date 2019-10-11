from typing import Tuple, Dict, List, Type, TypeVar, Any, cast
import os.path
import numpy as np
import json
from jinja2 import Environment, PackageLoader, select_autoescape
import textwrap

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
import nomad_meta_info

from nomad import utils
from nomad.metainfo import Definition, Package, Category, Section, Quantity, SubSection, Reference, units


def load_legacy_metainfo(
        package_names: List[str] = None) \
        -> Tuple[Dict[str, InfoKindEl], Dict[str, List[InfoKindEl]]]:
    """ Loads the old metainfo and returns them by package, and by kind. """

    meta_info_path = os.path.normpath(os.path.join(
        os.path.dirname(nomad_meta_info.__file__), 'all.nomadmetainfo.json'))

    meta_info_env, _ = loadJsonFile(filePath=meta_info_path)

    defs: Dict[str, InfoKindEl] = {}
    packages: Dict[str, List[InfoKindEl]] = {}
    for definition in meta_info_env.infoKindEls():
        defs[definition.name] = definition
        if package_names is None or definition.package in package_names:
            packages.setdefault(definition.package, []).append(definition)

    return defs, packages


legacy_defs, legacy_packages = load_legacy_metainfo(['common.nomadmetainfo.json', 'public.nomadmetainfo.json', 'vasp.nomadmetainfo.json'])

logger = utils.get_logger(__name__)
all_defs: Dict[str, Definition] = dict()

T = TypeVar('T', bound=Definition)


def convert_package(legacy_definitions: List[InfoKindEl], **kwargs) -> Package:
    package = Package(**kwargs)

    def flux_box(legacy_name: str, section_cls: Type[T], is_new: bool = False) -> T:
        if legacy_def.name in all_defs and is_new:
            logger.error(
                'double definition in legacy metainfo',
                def_name=legacy_def.name, def_type='section')

        definition: Definition = package.all_definitions.get(legacy_name)
        if definition is None:
            definition = package.m_create(
                section_cls, name=legacy_name, description=legacy_def.description)

        if is_new:
            all_defs[legacy_def.name] = definition

        return cast(T, definition)

    definition: Definition = None
    for legacy_def in legacy_definitions:
        if legacy_def.kindStr == 'type_abstract_document_content':
            definition = flux_box(legacy_def.name, Category, is_new=True)

        elif legacy_def.kindStr == 'type_section':
            definition = flux_box(legacy_def.name, Section, is_new=True)

        elif legacy_def.kindStr in ['type_dimension', 'type_document_content']:
            definition = Quantity(
                name=legacy_def.name, description=legacy_def.description)
            referenced_sections = legacy_def.extra_args.get('referencedSections')
            if referenced_sections is not None and len(referenced_sections) > 0:
                if len(referenced_sections) == 1:
                    definition.type = Reference(flux_box(referenced_sections[0], Section))

                else:
                    logger.error('Could not map non higher dimensional reference quantity %s.' % definition.name)
                    definition.type = np.dtype(int)

            elif legacy_def.kindStr == 'type_dimension':
                definition.type = int
            elif legacy_def.dtypeStr == 'D':
                definition.type = Any
            elif legacy_def.dtypeStr == 'C':
                definition.type = str
            elif legacy_def.dtypeStr == 'r':
                definition.type = int
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
            logger.error(
                'unknown kindStr', def_name=legacy_def.name, kind_str=legacy_def.kindStr)

        # superNames
        for legacy_super_name in legacy_def.superNames:
            legacy_super_def = legacy_defs.get(legacy_super_name)
            if legacy_super_def is None:
                logger.error(
                    'super name does not exist', def_name=legacy_def.name,
                    super_name=legacy_super_name)

            if legacy_super_def.kindStr == 'type_section':
                parent_def = flux_box(legacy_super_name, Section)
                if isinstance(definition, Section):
                    sub_section = parent_def.m_create(
                        SubSection, name=legacy_def.name, sub_section=definition)
                    sub_section.repeats = legacy_def.repeats is not None and legacy_def.repeats

                elif isinstance(definition, Quantity):
                    parent_def.m_add_sub_section(Section.quantities, definition)

            elif legacy_super_def.kindStr == 'type_abstract_document_content':
                category = flux_box(legacy_super_name, Category)
                definition.categories += [category]

            else:
                logger.error(
                    'super name is neither section nor category',
                    def_name=legacy_def.name, super_name=legacy_super_name)
    return package


common_pkg = convert_package(
    legacy_packages['common.nomadmetainfo.json'] + legacy_packages['public.nomadmetainfo.json'],
    name='common')

vasp_pkg = convert_package(legacy_packages['vasp.nomadmetainfo.json'], name='vasp')

for error in common_pkg.m_all_validate() + vasp_pkg.m_all_validate():
    print(error)

json.dumps([common_pkg.m_to_dict(), vasp_pkg.m_to_dict()], indent=2)


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


env = Environment(
    loader=PackageLoader('nomad.metainfo', 'templates'),
    autoescape=select_autoescape(['python']))
env.globals.update(
    format_description=format_description,
    format_type=format_type,
    format_unit=format_unit)

with open(os.path.join(os.path.dirname(__file__), 'common.py'), 'wt') as f:
    f.write(env.get_template('package.j2').render(pkg=common_pkg))
