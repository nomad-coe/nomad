from typing import Tuple, Dict, List, Type, TypeVar
import os.path

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
import nomad_meta_info

from nomad import utils
from nomad.metainfo import Definition, Package, Category, Section, Quantity, SubSection


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


legacy_defs, legacy_packages = load_legacy_metainfo(['common.nomadmetainfo.json', 'public.nomadmetainfo.json'])

logger = utils.get_logger(__name__)
all_defs: Dict[str, Definition] = dict()

T = TypeVar('T', bound=Definition)
dtype_strs = set()


def convert_package(legacy_definitions: List[InfoKindEl], **kwargs) -> Package:
    package = Package(**kwargs)

    def flux_box(legacy_name: str, section_cls: Type[T], is_new: bool = False) -> T:
        if legacy_def.name in all_defs and is_new:
            logger.error(
                'double definition in legacy metainfo',
                def_name=legacy_def.name, def_type='section')

        definition = package.all_definitions.get(legacy_name)
        if definition is None:
            definition = package.m_create(section_cls, name=legacy_name)

        if is_new:
            all_defs[legacy_def.name] = definition

        return definition

    for legacy_def in legacy_definitions:
        if legacy_def.kindStr == 'type_abstract_document_content':
            definition = flux_box(legacy_def.name, Category, is_new=True)

        elif legacy_def.kindStr == 'type_section':
            definition = flux_box(legacy_def.name, Section, is_new=True)

        elif legacy_def.kindStr in ['type_dimension', 'type_document_content']:
            definition = Quantity(name=legacy_def.name, type=int)
            # map shape, map type
            dtype_strs.add(legacy_def.dtypeStr)

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
                    parent_def.m_create(
                        SubSection, name=legacy_def.name, sub_section=definition)

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

# print(common_pkg.m_to_json(indent=2))
print(dtype_strs)
