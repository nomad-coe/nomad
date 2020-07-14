# pylint: skip-file
# type: ignore
from nomad import metainfo
from nomad.datamodel.metainfo import public, common

# Access the quantities of a section definition
for quantity in public.section_method.m_def.quantities:
    print(quantity.name, quantity.type, quantity.shape, quantity.unit)


# Access all the quantities of a section definition, including those added by other packages
import vaspparser.metainfo.vasp  # noqa
for quantity in public.section_method.m_def.all_quantities.values():
    print(quantity.name, quantity.type, quantity.shape, quantity.unit)


# Access sub-sections and their definitions
for sub_section in public.section_run.m_def.sub_sections:
    print(sub_section.name)  # access the name of the sub section definition
    print(sub_section.sub_section.name)  # access the name of the section definition that the sub section refers to


# Go through all sub sections recursively
def visit_section(section, indent=0):
    print(' ' * indent + section.name)
    for sub_section in section.all_sub_sections.values():
        visit_section(sub_section.sub_section, indent + 2)


visit_section(public.section_run.m_def)


# Look at the EntryArchive, e.g. where section_metadata (and everything else) is a subsection
from nomad.datamodel import EntryArchive  # noqa
visit_section(EntryArchive.m_def)


# To get everything within a metainfo package (i.e. what was former in a .nomadmetainfo.json file) as JSON/dict data:
import json  # noqa
import nomad.datamodel.datamodel  # noqa

print(json.dumps(nomad.datamodel.datamodel.m_package.m_to_dict(), indent=2))
print(json.dumps(public.m_package.m_to_dict(), indent=2))


# Using an environment that manages multiple packages and provides utility functions
# to find definitions by name.
from nomad.datamodel.metainfo import m_env  # noqa, contains all common, public, general metainfo
from vaspparser.metainfo import m_env as vasp_m_env  # noqa, contains also the vasp specific definitions
print(m_env.packages)
# Resolve definition by name
print(m_env.resolve_definitions('number_of_atoms', metainfo.Quantity))
# Traverse all definitions:
for definition in m_env.m_all_contents():
    print(definition)


# Dimensions are either numbers or rangens (e.g. 3, 1..3, 0..*) or references to
# shapeless, unitless, integer quantities (usually) of the same section.
# These quantities are not specifically designated as dimensions, because they represent
# quantities in their own right and are often used on their own.
# Dimensions of a specific quantity:
quantity = public.section_system.atom_labels
for dim in quantity.shape:
    if isinstance(dim, str):
        section = quantity.m_parent
        print('%s[%s]: %s' % (quantity.name, dim, m_env.resolve_definitions(dim, metainfo.Quantity)))

# All quantities used as dimensions in a package:
for definition in public.m_package.m_all_contents():
    if definition.m_def == metainfo.Quantity.m_def:
        for dim in definition.shape:
            if isinstance(dim, str) and '..' not in dim:
                print('%s[%s]: %s' % (quantity.name, dim, m_env.resolve_definitions(dim, metainfo.Quantity)))


# Categories are special classes, similar to sections and they Python definition is a
# subclass of MCategory or MSection:
print(public.atom_forces_type, issubclass(public.atom_forces_type, metainfo.MCategory))
print(public.section_system, issubclass(public.section_system, metainfo.MSection))
# Or the definition of the definition is Category or Section respectively:
print(public.atom_forces_type, public.atom_forces_type.m_def == metainfo.Category.m_def)
print(public.section_system, public.section_system.m_def == metainfo.Section.m_def)
# Get all sections and categories definitions in a package:
print(public.m_package.category_definitions)
print(public.m_package.section_definitions)
# Access the categories of a metainfo definition, e.g. quantity
print(public.section_single_configuration_calculation.energy_total.categories)


print(m_env.resolve_definition('EntryMetadata', metainfo.Section).all_quantities)
print(m_env.resolve_definition('Bulk', metainfo.Section).all_quantities)
print(m_env.resolve_definition('OptimadeEntry', metainfo.Section).all_quantities)