# pylint: skip-file
# type: ignore
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
