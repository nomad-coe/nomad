# This example reads all Python define metainfo and computes
# basic metrics on how many sections and properties are defined
# per package

from nomad.metainfo import Section, Property, Quantity, Package
from nomad.datamodel import all_metainfo_packages

all_sections = 0
all_properties = 0
all_quantities = 0
sections = dict()
properties = dict()
quantities = dict()
definitions = set()
metainfo = all_metainfo_packages()

for definition, _, _, _ in metainfo.m_traverse():
    if definition in definitions:
        continue
    definitions.add(definition)

    package = definition
    while package and not isinstance(package, Package):
        package = package.m_parent
    if package is None:
        continue

    if isinstance(definition, Section):
        all_sections += 1
        sections[package] = sections.get(package, 0) + 1
    if isinstance(definition, Property):
        all_properties += 1
        properties[package] = properties.get(package, 0) + 1
    if isinstance(definition, Quantity):
        all_quantities += 1
        quantities[package] = quantities.get(package, 0) + 1


for package in sections.keys():
    print(
        f'{package.name}: {sections[package]}, {properties.get(package, 0)}, {quantities.get(package, 0)}'
    )

print(f'SUM: {all_sections}, {all_properties}, {all_quantities}')
