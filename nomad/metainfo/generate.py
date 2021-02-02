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
This module contains functionality to generate metainfo source-code from metainfo
definitions.
'''

import numpy as np

from nomad import utils
from nomad.metainfo import Definition, SubSection, Package, Quantity, Section, Reference, MEnum

logger = utils.get_logger(__name__)


def generate_metainfo_code(metainfo_pkg: Package, python_package_path: str):
    '''
    Generates python code with metainfo definitions from the given :class:`Package`.

    Arguments:
        metainfo_pkg: The metainfo package.
        python_package_path:
            The file path for the python module file that should be generated.
    '''
    from jinja2 import Environment as JinjaEnvironment, PackageLoader, select_autoescape
    import textwrap

    def format_description(description, indent=0, width=90):
        paragraphs = [paragraph.strip() for paragraph in description.split('\n\n')]

        def format_paragraph(paragraph, first):
            lines = textwrap.wrap(text=paragraph, width=width - indent * 4)
            lines = [line.replace('\\', '\\\\') for line in lines]
            return textwrap.indent(
                '\n'.join(lines), ' ' * 4 * indent, lambda x: not (first and x.startswith(lines[0])))

        return '\n\n'.join([
            format_paragraph(p, i == 0)
            for i, p in enumerate(paragraphs) if p != ''])

    def format_type(pkg, mi_type):
        if isinstance(mi_type, np.dtype):
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

        if isinstance(mi_type, MEnum):
            return 'MEnum(%s)' % ', '.join(["'%s'" % item for item in mi_type])

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

    def format_aliases(pkg):
        result = ''
        for definition in pkg.category_definitions + pkg.section_definitions:
            for alias in definition.aliases:
                result += '%s = %s\n' % (alias, definition.name)

        if result != '':
            return '\n\n\n%s' % result

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
        fromat_package_import=fromat_package_import,
        format_aliases=format_aliases)

    with open(python_package_path, 'wt') as f:
        code = env.get_template('package_new.j2').render(pkg=metainfo_pkg)
        code = '\n'.join([
            line.rstrip() if line.strip() != '' else ''
            for line in code.split('\n')])
        f.write(code)


if __name__ == '__main__':
    # Simple use case that merges old common/public defs
    import json

    from nomad.metainfo import Category
    from nomad.datamodel.metainfo.public_old import m_package
    from nomad.datamodel.metainfo.common_old import m_package as common_pkg

    for section in common_pkg.section_definitions:  # pylint: disable=not-an-iterable
        if section.extends_base_section:
            base_section = section.base_sections[0]
            for name, attr in section.section_cls.__dict__.items():
                if isinstance(attr, Quantity):
                    base_section.m_add_sub_section(Section.quantities, attr.m_copy(deep=True))
                elif isinstance(attr, SubSection):
                    base_section.m_add_sub_section(Section.sub_sections, attr.m_copy(deep=True))
        else:
            m_package.m_add_sub_section(Package.section_definitions, section)

    for category in common_pkg.category_definitions:  # pylint: disable=not-an-iterable
        m_package.m_add_sub_section(Package.category_definitions, category)

    for definition in m_package.section_definitions + m_package.category_definitions:
        old_name = definition.name
        new_name = ''.join([item[0].title() + item[1:] for item in old_name.split('_')])
        if new_name.startswith('Section'):
            new_name = new_name[7:]
        if new_name != old_name:
            definition.aliases = [old_name]
            definition.name = new_name

    unused_category = m_package.m_create(Category)
    unused_category.name = 'Unused'
    unused_category.description = 'This metainfo definition is not used by NOMAD data.'

    with open('local/metainfostats.json', 'rt') as f:
        stats = json.load(f)

    unused = []
    for (definition, _, _) in m_package.m_traverse():
        if isinstance(definition, (SubSection, Quantity)):
            if definition.name not in stats:
                unused.append(definition)

    for definition in unused:
        if unused_category not in definition.categories:
            definition.categories += [unused_category]

    generate_metainfo_code(m_package, 'nomad/datamodel/metainfo/common_dft.py')
