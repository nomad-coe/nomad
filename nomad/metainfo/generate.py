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
from nomad.metainfo import Definition, Package, Reference, MEnum

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
        code = env.get_template('package.j2').render(pkg=metainfo_pkg)
        code = '\n'.join([
            line.rstrip() if line.strip() != '' else ''
            for line in code.split('\n')])
        f.write(code)


if __name__ == '__main__':
    # Simple use case that re-generates the common_dft package
    from nomad.datamodel.metainfo.common_dft import m_package

    generate_metainfo_code(m_package, 'nomad/datamodel/metainfo/common_dft.py')
