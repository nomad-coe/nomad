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

import sys
import os
import click

from nomad.cli.cli import cli


@cli.group(help='Commands related to the nomad source code.')
def dev():
    pass


@dev.command(help='Runs tests and linting of the nomad python source code. Useful before committing code.')
@click.option('--skip-tests', help='Do no tests, just do code checks.', is_flag=True)
@click.option('-x', '--exitfirst', help='Stop testing after first failed test case.', is_flag=True)
def qa(skip_tests: bool, exitfirst: bool):
    os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
    ret_code = 0
    if not skip_tests:
        click.echo('Run tests ...')
        ret_code += os.system('python -m pytest -sv%s tests' % ('x' if exitfirst else ''))
    click.echo('Run code style checks ...')
    ret_code += os.system('python -m pycodestyle --ignore=E501,E701,E731 nomad tests')
    click.echo('Run linter ...')
    ret_code += os.system('python -m pylint --load-plugins=pylint_mongoengine,nomad/metainfo/pylint_plugin nomad tests')
    click.echo('Run static type checks ...')
    ret_code += os.system('python -m mypy --ignore-missing-imports --follow-imports=silent --no-strict-optional nomad tests')

    sys.exit(ret_code)


@dev.command(help='Runs tests and linting of the nomad gui source code. Useful before committing code.')
def gui_qa():
    click.echo('Run gui code linting ...')
    os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../gui')))
    ret_code = 0
    ret_code += os.system('yarn run eslint \'src/**/*.js\'')
    sys.exit(ret_code)


@dev.command(help='Generates a JSON with all metainfo.')
def metainfo():
    import json
    from nomad.metainfo import Package, Environment

    # TODO the __init_metainfo__() should not be necessary and automatically performed
    # Also load and initialize the datamodel definitions
    import nomad.metainfo.metainfo
    import nomad.datamodel.datamodel
    import nomad.datamodel.dft
    import nomad.datamodel.ems
    import nomad.datamodel.optimade
    import nomad.datamodel.encyclopedia
    nomad.metainfo.metainfo.m_package.__init_metainfo__()
    nomad.datamodel.datamodel.m_package.__init_metainfo__()
    nomad.datamodel.dft.m_package.__init_metainfo__()  # pylint: disable=no-member
    nomad.datamodel.ems.m_package.__init_metainfo__()  # pylint: disable=no-member
    nomad.datamodel.optimade.m_package.__init_metainfo__()  # pylint: disable=no-member
    nomad.datamodel.encyclopedia.m_package.__init_metainfo__()

    # Ensure all parser metainfo is loaded
    from nomad.parsing.parsers import parsers
    for parser in parsers:
        _ = parser.metainfo_env

    export = Environment()
    for package in Package.registry.values():
        export.m_add_sub_section(Environment.packages, package)

    print(json.dumps(export.m_to_dict(with_meta=True), indent=2))


@dev.command(help='Generates a JSON with all search quantities.')
def search_quantities():
    from nomad import search
    # Due to this import, the parsing module will register all code_names based on parser
    # implementations.
    from nomad.parsing.parsers import parser_dict  # pylint: disable=unused-import
    import json

    def to_dict(search_quantity):
        result = {
            'name': search_quantity.qualified_name,
            'description': search_quantity.description,
            'many': search_quantity.many,
        }

        if search_quantity.statistic_fixed_size is not None:
            result['statistic_size'] = search_quantity.statistic_fixed_size
        if search_quantity.statistic_values is not None:
            result['statistic_values'] = search_quantity.statistic_values

        return result

    export = {
        search_quantity.qualified_name: to_dict(search_quantity)
        for search_quantity in search.search_quantities.values()
    }
    print(json.dumps(export, indent=2))


@dev.command(help='Generates source-code for the new metainfo from .json files of the old.')
@click.argument('path', nargs=-1)
def legacy_metainfo(path):
    from nomad.metainfo.legacy import convert, generate_metainfo_code

    if len(path) == 0:
        path = [
            'abinit.nomadmetainfo.json',
            'aptfim.nomadmetainfo.json',
            'atk.nomadmetainfo.json',
            'band.nomadmetainfo.json',
            'bigdft.nomadmetainfo.json',
            'castep.nomadmetainfo.json',
            'cp2k.nomadmetainfo.json',
            'cpmd.nomadmetainfo.json',
            'crystal.nomadmetainfo.json',
            'dl_poly.nomadmetainfo.json',
            'dmol3.nomadmetainfo.json',
            'eels.nomadmetainfo.json',
            'elastic.nomadmetainfo.json',
            'elk.nomadmetainfo.json',
            'exciting.nomadmetainfo.json',
            'fhi_aims.nomadmetainfo.json',
            'fleur.nomadmetainfo.json',
            'gamess.nomadmetainfo.json',
            'gaussian.nomadmetainfo.json',
            'gpaw.nomadmetainfo.json',
            'gulp.nomadmetainfo.json',
            'lib_atoms.nomadmetainfo.json',
            'molcas.nomadmetainfo.json',
            'mpes.nomadmetainfo.json',
            'nwchem.nomadmetainfo.json',
            'octopus.nomadmetainfo.json',
            'onetep.nomadmetainfo.json',
            'orca.nomadmetainfo.json',
            'phonopy.nomadmetainfo.json',
            'photoemission.nomadmetainfo.json',
            'qbox.nomadmetainfo.json',
            'quantum_espresso.nomadmetainfo.json',
            'siesta.nomadmetainfo.json',
            'turbomole.nomadmetainfo.json',
            'vasp.nomadmetainfo.json',
            'wien2k.nomadmetainfo.json',
            'dft.nomadmetainfo.json',
            'ems.nomadmetainfo.json']

    for element in path:
        env = convert(element)
        generate_metainfo_code(env)


@dev.command(help='Generates a JSON file that compiles all the parser metadata from each parser project.')
def parser_metadata():
    import inspect
    import re
    import json
    import yaml

    from nomad.parsing import LegacyParser, FairdiParser
    from nomad.parsing.parsers import parser_dict

    parsers_metadata = {}
    for parser in parser_dict.values():
        if isinstance(parser, LegacyParser):
            parser_class = parser.parser_class
        elif isinstance(parser, FairdiParser):
            parser_class = parser.__class__
        else:
            continue
        parser_code_file = inspect.getfile(parser_class)

        path_match = re.match(r'(.*/dependencies/parsers/[^/]+)/.*', parser_code_file)
        if path_match:
            parser_metadata_file = os.path.join(path_match.group(1), 'metadata.yaml')
            if os.path.exists(parser_metadata_file):
                with open(parser_metadata_file) as f:
                    parser_metadata = yaml.load(f, Loader=yaml.FullLoader)
                parsers_metadata[parser.code_name] = parser_metadata

    parsers_metadata = {
        key: parsers_metadata[key]
        for _, key in sorted([(key.lower(), key) for key in parsers_metadata], key=lambda x: x[0])}

    print(json.dumps(parsers_metadata, indent=2))


@dev.command(help='Generate toolkit tutorial metadata from anaytics submodules.')
def toolkit_metadata():
    import requests
    import re
    import json
    modules = requests.get(
        'https://gitlab.mpcdf.mpg.de/api/v4/projects/3161/repository/files/.gitmodules/raw?ref=master').text

    tutorials = []
    lines = modules.split('\n')
    for line in lines:
        match = re.match(r'\s*url = (.*)$', line)
        if match:
            url = match.group(1).replace('.git', '') + '/-/raw/master/metainfo.json'
            response = requests.get(url)
            if response.status_code != 200:
                print('Could not get metadata for %s' % match.group(1), file=sys.stderr)
                continue
            try:
                tutorials.append(response.json())
            except Exception:
                print('Could not get metadata for %s. Project is probably not public.' % match.group(1), file=sys.stderr)

    print(json.dumps(dict(tutorials=tutorials), indent=2))
