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


@dev.command(help='Runs tests and linting of the nomad source code. Useful before committing code.')
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
            'skeleton.nomadmetainfo.json',
            'turbomole.nomadmetainfo.json',
            'vasp.nomadmetainfo.json',
            'wien2k.nomadmetainfo.json',
            'dft.nomadmetainfo.json',
            'ems.nomadmetainfo.json']

    for element in path:
        env = convert(element)
        generate_metainfo_code(env)
