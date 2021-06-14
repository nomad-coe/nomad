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
    ret_code += os.system('python -m pylint --load-plugins=pylint_mongoengine,nomad.metainfo.pylint_plugin nomad tests')
    click.echo('Run static type checks ...')
    ret_code += os.system('python -m mypy --ignore-missing-imports --follow-imports=silent --no-strict-optional nomad tests')

    sys.exit(ret_code)


@dev.command(help='Runs tests and linting of the nomad gui source code. Useful before committing code.')
@click.option('--skip-tests', help='Do no tests, just do code checks.', is_flag=True)
def gui_qa(skip_tests: bool):
    os.chdir(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../gui')))
    ret_code = 0

    if not skip_tests:
        click.echo('Run gui testing ...')
        ret_code += os.system('yarn run test')
    click.echo('Run gui code linting ...')
    ret_code += os.system('yarn run lint')

    sys.exit(ret_code)


@dev.command(help='Generates a JSON with all metainfo.')
def metainfo():
    import json
    print(json.dumps(metainfo_undecorated(), indent=2))


def metainfo_undecorated():
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

    # TODO similar to before, due to lazyloading, we need to explicily access parsers
    # to actually import all parsers and indirectly all metainfo packages
    from nomad.parsing import parsers
    parsers.parsers

    export = Environment()
    for package in Package.registry.values():
        export.m_add_sub_section(Environment.packages, package)

    return export.m_to_dict(with_meta=True)


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


@dev.command(help='Generates a JSON file that compiles all the parser metadata from each parser project.')
def parser_metadata():
    import inspect
    import re
    import json
    import yaml

    from nomad.parsing import Parser
    from nomad.parsing.parsers import parser_dict

    parsers_metadata = {}
    for parser in parser_dict.values():
        if isinstance(parser, Parser):
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


@dev.command(help=(
    'Updates parser`s README files by combining a general template with  '
    'a parser`s metadata YAML file.'))
@click.option('--parser', help='Only updated the README of the given parsers subdirctory.')
def update_parser_readmes(parser):
    from glob import glob
    import re
    import yaml

    os.chdir(os.path.join(os.path.dirname(__file__), '../..'))

    # filenames
    local_fn = 'README.md'
    generic_fn = './README.parsers.md'
    parser_path = './dependencies/parsers/'

    for num, ddir in enumerate(sorted(glob(parser_path + '*/')), start=1):
        if parser is not None and parser != ddir.split(os.sep)[-2]:
            print(f'Skip {ddir}')
            continue

        _, parser_dirname = os.path.split(ddir)
        print('{} Working on {}' .format(num, parser_dirname))

        # Open general template
        with open(generic_fn, 'r') as generic:  # read only
            body = generic.read()

        # Open local YAML metadata file
        local_mdata = ddir + 'metadata.yaml'
        with open(local_mdata, 'r') as mdata_f:
            mdata = yaml.load(mdata_f, Loader=yaml.FullLoader)
            if mdata.get('codeName', '').strip() == '':
                mdata['codeName'] = os.path.basename(os.path.dirname(ddir))
            if 'preamble' not in mdata:
                mdata['preamble'] = ''

        # Find & Replace Parser`s metadata on its README file
        local_readme = ddir + local_fn
        with open(local_readme, 'w') as local:
            for key in mdata.keys():
                ignore = ['codeLabelStyle', 'parserDirName']
                if key in ignore:
                    continue

                find = r'\$' + key + r'\$'
                replace = mdata[key]
                print(f'\tReplacing {key} with {replace}')

                if key == 'parserSpecific':
                    if mdata[key] != '':
                        replace = r'## Parser Specific\n' + replace

                body = re.sub(find, replace, body)

            # Extra: to replace the codeName (there's no YAML key for this)
            key = 'codeName'
            print('\tReplacing', key)
            find = r'\$' + key + r'\$'
            replace = parser_dirname
            body = re.sub(find, replace, body)

            # Extra: to replace the comment at the top of the gereral template
            find = (
                r'\*\*\*Note:\*\* This is a general README file for NOMAD parsers, '
                r'consult the README of specific parser projects for more detailed '
                r'information!\*\n\n')
            replace = ''
            print('\tReplacing the top comment')
            # print('\t', find, ' -> ', replace)
            body = re.sub(find, replace, body)

            # save file
            local.seek(0)  # go to the top
            local.write(body.strip())
            local.truncate()


@dev.command(help='Creates a Javascript source file containing the required unit conversion factors.')
@click.pass_context
def units(ctx):
    import re
    import json
    from nomad.units import ureg

    # Mapping from unit name to dimension
    unit_map = {
        # Time
        "second": {
            "dimension": "time",
            "label": "Second",
            "abbreviation": "s",
        },
        "atomic_unit_of_time": {
            "dimension": "time",
            "label": "Atomic unit of time",
            "abbreviation": "atomic_unit_of_time",
        },
        # Length
        "meter": {
            "dimension": "length",
            "label": "Meter",
            "abbreviation": "m",
        },
        "bohr": {
            "dimension": "length",
            "label": "Bohr",
            "abbreviation": "bohr",
        },
        "angstrom": {
            "dimension": "length",
            "label": "Ångstrom",
            "abbreviation": "Å",
        },
        # Mass
        "kilogram": {
            "dimension": "mass",
            "label": "Kilogram",
            "abbreviation": "kg",
        },
        "electron_mass": {
            "dimension": "mass",
            "label": "Electron mass",
            "abbreviation": "mₑ",
        },
        "unified_atomic_mass_unit": {
            "dimension": "mass",
            "label": "Unified atomic mass unit",
            "abbreviation": "u",
        },
        # Current
        "ampere": {
            "dimension": "current",
            "label": "Ampere",
            "abbreviation": "A",
        },
        "atomic_unit_of_current": {
            "dimension": "current",
            "label": "Atomic unit of current",
            "abbreviation": "atomic_unit_of_current",
        },
        # Substance
        "mole": {
            "dimension": "substance",
            "label": "Mole",
            "abbreviation": "mole",
        },
        # Luminosity
        "candela": {
            "dimension": "luminosity",
            "label": "Candela",
            "abbreviation": "cd",
        },
        # Temperature
        "kelvin": {
            "dimension": "temperature",
            "label": "Kelvin",
            "abbreviation": "K",
        },
        "celsius": {
            "dimension": "temperature",
            "label": "Celsius",
            "abbreviation": "°C",
        },
        "fahrenheit": {
            "dimension": "temperature",
            "label": "Fahrenheit",
            "abbreviation": "°F",
        },
        "atomic_unit_of_temperature": {
            "dimension": "temperature",
            "label": "Atomic unit of temperature",
            "abbreviation": "atomic_unit_of_temperature",
        },
        # Force
        "newton": {
            "dimension": "force",
            "label": "Newton",
            "abbreviation": "N",
        },
        "atomic_unit_of_force": {
            "dimension": "force",
            "label": "Atomic unit of force",
            "abbreviation": "atomic_unit_of_force",
        },
        # Pressure
        "pascal": {
            "dimension": "pressure",
            "label": "Pascal",
            "abbreviation": "Pa"
        },
        # Energy
        "joule": {
            "dimension": "energy",
            "label": "Joule",
            "abbreviation": "J",
        },
        "electron_volt": {
            "dimension": "energy",
            "label": "Electron volt",
            "abbreviation": "eV",
        },
        "hartree": {
            "dimension": "energy",
            "label": "Hartree",
            "abbreviation": "Ha",
        },
        # Power
        "watt": {
            "dimension": "power",
            "label": "Watt",
            "abbreviation": "W",
        },
        # Frequency
        "hertz": {
            "dimension": "frequency",
            "label": "Hertz",
            "abbreviation": "Hz",
        },
        # Electric potential
        "volt": {
            "dimension": "electric_potential",
            "label": "Volt",
            "abbreviation": "V",
        },
        # Capacitance
        "farad": {
            "dimension": "capacitance",
            "label": "Farad",
            "abbreviation": "F",
        },
        # Charge
        "coulomb": {
            "dimension": "charge",
            "label": "Coulomb",
            "abbreviation": "C",
        },
        "elementary_charge": {
            "dimension": "charge",
            "label": "Elementary charge",
            "abbreviation": "e",
        },
        # Magnetic field
        "tesla": {
            "dimension": "magnetic_field",
            "label": "Tesla",
            "abbreviation": "T",
        },
        # Magnetic flux
        "weber": {
            "dimension": "magnetic_flux",
            "label": "Weber",
            "abbreviation": "Wb",
        },
        # Inductance
        "henry": {
            "dimension": "inductance",
            "label": "Henry",
            "abbreviation": "H",
        },
        # dimensionless
        "dimensionless": {
            "dimension": "dimensionless",
            "label": "Dimensionless",
            "abbreviation": "",
        },
    }

    # Units that are supported
    unit_table = {
        # Base units
        "time": {
            "units": [
                "second",
                "atomic_unit_of_time",
            ],
            "multipliers": {},
        },
        "length": {
            "units": [
                "meter",
                "bohr",
                "angstrom",
            ],
            "multipliers": {},
        },
        "mass": {
            "units": [
                "kilogram",
                "electron_mass",
                "unified_atomic_mass_unit",
            ],
            "multipliers": {},
        },
        "current": {
            "units": [
                "ampere",
                "atomic_unit_of_current",
            ],
            "multipliers": {},
        },
        "substance": {
            "units": [
                "mole",
            ],
            "multipliers": {},
        },
        "luminosity": {
            "units": [
                "candela",
            ],
            "multipliers": {},
        },
        "temperature": {
            "units": [
                "kelvin",
                "celsius",
                "fahrenheit",
                "atomic_unit_of_temperature",
            ],
            "multipliers": {},
            "constants": {},
        },
        # Derived units
        "force": {
            "units": [
                "newton",
                "atomic_unit_of_force",
            ],
            "multipliers": {},
        },
        "pressure": {
            "units": [
                "pascal",
            ],
            "multipliers": {},
        },
        "energy": {
            "units": [
                "joule",
                "electron_volt",
                "hartree",
            ],
            "multipliers": {},
        },
        "power": {
            "units": [
                "watt",
            ],
            "multipliers": {},
        },
        "frequency": {
            "units": [
                "hertz",
            ],
            "multipliers": {},
        },
        "electric_potential": {
            "units": [
                "volt",
            ],
            "multipliers": {},
        },
        "capacitance": {
            "units": [
                "farad",
            ],
            "multipliers": {},
        },
        "charge": {
            "units": [
                "coulomb",
                "elementary_charge",
            ],
            "multipliers": {},
        },
        "magnetic_field": {
            "units": [
                "tesla",
            ],
            "multipliers": {},
        },
        "magnetic_flux": {
            "units": [
                "weber",
            ],
            "multipliers": {},
        },
        "inductance": {
            "dimension": "inductance",
            "units": [
                "henry",
            ],
            "multipliers": {},
        },
        "dimensionless": {
            "dimension": "dimensionless",
            "units": [
                "dimensionless",
            ],
            "multipliers": {},
        },
    }

    # Unit systems
    unit_systems = {
        "SI": {
            "label": "SI",
            "description": "International System of Units (SI)",
            "units": {
                "time": "second",
                "length": "meter",
                "mass": "kilogram",
                "current": "ampere",
                "substance": "mole",
                "luminosity": "candela",
                "temperature": "kelvin",
                "force": "newton",
                "pressure": "pascal",
                "energy": "joule",
                "power": "watt",
                "frequency": "hertz",
                "electric_potential": "volt",
                "charge": "coulomb",
            },
        },
        "AU": {
            "label": "Atomic units",
            "description": "Hartree atomic units",
            "units": {
                "time": "atomic_unit_of_time",
                "length": "bohr",
                "mass": "electron_mass",
                "current": "atomic_unit_of_current",
                "temperature": "atomic_unit_of_temperature",
                "force": "atomic_unit_of_force",
                "energy": "hartree",
            }
        }
    }

    # Precompute conversion factors and possible shifts
    for value in unit_table.values():
        units = value["units"]
        for i_unit in units:
            for j_unit in units:
                # Create dictionaries if not present
                multipliers = value["multipliers"]
                if i_unit not in multipliers:
                    multipliers[i_unit] = {}
                if j_unit not in multipliers[i_unit]:
                    multipliers[i_unit][j_unit] = {}

                # Check if there is a constant shift: y = ax + b -> y(0) = b.
                # Uses delta units for temperatures since only they are
                # multiplicative.
                y_0 = 0 * getattr(ureg, "delta_" + i_unit, getattr(ureg, i_unit))
                b = y_0.to(getattr(ureg, j_unit)).magnitude

                # Solving the multiplication factor with:
                # y(1) = a + b -> a = y(1) - b
                # Uses delta units for temperatures since only they are
                # multiplicative. Causes minor numerical accuracy issues in the
                # factor due to floating point precision
                y_1 = 1 * getattr(ureg, "delta_" + i_unit, getattr(ureg, i_unit))
                a = y_1.to(getattr(ureg, "delta_" + j_unit, getattr(ureg, j_unit))).magnitude
                multipliers[i_unit][j_unit] = a

                # Create dictionaries if not present
                if b != 0:
                    constants = value["constants"]
                    if i_unit not in constants:
                        constants[i_unit] = {}
                    if j_unit not in constants[i_unit]:
                        constants[i_unit][j_unit] = {}
                    constants[i_unit][j_unit] = b

    # Check that all defined units are present in the conversion table
    for unit, info in unit_map.items():
        dimension = info["dimension"]
        unit_list = unit_table[dimension]["units"]
        assert unit in unit_list, "Could not find '{}' in the unit table under the dimension '{}'.".format(unit, dimension)

    # Check unit system correctness
    for system in unit_systems.values():
        for dimension, unit in system["units"].items():
            info = unit_map[unit]
            assert dimension == info["dimension"]

    # Go through the metainfo and check that all units are defined
    all_metainfo = metainfo_undecorated()
    units = set()
    packages = all_metainfo["packages"]
    for package in packages:
        sections = package["section_definitions"]
        for section in sections:
            quantities = section.get("quantities", [])
            for quantity in quantities:
                unit = quantity.get("unit")
                if unit is not None:
                    parts = unit.split()
                    for part in parts:
                        is_operator = part in {"/", "**", "*"}
                        is_number = True
                        try:
                            int(part)
                        except Exception:
                            is_number = False
                        if not is_operator and not is_number:
                            units.add(part)
    for unit in units:
        assert unit in unit_map, "The unit '{}' is not defined in the unit definitions.".format(unit)

    # Print unit conversion table and unit systems as a Javascript source file
    output = """/*
 * Copyright The NOMAD Authors.
 *
 * This file is part of NOMAD. See https://nomad-lab.eu for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
// Generated by NOMAD CLI. Do not edit manually.
"""
    output += "export const unitMap = "
    json_string = json.dumps(unit_map, indent=2)
    json_string = re.sub(r'(?<!: )"(\S*?)":', '\\1:', json_string)
    json_string = json_string.replace("\"", "'")
    output += json_string
    output += "\nexport const conversionMap = "
    json_string = json.dumps(unit_table, indent=2)
    json_string = re.sub(r'(?<!: )"(\S*?)":', '\\1:', json_string)
    json_string = json_string.replace("\"", "'")
    output += json_string
    output += "\nexport const unitSystems = "
    json_string = json.dumps(unit_systems, indent=2)
    json_string = re.sub(r'(?<!: )"(\S*?)":', '\\1:', json_string)
    json_string = json_string.replace("\"", "'")
    output += json_string
    print(output)
