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

"""
The *parsing* module is an interface for the existing NOMAD-coe parsers.
The parser code is used via :mod:`nomad.dependencies`. This module redefines
some of the old NOMAD-coe python-common functionality to create a more coherent
interface to the parsers.

Assumption about parsers
------------------------
For now, we make a few assumption about parsers
- they always work on the same *meta-info* version
- they have no conflicting python requirments
- they can be loaded at the same time and can be used within the same python process
- they are uniquely identified by a GIT URL and publicly accessible
- their version is uniquly identified by a GIT commit SHA

Each parser is defined via an instance of :class:`Parser`.

.. autoclass:: nomad.parsing.Parser
    :members:

The implementation :class:`LegacyParser` is used for most NOMAD-coe parsers.

.. autoclass:: nomad.parsing.LegacyParser

The parser definitions are available via the following two variables.

.. autodata:: nomad.parsing.parsers
.. autodata:: nomad.parsing.parser_dict

Parsers are reused for multiple caclulations.

Parsers in NOMAD-coe use a *backend* to create output. There are different NOMAD-coe
basends. In nomad@FAIR, we only currently only use a single backed. A version of
NOMAD-coe's *LocalBackend*. It stores all parser results in memory. The following
classes provide a interface definition for *backends* as an ABC and a concrete implementation
based on NOMAD-coe's *python-common* module.

.. autoclass:: nomad.parsing.AbstractParserBackend
    :members:
.. autoclass:: nomad.parsing.LocalBackend
    :members:

"""

from nomad.parsing.backend import AbstractParserBackend, LocalBackend, LegacyLocalBackend, JSONStreamWriter, BadContextURI, WrongContextState
from nomad.parsing.parser import Parser, LegacyParser
from nomad.parsing.artificial import TemplateParser, GenerateRandomParser
from nomad.dependencies import dependencies_dict as dependencies

parsers = [
    GenerateRandomParser(),
    TemplateParser(),
    LegacyParser(
        python_git=dependencies['parsers/vasp'],
        parser_class_name='vaspparser.VASPRunParserInterface',
        main_file_re=r'^.*\.xml(\.[^\.]*)?$',
        main_contents_re=(
            r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
            r'?\s*<modeling>'
            r'?\s*<generator>'
            r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
            r'?')
    ),
    LegacyParser(
        python_git=dependencies['parsers/exciting'],
        parser_class_name='excitingparser.ExcitingParser',
        main_file_re=r'^.*/INFO\.OUT?',
        main_contents_re=(
            r'^\s*=================================================+\s*'
            r'\s*\|\s*EXCITING\s+\S+\s+started\s*='
            r'\s*\|\s*version hash id:\s*\S*\s*=')
    ),
    LegacyParser(
        python_git=dependencies['parsers/fhi-aims'],
        parser_class_name='fhiaimsparser.FHIaimsParser',
        main_file_re=r'^.*\.out$',
        main_contents_re=(
            r'^(.*\n)*'
            r'?\s*Invoking FHI-aims \.\.\.'
            r'?\s*Version')
    ),
    LegacyParser(
        python_git=dependencies['parsers/cp2k'],
        parser_class_name='cp2kparser.CP2KParser',
        main_file_re=r'^.*\.out$',  # This looks for files with .out
        main_contents_re=(
            r'\*\*\*\* \*\*\*\* \*\*\*\*\*\*  \*\*  PROGRAM STARTED AT\s.*\n'
            r' \*\*\*\*\* \*\* \*\*\*  \*\*\* \*\*   PROGRAM STARTED ON\s*.*\n'
            r' \*\*    \*\*\*\*   \*\*\*\*\*\*    PROGRAM STARTED BY .*\n'
            r' \*\*\*\*\* \*\*    \*\* \*\* \*\*   PROGRAM PROCESS ID .*\n'
            r'  \*\*\*\* \*\*  \*\*\*\*\*\*\*  \*\*  PROGRAM STARTED IN .*\n')
    ),
]
""" Instanciation and constructor based config of all parsers. """

parser_dict = {parser.name: parser for parser in parsers}  # type: ignore
""" A dict to access parsers by name. Usually 'parsers/<...>', e.g. 'parsers/vasp'. """
