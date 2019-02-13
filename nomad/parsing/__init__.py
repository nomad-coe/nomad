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
This module redefines some of the old NOMAD-coe python-common functionality to create a
more coherent interface to the parsers.

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
basends. In nomad@FAIRDI, we only currently only use a single backed. A version of
NOMAD-coe's *LocalBackend*. It stores all parser results in memory. The following
classes provide a interface definition for *backends* as an ABC and a concrete implementation
based on NOMAD-coe's *python-common* module.

.. autoclass:: nomad.parsing.AbstractParserBackend
    :members:
.. autoclass:: nomad.parsing.LocalBackend
    :members:

"""

from nomad.parsing.backend import AbstractParserBackend, LocalBackend, LegacyLocalBackend, JSONStreamWriter, BadContextURI, WrongContextState
from nomad.parsing.parser import Parser, LegacyParser, VaspOutcarParser
from nomad.parsing.artificial import TemplateParser, GenerateRandomParser, ChaosParser


parsers = [
    GenerateRandomParser(),
    TemplateParser(),
    ChaosParser(),
    LegacyParser(
        name='parsers/vasp',
        parser_class_name='vaspparser.VASPRunParserInterface',
        main_file_re=r'^.*\.xml(\.[^\.]*)?$',
        main_contents_re=(
            r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
            r'?\s*<modeling>'
            r'?\s*<generator>'
            r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
            r'?')
    ),
    VaspOutcarParser(
        name='parsers/vasp',
        parser_class_name='vaspparser.VaspOutcarParser',
        main_file_re=r'^OUTCAR(\.[^\.]*)?$',
        main_contents_re=(r'^\svasp\..*$')
    ),
    LegacyParser(
        name='parsers/exciting',
        parser_class_name='excitingparser.ExcitingParser',
        main_file_re=r'^.*/INFO\.OUT?',
        main_contents_re=(
            r'^\s*=================================================+\s*'
            r'\s*\|\s*EXCITING\s+\S+\s+started\s*='
            r'\s*\|\s*version hash id:\s*\S*\s*=')
    ),
    LegacyParser(
        name='parsers/fhi-aims',
        parser_class_name='fhiaimsparser.FHIaimsParser',
        main_file_re=r'^.*\.out$',
        main_contents_re=(
            r'^(.*\n)*'
            r'?\s*Invoking FHI-aims \.\.\.'
            r'?\s*Version')
    ),
    LegacyParser(
        name='parsers/cp2k',
        parser_class_name='cp2kparser.CP2KParser',
        main_file_re=r'^.*\.out$',  # This looks for files with .out
        main_contents_re=(
            r'\*\*\*\* \*\*\*\* \*\*\*\*\*\*  \*\*  PROGRAM STARTED AT\s.*\n'
            r' \*\*\*\*\* \*\* \*\*\*  \*\*\* \*\*   PROGRAM STARTED ON\s*.*\n'
            r' \*\*    \*\*\*\*   \*\*\*\*\*\*    PROGRAM STARTED BY .*\n'
            r' \*\*\*\*\* \*\*    \*\* \*\* \*\*   PROGRAM PROCESS ID .*\n'
            r'  \*\*\*\* \*\*  \*\*\*\*\*\*\*  \*\*  PROGRAM STARTED IN .*\n'
        )
    ),
    LegacyParser(
        name='parsers/crystal',
        parser_class_name='crystalparser.CrystalParser',
        main_file_re=r'^.*\.out$',
        main_contents_re=(
            r'\s*[\*]{22,}'  # Looks for '*' 22 times or more in a row.
            r'\s*\*\s{20,}\*'  # Looks for a '*' sandwhiched by whitespace.
            r'\s*\*\s{10,}CRYSTAL(?P<majorVersion>[\d]+)\s{10,}\*'
            r'\s*\*\s{10,}public \: (?P<minorVersion>[\d\.]+) \- .*\*'
        )
    ),
    # The main contents regex of CPMD was causing a catostrophic backtracking issue
    # when searching through the first 500 bytes of main files. We decided
    # to use only a portion of the regex to avoid that issue.
    LegacyParser(
        name='parsers/cpmd',
        parser_class_name='cpmdparser.CPMDParser',
        main_file_re=r'^.*\.out$',
        main_contents_re=(
            # r'\s+\*\*\*\*\*\*  \*\*\*\*\*\*    \*\*\*\*  \*\*\*\*  \*\*\*\*\*\*\s*'
            # r'\s+\*\*\*\*\*\*\*  \*\*\*\*\*\*\*   \*\*\*\*\*\*\*\*\*\*  \*\*\*\*\*\*\*\s+'
            r'\*\*\*       \*\*   \*\*\*  \*\* \*\*\*\* \*\*  \*\*   \*\*\*'
            # r'\s+\*\*        \*\*   \*\*\*  \*\*  \*\*  \*\*  \*\*    \*\*\s+'
            # r'\s+\*\*        \*\*\*\*\*\*\*   \*\*      \*\*  \*\*    \*\*\s+'
            # r'\s+\*\*\*       \*\*\*\*\*\*    \*\*      \*\*  \*\*   \*\*\*\s+'
            # r'\s+\*\*\*\*\*\*\*  \*\*        \*\*      \*\*  \*\*\*\*\*\*\*\s+'
            # r'\s+\*\*\*\*\*\*  \*\*        \*\*      \*\*  \*\*\*\*\*\*\s+'
        )
    ),
    LegacyParser(
        name='parsers/nwchem',
        parser_class_name='nwchemparser.NWChemParser',
        main_file_re=r'^.*\.out$',
        main_contents_re=(
            r'\s+Northwest Computational Chemistry Package \(NWChem\) \d+\.\d+'
            r'\s+------------------------------------------------------'
            r'\s+Environmental Molecular Sciences Laboratory'
            r'\s+Pacific Northwest National Laboratory'
            r'\s+Richland, WA 99352'
        )
    ),
    LegacyParser(
        name='parsers/bigdft',
        parser_class_name='bigdftparser.BigDFTParser',
        main_file_re=r'^.*\.out$',
        main_contents_re=(
            r'__________________________________ A fast and precise DFT wavelet code\s*'
            r'\|     \|     \|     \|     \|     \|\s*'
            r'\|     \|     \|     \|     \|     \|      BBBB         i       gggggg\s*'
            r'\|_____\|_____\|_____\|_____\|_____\|     B    B               g\s*'
            r'\|     \|  :  \|  :  \|     \|     \|    B     B        i     g\s*'
            r'\|     \|-0\+--\|-0\+--\|     \|     \|    B    B         i     g        g\s*'
            r'\|_____\|__:__\|__:__\|_____\|_____\|___ BBBBB          i     g         g\s*'
            r'\|  :  \|     \|     \|  :  \|     \|    B    B         i     g         g\s*'
            r'\|--\+0-\|     \|     \|-0\+--\|     \|    B     B     iiii     g         g\s*'
            r'\|__:__\|_____\|_____\|__:__\|_____\|    B     B        i      g        g\s*'
            r'\|     \|  :  \|  :  \|     \|     \|    B BBBB        i        g      g\s*'
            r'\|     \|-0\+--\|-0\+--\|     \|     \|    B        iiiii          gggggg\s*'
            r'\|_____\|__:__\|__:__\|_____\|_____\|__BBBBB\s*'
            r'\|     \|     \|     \|  :  \|     \|                           TTTTTTTTT\s*'
            r'\|     \|     \|     \|--\+0-\|     \|  DDDDDD          FFFFF        T\s*'
            r'\|_____\|_____\|_____\|__:__\|_____\| D      D        F        TTTT T\s*'
            r'\|     \|     \|     \|  :  \|     \|D        D      F        T     T\s*'
            r'\|     \|     \|     \|--\+0-\|     \|D         D     FFFF     T     T\s*'
            r'\|_____\|_____\|_____\|__:__\|_____\|D___      D     F         T    T\s*'
            r'\|     \|     \|  :  \|     \|     \|D         D     F          TTTTT\s*'
            r'\|     \|     \|--\+0-\|     \|     \| D        D     F         T    T\s*'
            r'\|_____\|_____\|__:__\|_____\|_____\|          D     F        T     T\s*'
            r'\|     \|     \|     \|     \|     \|         D               T    T\s*'
            r'\|     \|     \|     \|     \|     \|   DDDDDD       F         TTTT\s*'
            r'\|_____\|_____\|_____\|_____\|_____\|______                    www\.bigdft\.org'
        )
    ),
    LegacyParser(
        name='parsers/wien2k',
        parser_class_name='wien2kparser.Wien2kParser',
        main_file_re=r'^.*\.scf$',  # This looks for files with .scf
        main_contents_re=r':ITE[0-9]+:  1. ITERATION'
    ),
    LegacyParser(
        name='parsers/gaussian',
        parser_class_name='gaussianparser.GaussianParser',
        main_file_re=r'^.*\.out$',
        main_contents_re=(
            r'\s*Cite this work as:'
            r'\s*Gaussian [0-9]+, Revision [A-Za-z0-9.]*,'
            r'\s\*\*\*\*\*\*\*\*\*\*\*\**'
            r'\s*Gaussian\s*(?P<program_version>[0-9]+):\s*(?P<x_gaussian_program_implementation>[A-Za-z0-9-.]+)\s*(?P<x_gaussian_program_release_date>[0-9][0-9]?\-[A-Z][a-z][a-z]\-[0-9]+)'
            r'\s*(?P<x_gaussian_program_execution_date>[0-9][0-9]?\-[A-Z][a-z][a-z]\-[0-9]+)')
    ),
    LegacyParser(
        name='parsers/quantumespresso',
        parser_class_name='quantumespressoparser.QuantumEspressoParserPWSCF',
        main_file_re=r'^.*\.out$',  # It's not clear what type of extensions we should handler (.log?), either *star or .log,
        main_contents_re=r"^PWSCF$"
    )
]

""" Instanciation and constructor based config of all parsers. """

parser_dict = {parser.name: parser for parser in parsers}  # type: ignore
""" A dict to access parsers by name. Usually 'parsers/<...>', e.g. 'parsers/vasp'. """
