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
- their version is uniquely identified by a GIT commit SHA

Each parser is defined via an instance of :class:`Parser`.

.. autoclass:: nomad.parsing.Parser
    :members:

The implementation :class:`LegacyParser` is used for most NOMAD-coe parsers.

.. autoclass:: nomad.parsing.LegacyParser

The parser definitions are available via the following two variables.

.. autodata:: nomad.parsing.parsers
.. autodata:: nomad.parsing.parser_dict

Parsers are reused for multiple caclulations.

Parsers and calculation files are matched via regular expressions.

.. autofunc:: nomad.parsing.match_parser

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
from typing import Callable, IO
import magic
import gzip
import bz2

from nomad import files

from nomad.parsing.backend import AbstractParserBackend, LocalBackend, LegacyLocalBackend, JSONStreamWriter, BadContextURI, WrongContextState
from nomad.parsing.parser import Parser, LegacyParser, VaspOutcarParser
from nomad.parsing.artificial import TemplateParser, GenerateRandomParser, ChaosParser


_compressions = {
    b'\x1f\x8b\x08': ('gz', gzip.open),
    b'\x42\x5a\x68': ('bz2', bz2.open)
}


def match_parser(mainfile: str, upload_files: files.StagingUploadFiles) -> 'Parser':
    """
    Performs parser matching. This means it take the given mainfile and potentially
    opens it with the given callback and tries to identify a parser that can parse
    the file.

    This is determined by filename (e.g. *.out), mime type (e.g. text/*, application/xml),
    and beginning file contents.

    Arguments:
        mainfile: The upload relative path to the mainfile
        open: A function that allows to open a stream to the file

    Returns: The parser, or None if no parser could be matched.
    """
    with upload_files.raw_file(mainfile, 'rb') as f:
        compression, open_compressed = _compressions.get(f.read(3), (None, open))

    mainfile_path = upload_files.raw_file_object(mainfile).os_path
    with open_compressed(mainfile_path, 'rb') as f:
        buffer = f.read(2048)

    mime_type = magic.from_buffer(buffer, mime=True)
    if mime_type.startswith('application') and not mime_type.endswith('xml'):
        return None

    for parser in parsers:
        if parser.is_mainfile(mainfile_path, mime_type, buffer.decode('utf-8'), compression):
            # TODO: deal with multiple possible parser specs
            return parser

    return None


parsers = [
    GenerateRandomParser(),
    TemplateParser(),
    ChaosParser(),
    LegacyParser(
        name='parsers/vasp',
        parser_class_name='vaspparser.VASPRunParserInterface',
        mainfile_mime_re=r'(application/xml)|(text/.*)',
        mainfile_contents_re=(
            r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
            r'?\s*<modeling>'
            r'?\s*<generator>'
            r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
            r'?'),
        supported_compressions=['gz', 'bz2']
    ),
    VaspOutcarParser(
        name='parsers/vasp',
        parser_class_name='vaspparser.VaspOutcarParser',
        mainfile_name_re=r'(.*/)?OUTCAR(\.[^\.]*)?',
        mainfile_contents_re=(r'^\svasp\.')
    ),
    LegacyParser(
        name='parsers/exciting',
        parser_class_name='excitingparser.ExcitingParser',
        mainfile_name_re=r'^.*/INFO\.OUT?',
        mainfile_contents_re=(
            r'^\s*=================================================+\s*'
            r'\s*\|\s*EXCITING\s+\S+\s+started\s*='
            r'\s*\|\s*version hash id:\s*\S*\s*=')
    ),
    LegacyParser(
        name='parsers/fhi-aims',
        parser_class_name='fhiaimsparser.FHIaimsParser',
        mainfile_contents_re=(
            r'^(.*\n)*'
            r'?\s*Invoking FHI-aims \.\.\.'
            r'?\s*Version')
    ),
    LegacyParser(
        name='parsers/cp2k',
        parser_class_name='cp2kparser.CP2KParser',
        mainfile_contents_re=(
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
        mainfile_contents_re=(
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
        mainfile_contents_re=(
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
        mainfile_contents_re=(
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
        mainfile_contents_re=(
            # r'__________________________________ A fast and precise DFT wavelet code\s*'
            # r'\|     \|     \|     \|     \|     \|\s*'
            # r'\|     \|     \|     \|     \|     \|      BBBB         i       gggggg\s*'
            # r'\|_____\|_____\|_____\|_____\|_____\|     B    B               g\s*'
            # r'\|     \|  :  \|  :  \|     \|     \|    B     B        i     g\s*'
            # r'\|     \|-0\+--\|-0\+--\|     \|     \|    B    B         i     g        g\s*'
            r'\|_____\|__:__\|__:__\|_____\|_____\|___ BBBBB          i     g         g\s*'
            # r'\|  :  \|     \|     \|  :  \|     \|    B    B         i     g         g\s*'
            # r'\|--\+0-\|     \|     \|-0\+--\|     \|    B     B     iiii     g         g\s*'
            # r'\|__:__\|_____\|_____\|__:__\|_____\|    B     B        i      g        g\s*'
            # r'\|     \|  :  \|  :  \|     \|     \|    B BBBB        i        g      g\s*'
            # r'\|     \|-0\+--\|-0\+--\|     \|     \|    B        iiiii          gggggg\s*'
            # r'\|_____\|__:__\|__:__\|_____\|_____\|__BBBBB\s*'
            # r'\|     \|     \|     \|  :  \|     \|                           TTTTTTTTT\s*'
            # r'\|     \|     \|     \|--\+0-\|     \|  DDDDDD          FFFFF        T\s*'
            # r'\|_____\|_____\|_____\|__:__\|_____\| D      D        F        TTTT T\s*'
            # r'\|     \|     \|     \|  :  \|     \|D        D      F        T     T\s*'
            # r'\|     \|     \|     \|--\+0-\|     \|D         D     FFFF     T     T\s*'
            # r'\|_____\|_____\|_____\|__:__\|_____\|D___      D     F         T    T\s*'
            # r'\|     \|     \|  :  \|     \|     \|D         D     F          TTTTT\s*'
            # r'\|     \|     \|--\+0-\|     \|     \| D        D     F         T    T\s*'
            # r'\|_____\|_____\|__:__\|_____\|_____\|          D     F        T     T\s*'
            # r'\|     \|     \|     \|     \|     \|         D               T    T\s*'
            # r'\|     \|     \|     \|     \|     \|   DDDDDD       F         TTTT\s*'
            # r'\|_____\|_____\|_____\|_____\|_____\|______                    www\.bigdft\.org'
        )
    ),
    LegacyParser(
        name='parsers/wien2k',
        parser_class_name='wien2kparser.Wien2kParser',
        mainfile_contents_re=r':LABEL\d+: using WIEN2k_\d+\.\d+'
    ),
    LegacyParser(
        name='parsers/band',
        parser_class_name='bandparser.BANDParser',
        mainfile_contents_re=r' +\* +Amsterdam Density Functional +\(ADF\)'),
    LegacyParser(
        name='parsers/gaussian',
        parser_class_name='gaussianparser.GaussianParser',
        mainfile_contents_re=(
            r'\s*Cite this work as:'
            r'\s*Gaussian [0-9]+, Revision [A-Za-z0-9.]*,'
            r'\s\*\*\*\*\*\*\*\*\*\*\*\**'
            r'\s*Gaussian\s*([0-9]+):\s*([A-Za-z0-9-.]+)\s*([0-9][0-9]?\-[A-Z][a-z][a-z]\-[0-9]+)'
            r'\s*([0-9][0-9]?\-[A-Z][a-z][a-z]\-[0-9]+)')
    ),
    LegacyParser(
        name='parsers/quantumespresso',
        parser_class_name='quantumespressoparser.QuantumEspressoParserPWSCF',
        mainfile_contents_re=(
            r'^\s*Program (\S+)\s+v\.(\S+)(?:\s+\(svn\s+rev\.\s+'
            r'(\d+)\s*\))?\s+starts[^\n]+'
            r'(?:\s*\n?)*This program is part of the open-source Quantum')
    ),
    LegacyParser(
        name='parsers/abinit',
        parser_class_name='abinitparser.AbinitParser',
        mainfile_contents_re=(r'^\n\.Version\s*[0-9.]*\s*of ABINIT\s*')
    ),
    # LegacyParser(
    #     name='parsers/orca',
    #     parser_class_name='orcaparser.OrcaParser',
    #     mainfile_contents_re=(
    #         r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
    #         r'\s+\* O   R   C   A \*\s*'
    #         r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
    #         r'\s*'
    #         r'\s*--- An Ab Initio, DFT and Semiempirical electronic structure package ---\s*')
    # ),
    LegacyParser(
        name='parsers/castep',
        parser_class_name='castepparser.CastepParser',
        mainfile_contents_re=(r'\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*')
    )
]


""" Instanciation and constructor based config of all parsers. """

parser_dict = {parser.name: parser for parser in parsers}  # type: ignore
""" A dict to access parsers by name. Usually 'parsers/<...>', e.g. 'parsers/vasp'. """
