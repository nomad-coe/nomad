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

'''
The *parsing* module is an interface for the existing NOMAD-coe parsers.
This module redefines some of the old NOMAD-coe python-common functionality to create a
more coherent interface to the parsers.

Assumption about parsers
------------------------
For now, we make a few assumption about parsers
- they always work on the same *meta-info* version
- they have no conflicting python requirements
- they can be loaded at the same time and can be used within the same python process
- they are uniquely identified by a GIT URL and publicly accessible
- their version is uniquely identified by a GIT commit SHA

Each parser is defined via an instance of :class:`Parser`. The implementation :class:`LegacyParser` is used for most NOMAD-coe parsers.

.. autoclass:: nomad.parsing.Parser
    :members:

The are sub-classes for parsers with special purposes.

.. autoclass:: nomad.parsing.Parser
.. autoclass:: nomad.parsing.MatchingParser
.. autoclass:: nomad.parsing.MissingParser
.. autoclass:: nomad.parsing.BrokenParser
.. autoclass:: nomad.parsing.TemplateParser
.. autoclass:: nomad.parsing.GenerateRandomParser
.. autoclass:: nomad.parsing.ChaosParser
.. autoclass:: nomad.parsing.EmptyParser


The implementation :class:`LegacyParser` is used for most NOMAD-coe parsers.

.. autoclass:: nomad.parsing.LegacyParser


The parser definitions are available via the following two variables.

.. autodata:: nomad.parsing.parsers
.. autodata:: nomad.parsing.parser_dict

Parsers are reused for multiple calculations.

Parsers and calculation files are matched via regular expressions.

.. autofunction:: nomad.parsing.match_parser

Parsers in NOMAD-coe use a *backend* to create output. There are different NOMAD-coe
basends. In nomad@FAIRDI, we only currently only use a single backed. The following
classes provide a interface definition for *backends* as an ABC and a concrete implementation
based on nomad@fairdi's metainfo:

.. autoclass:: nomad.parsing.AbstractParserBackend
    :members:
.. autoclass:: nomad.parsing.Backend
    :members:
'''

from typing import Callable, IO, Union, Dict
import os.path

from nomad import config, datamodel

from nomad.parsing.legacy import (
    AbstractParserBackend, Backend, BackendError, BadContextUri, LegacyParser, VaspOutcarParser)
from nomad.parsing.parser import Parser, BrokenParser, MissingParser, MatchingParser
from nomad.parsing.artificial import (
    TemplateParser, GenerateRandomParser, ChaosParser, EmptyParser)

try:
    # these packages are not available without parsing extra, which is ok, if the
    # parsers are only initialized to load their metainfo definitions
    import magic
    import gzip
    import bz2
    import lzma

    _compressions = {
        b'\x1f\x8b\x08': ('gz', gzip.open),
        b'\x42\x5a\x68': ('bz2', bz2.open),
        b'\xfd\x37\x7a': ('xz', lzma.open)
    }

    encoding_magic = magic.Magic(mime_encoding=True)

except ImportError:
    pass


def match_parser(mainfile_path: str, strict=True) -> 'Parser':
    '''
    Performs parser matching. This means it take the given mainfile and potentially
    opens it with the given callback and tries to identify a parser that can parse
    the file.

    This is determined by filename (e.g. *.out), mime type (e.g. text/*, application/xml),
    and beginning file contents.

    Arguments:
        mainfile_path: Path to the mainfile
        strict: Only match strict parsers, e.g. no artificial parsers for missing or empty entries.

    Returns: The parser, or None if no parser could be matched.
    '''
    mainfile = os.path.basename(mainfile_path)
    if mainfile.startswith('.') or mainfile.startswith('~'):
        return None

    with open(mainfile_path, 'rb') as f:
        compression, open_compressed = _compressions.get(f.read(3), (None, open))

    with open_compressed(mainfile_path, 'rb') as cf:  # type: ignore
        buffer = cf.read(config.parser_matching_size)

    mime_type = magic.from_buffer(buffer, mime=True)

    decoded_buffer = None
    encoding = None
    try:  # Try to open the file as a string for regex matching.
        decoded_buffer = buffer.decode('utf-8')
    except UnicodeDecodeError:
        # This file is either binary or has wrong encoding
        encoding = encoding_magic.from_buffer(buffer)

        if config.services.force_raw_file_decoding:
            encoding = 'iso-8859-1'

        if encoding in ['iso-8859-1']:
            try:
                decoded_buffer = buffer.decode(encoding)
            except Exception:
                pass

    for parser in parsers:
        if strict and isinstance(parser, (MissingParser, EmptyParser)):
            continue

        if parser.is_mainfile(mainfile_path, mime_type, buffer, decoded_buffer, compression):
            # potentially convert the file
            if encoding in ['iso-8859-1']:
                try:
                    with open(mainfile_path, 'rb') as binary_file:
                        content = binary_file.read().decode(encoding)
                except Exception:
                    pass
                else:
                    with open(mainfile_path, 'wt') as text_file:
                        text_file.write(content)

            # TODO: deal with multiple possible parser specs
            return parser

    return None


parsers = [
    GenerateRandomParser(),
    TemplateParser(),
    ChaosParser(),
    LegacyParser(
        name='parsers/phonopy', code_name='Phonopy', code_homepage='https://phonopy.github.io/phonopy/',
        parser_class_name='phonopyparser.PhonopyParserWrapper',
        # mainfile_contents_re=r'',  # Empty regex since this code calls other DFT codes.
        mainfile_name_re=(r'.*/phonopy-FHI-aims-displacement-0*1/control.in$')
    ),
    LegacyParser(
        name='parsers/vasp', code_name='VASP', code_homepage='https://www.vasp.at/',
        parser_class_name='vaspparser.VASPRunParser',
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_contents_re=(
            r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
            r'?\s*<modeling>'
            r'?\s*<generator>'
            r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
            r'?'),
        supported_compressions=['gz', 'bz2', 'xz']
    ),
    VaspOutcarParser(
        name='parsers/vasp-outcar', code_name='VASP', code_homepage='https://www.vasp.at/',
        parser_class_name='vaspparser.VaspOutcarParser',
        mainfile_name_re=r'(.*/)?OUTCAR(\.[^\.]*)?',
        mainfile_contents_re=(r'^\svasp\.')
    ),
    LegacyParser(
        name='parsers/exciting', code_name='exciting', code_homepage='http://exciting-code.org/',
        parser_class_name='excitingparser.ExcitingParser',
        mainfile_name_re=r'^.*.OUT(\.[^/]*)?$',
        mainfile_contents_re=(r'EXCITING.*started')
    ),
    LegacyParser(
        name='parsers/fhi-aims', code_name='FHI-aims', code_homepage='https://aimsclub.fhi-berlin.mpg.de/',
        parser_class_name='fhiaimsparser.FHIaimsParser',
        mainfile_contents_re=(
            r'^(.*\n)*'
            r'?\s*Invoking FHI-aims \.\.\.'
            # r'?\s*Version'
        )
    ),
    LegacyParser(
        name='parsers/cp2k', code_name='CP2K', code_homepage='https://www.cp2k.org/',
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
        name='parsers/crystal', code_name='Crystal', code_homepage='https://www.crystal.unito.it/',
        parser_class_name='crystalparser.CrystalParser',
        mainfile_contents_re=(
            r'(CRYSTAL\s*\n\d+ \d+ \d+)|(CRYSTAL will run on \d+ processors)|'
            r'(\s*\*\s*CRYSTAL[\d]+\s*\*\s*\*\s*(public|Release) \: [\d\.]+.*\*)|'
            r'(Executable:\s*[/_\-a-zA-Z0-9]*MPPcrystal)'
        )
    ),
    # The main contents regex of CPMD was causing a catostrophic backtracking issue
    # when searching through the first 500 bytes of main files. We decided
    # to use only a portion of the regex to avoid that issue.
    LegacyParser(
        name='parsers/cpmd', code_name='CPMD', code_homepage='https://www.lcrc.anl.gov/for-users/software/available-software/cpmd/',
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
        name='parsers/nwchem', code_name='NWChem', code_homepage='http://www.nwchem-sw.org/',
        parser_class_name='nwchemparser.NWChemParser',
        mainfile_contents_re=(
            r'Northwest Computational Chemistry Package \(NWChem\) (\d+\.)+\d+'
        )
    ),
    LegacyParser(
        name='parsers/bigdft', code_name='BigDFT', code_homepage='http://bigdft.org/',
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
        name='parsers/wien2k', code_name='WIEN2k', code_homepage='http://www.wien2k.at/',
        parser_class_name='wien2kparser.Wien2kParser',
        mainfile_contents_re=r'\s*---------\s*:ITE[0-9]+:\s*[0-9]+\.\s*ITERATION\s*---------'
    ),
    LegacyParser(
        name='parsers/band', code_name='BAND', code_homepage='https://www.scm.com/product/band_periodicdft/',
        parser_class_name='bandparser.BANDParser',
        mainfile_contents_re=r' +\* +Amsterdam Density Functional +\(ADF\)'),
    LegacyParser(
        name='parsers/gaussian', code_name='Gaussian', code_homepage='http://gaussian.com/',
        parser_class_name='gaussianparser.GaussianParser',
        mainfile_mime_re=r'.*',
        mainfile_contents_re=(
            r'\s*Cite this work as:'
            r'\s*Gaussian [0-9]+, Revision [A-Za-z0-9\.]*,')
    ),
    LegacyParser(
        name='parsers/quantumespresso', code_name='Quantum Espresso', code_homepage='https://www.quantum-espresso.org/',
        parser_class_name='quantumespressoparser.QuantumEspressoParserPWSCF',
        mainfile_contents_re=(
            r'(Program PWSCF.*starts)|'
            r'(Current dimensions of program PWSCF are)')
        #    r'^(.*\n)*'
        #    r'\s*Program (\S+)\s+v\.(\S+)(?:\s+\(svn\s+rev\.\s+'
        #    r'(\d+)\s*\))?\s+starts[^\n]+'
        #    r'(?:\s*\n?)*This program is part of the open-source Quantum')
    ),
    LegacyParser(
        name='parsers/abinit', code_name='ABINIT', code_homepage='https://www.abinit.org/',
        parser_class_name='abinitparser.AbinitParser',
        mainfile_contents_re=(r'^\n*\.Version\s*[0-9.]*\s*of ABINIT\s*')
    ),
    LegacyParser(
        name='parsers/orca', code_name='ORCA', code_homepage='https://orcaforum.kofo.mpg.de/',
        parser_class_name='orcaparser.OrcaParser',
        mainfile_contents_re=(
            r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s+\* O   R   C   A \*\s*'
            r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s*'
            r'\s*--- An Ab Initio, DFT and Semiempirical electronic structure package ---\s*')
    ),
    LegacyParser(
        name='parsers/castep', code_name='CASTEP', code_homepage='http://www.castep.org/',
        parser_class_name='castepparser.CastepParser',
        mainfile_contents_re=(r'\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*')
    ),
    LegacyParser(
        name='parsers/dl-poly', code_name='DL_POLY', code_homepage='https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx',
        parser_class_name='dlpolyparser.DlPolyParserWrapper',
        mainfile_contents_re=(r'\*\* DL_POLY \*\*')
    ),
    LegacyParser(
        name='parsers/lib-atoms', code_name='libAtoms', code_homepage='https://libatoms.github.io/',
        parser_class_name='libatomsparser.LibAtomsParserWrapper',
        mainfile_contents_re=(r'\s*<GAP_params\s')
    ),
    LegacyParser(
        name='parsers/octopus', code_name='Octopus', code_homepage='https://octopus-code.org/',
        parser_class_name='octopusparser.OctopusParserWrapper',
        mainfile_contents_re=(r'\|0\) ~ \(0\) \|')
        # We decided to use the octopus eyes instead of
        # r'\*{32} Grid \*{32}Simulation Box:' since it was so far down in the file.
    ),
    # match gpaw2 first, other .gpw files are then considered to be "gpaw1"
    LegacyParser(
        name='parsers/gpaw2', code_name='GPAW', code_homepage='https://wiki.fysik.dtu.dk/gpaw/',
        parser_class_name='gpawparser.GPAWParser2Wrapper',
        mainfile_binary_header=b'GPAW',
        mainfile_name_re=(r'^.*\.(gpw2|gpw)$'),
        mainfile_mime_re=r'application/(x-tar|octet-stream)'
    ),
    LegacyParser(
        name='parsers/gpaw', code_name='GPAW', code_homepage='https://wiki.fysik.dtu.dk/gpaw/',
        parser_class_name='gpawparser.GPAWParserWrapper',
        mainfile_name_re=(r'^.*\.gpw$'),
        mainfile_mime_re=r'application/(x-tar|octet-stream)'
    ),
    LegacyParser(
        name='parsers/atk', code_name='ATK', code_homepage='https://www.synopsys.com/silicon/quantumatk.html',
        parser_class_name='atkparser.ATKParserWrapper',
        # mainfile_contents_re=r'',  # We can't read .gpw as txt - of UlmGPAW|AFFormatGPAW'
        mainfile_name_re=r'^.*\.nc',
        # The previously used mime type r'application/x-netcdf' wasn't found by magic library.
        mainfile_mime_re=r'application/octet-stream'
    ),
    LegacyParser(
        name='parsers/gulp', code_name='gulp', code_homepage='http://gulp.curtin.edu.au/gulp/',
        parser_class_name='gulpparser.GULPParser',
        mainfile_contents_re=(
            r'\s*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*'
            r'\*\*\*\*\*\*\*\*\*\*\*\*\*\s*'
            r'\s*\*\s*GENERAL UTILITY LATTICE PROGRAM\s*\*\s*')
    ),
    LegacyParser(
        name='parsers/siesta', code_name='Siesta', code_homepage='https://departments.icmab.es/leem/siesta/',
        parser_class_name='siestaparser.SiestaParser',
        mainfile_contents_re=(
            r'(Siesta Version: siesta-|SIESTA [0-9]\.[0-9]\.[0-9])|'
            r'(\*\s*WELCOME TO SIESTA\s*\*)')
    ),
    LegacyParser(
        name='parsers/elk', code_name='elk', code_homepage='http://elk.sourceforge.net/',
        parser_class_name='elkparser.ElkParser',
        mainfile_contents_re=r'\| Elk version [0-9.a-zA-Z]+ started \|'
    ),
    LegacyParser(
        name='parsers/elastic', code_name='elastic', code_homepage='http://exciting-code.org/elastic',
        parser_class_name='elasticparser.ElasticParser',
        mainfile_contents_re=r'\s*Order of elastic constants\s*=\s*[0-9]+\s*'
    ),
    LegacyParser(
        name='parsers/gamess', code_name='GAMESS', code_homepage='https://www.msg.chem.iastate.edu/gamess/versions.html',
        parser_class_name='gamessparser.GamessParser',
        mainfile_contents_re=(
            r'\s*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s*\*\s*GAMESS VERSION =\s*(.*)\*\s*'
            r'\s*\*\s*FROM IOWA STATE UNIVERSITY\s*\*\s*')
    ),
    LegacyParser(
        name='parsers/turbomole', code_name='turbomole', code_homepage='https://www.turbomole.org/',
        parser_class_name='turbomoleparser.TurbomoleParser',
        mainfile_contents_re=(
            r'Copyright \(C\) [0-9]+ TURBOMOLE GmbH, Karlsruhe')
    ),
    LegacyParser(
        name='parsers/skeleton', code_name='skeleton', code_homepage=None,
        domain='ems',
        parser_class_name='skeletonparser.SkeletonParserInterface',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_contents_re=(r'skeleton experimental metadata format')
    ),
    LegacyParser(
        name='parsers/mpes', code_name='mpes', code_homepage='https://github.com/mpes-kit/mpes', domain='ems',
        parser_class_name='mpesparser.MPESParserInterface',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=(r'.*.meta'),
        mainfile_contents_re=(r'"data_repository_name": "zenodo.org"')
    ),
    LegacyParser(
        name='parsers/aptfim', code_name='mpes', code_homepage='https://github.com/mpes-kit/mpes', domain='ems',
        parser_class_name='aptfimparser.APTFIMParserInterface',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=(r'.*.aptfim')
    ),
    LegacyParser(
        name='parsers/eels', code_name='eels', code_homepage='https://eelsdb.eu/', domain='ems',
        parser_class_name='eelsparser.EelsParserInterface',
        mainfile_mime_re=r'text/.*',
        mainfile_name_re=(r'.*.txt'),
        mainfile_contents_re=(r'api_permalink = https://api\.eelsdb\.eu')
    ),
    LegacyParser(
        name='parsers/qbox', code_name='qbox', code_homepage='http://qboxcode.org/', domain='dft',
        parser_class_name='qboxparser.QboxParser',
        mainfile_mime_re=r'(application/xml)|(text/.*)',
        mainfile_contents_re=(r'http://qboxcode.org')
    ),
    LegacyParser(
        name='parsers/dmol', code_name='DMol3', code_homepage='http://dmol3.web.psi.ch/dmol3.html', domain='dft',
        parser_class_name='dmol3parser.Dmol3Parser',
        mainfile_name_re=r'.*\.outmol',
        mainfile_contents_re=r'Materials Studio DMol\^3'
    ),
    LegacyParser(
        name='parsers/fleur', code_name='fleur', code_homepage='https://www.flapw.de/', domain='dft',
        parser_class_name='fleurparser.FleurParser',
        mainfile_contents_re=r'This output is generated by fleur.'
    ),
    LegacyParser(
        name='parsers/molcas', code_name='MOLCAS', code_homepage='http://www.molcas.org/', domain='dft',
        parser_class_name='molcasparser.MolcasParser',
        mainfile_contents_re=r'M O L C A S'
    ),
    LegacyParser(
        name='parsers/onetep', code_name='ONETEP', code_homepage='https://www.onetep.org/', domain='dft',
        parser_class_name='onetepparser.OnetepParser',
        mainfile_contents_re=r'####### #     # ####### ####### ####### ######'
    )
]

empty_parsers = [
    EmptyParser(
        name='missing/octopus', code_name='Octopus', code_homepage='https://octopus-code.org/',
        domain='dft',
        mainfile_name_re=r'(inp)|(.*/inp)'
    ),
    EmptyParser(
        name='missing/crystal', code_name='Crystal', code_homepage='https://www.crystal.unito.it/index.php',
        domain='dft',
        mainfile_name_re=r'.*\.cryst\.out'
    ),
    EmptyParser(
        name='missing/wien2k', code_name='WIEN2k', code_homepage='http://www.wien2k.at/',
        domain='dft',
        mainfile_name_re=r'.*\.scf'
    ),
    EmptyParser(
        name='missing/fhi-aims', code_name='FHI-aims', code_homepage='https://aimsclub.fhi-berlin.mpg.de/',
        domain='dft',
        mainfile_name_re=r'.*\.fhiaims'
    )
]

if config.use_empty_parsers:
    # There are some entries with PIDs that have mainfiles which do not match what
    # the actual parsers expect. We use the EmptyParser to produce placeholder entries
    # to keep the PIDs. These parsers will not match for new, non migrated data.
    parsers.extend(empty_parsers)

parsers.append(BrokenParser())

''' Instantiation and constructor based config of all parsers. '''

parser_dict = {parser.name: parser for parser in parsers + empty_parsers}  # type: ignore
''' A dict to access parsers by name. Usually 'parsers/<...>', e.g. 'parsers/vasp'. '''

# renamed parsers
parser_dict['parser/broken'] = parser_dict['parsers/broken']
parser_dict['parser/fleur'] = parser_dict['parsers/fleur']
parser_dict['parser/molcas'] = parser_dict['parsers/molcas']
parser_dict['parser/octopus'] = parser_dict['parsers/octopus']
parser_dict['parser/onetep'] = parser_dict['parsers/onetep']

# register code names as possible statistic value to the dft datamodel
code_names = sorted(
    set([
        getattr(parser, 'code_name')
        for parser in parsers
        if parser.domain == 'dft' and getattr(parser, 'code_name', None) is not None and getattr(parser, 'code_name') != 'currupted mainfile']),
    key=lambda code_name: code_name.lower())
datamodel.DFTMetadata.code_name.a_search.statistic_values = code_names + [config.services.unavailable_value, config.services.not_processed_value]
