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

import os.path

from nomad import config, datamodel

from .parser import MissingParser, BrokenParser, Parser, ArchiveParser
from .artificial import EmptyParser, GenerateRandomParser, TemplateParser, ChaosParser

from eelsdbconverter import EELSApiJsonConverter
from mpesparser import MPESParser
from aptfimparser import APTFIMParser
from vaspparser import VASPParser
from phonopyparser import PhonopyParser
from elasticparser import ElasticParser
from lammpsparser import LammpsParser
from gromacsparser import GromacsParser
from crystalparser import CrystalParser
from fhiaimsparser import FHIAimsParser
from excitingparser import ExcitingParser
from abinitparser import AbinitParser
from quantumespressoparser import QuantumEspressoParser
from gaussianparser import GaussianParser
from gpawparser import GPAWParser
from octopusparser import OctopusParser
from orcaparser import OrcaParser
from cp2kparser import CP2KParser
from fhivibesparser import FHIVibesParser
from turbomoleparser import TurbomoleParser
from castepparser import CastepParser
from wien2kparser import Wien2kParser
from nwchemparser import NWChemParser
from lobsterparser import LobsterParser
from bandparser import BandParser
from amberparser import AmberParser
from asapparser import AsapParser
from bigdftparser import BigDFTParser
from cpmdparser import CPMDParser
from dftbplusparser import DFTBPlusParser
from dlpolyparser import DLPolyParser
from dmol3parser import Dmol3Parser
from elkparser import ElkParser
from fleurparser import FleurParser
from fploparser import FploParser
from gamessparser import GamessParser
from gromosparser import GromosParser
from gulpparser import GulpParser
from molcasparser import MolcasParser
from mopacparser import MopacParser
from namdparser import NAMDParser
from onetepparser import OnetepParser
from siestaparser import SiestaParser
from tinkerparser import TinkerParser
from charmmparser import CharmmParser
from libatomsparser import LibAtomsParser
from atkparser import ATKParser
from qboxparser import QboxParser
from openkimparser import OpenKIMParser
from openmxparser import OpenmxParser

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


def match_parser(mainfile_path: str, strict=True) -> Parser:
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
    PhonopyParser(),
    VASPParser(),
    ExcitingParser(),
    FHIAimsParser(),
    FHIVibesParser(),
    CP2KParser(),
    CrystalParser(),
    CPMDParser(),
    NWChemParser(),
    BigDFTParser(),
    Wien2kParser(),
    BandParser(),
    QuantumEspressoParser(),
    GaussianParser(),
    AbinitParser(),
    OrcaParser(),
    CastepParser(),
    DLPolyParser(),
    LibAtomsParser(),
    OctopusParser(),
    GPAWParser(),
    ATKParser(),
    GulpParser(),
    SiestaParser(),
    ElkParser(),
    ElasticParser(),
    GamessParser(),
    TurbomoleParser(),
    MPESParser(),
    APTFIMParser(),
    EELSApiJsonConverter(),
    QboxParser(),
    Dmol3Parser(),
    FleurParser(),
    MolcasParser(),
    OnetepParser(),
    OpenKIMParser(),
    TinkerParser(),
    LammpsParser(),
    AmberParser(),
    GromacsParser(),
    LobsterParser(),
    GromosParser(),
    NAMDParser(),
    CharmmParser(),
    DFTBPlusParser(),
    AsapParser(),
    FploParser(),
    MopacParser(),
    OpenmxParser(),
    ArchiveParser()
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
code_names = []
for parser in parsers:
    if parser.domain == 'dft' and \
            getattr(parser, 'code_name', None) is not None and \
            getattr(parser, 'code_name') != 'currupted mainfile' and \
            getattr(parser, 'code_name') != 'Template':
        code_names.append(getattr(parser, 'code_name'))
code_names = sorted(set(code_names), key=lambda code_name: code_name.lower())
datamodel.DFTMetadata.code_name.a_search.statistic_values = code_names + [
    config.services.unavailable_value, config.services.not_processed_value]
