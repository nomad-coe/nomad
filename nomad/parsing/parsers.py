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
from typing import Tuple, List, Dict
from collections.abc import Iterable
import pkgutil
from pathlib import Path

from nomad import config
from nomad.datamodel import EntryArchive, EntryMetadata, results
from nomad.datamodel.context import Context, ClientContext

from .parser import MissingParser, BrokenParser, Parser, ArchiveParser, MatchingParserInterface
from .artificial import EmptyParser, GenerateRandomParser, TemplateParser, ChaosParser
from .tabular import TabularDataParser
from nomad.datamodel.metainfo.eln.elabftw_parser import ELabFTWParser


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


def match_parser(mainfile_path: str, strict=True, parser_name: str = None) -> Tuple[Parser, List[str]]:
    '''
    Performs parser matching. This means it take the given mainfile and potentially
    opens it with the given callback and tries to identify a parser that can parse
    the file.

    This is determined by filename (e.g. *.out), mime type (e.g. text/*, application/xml),
    and beginning file contents.

    Arguments:
        mainfile_path: Path to the mainfile
        strict: Only match strict parsers, e.g. no artificial parsers for missing or empty entries.
        parser_name: Optional, to force the matching to test only a specific parser

    Returns:
        A tuple (`parser`, `mainfile_keys`). The `parser` is the matched parser, and
        `mainfile_keys` defines the keys to use for child entries, if any. If there are
        no child entries, `mainfile_keys` will be None. If no parser matches, we return
        (None, None).
    '''
    mainfile = os.path.basename(mainfile_path)
    if mainfile.startswith('.') or mainfile.startswith('~'):
        return None, None

    with open(mainfile_path, 'rb') as f:
        compression, open_compressed = _compressions.get(f.read(3), (None, open))

    with open_compressed(mainfile_path, 'rb') as cf:  # type: ignore
        buffer = cf.read(config.process.parser_matching_size)

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
    if parser_name:
        parser = parser_dict.get(parser_name)
        assert parser is not None, f'parser by the name `{parser_name}` does not exist'
        parsers_to_check = [parser]
    else:
        parsers_to_check = parsers
    for parser in parsers_to_check:
        if strict and isinstance(parser, (MissingParser, EmptyParser)):
            continue

        match_result = parser.is_mainfile(mainfile_path, mime_type, buffer, decoded_buffer, compression)
        if match_result:
            if isinstance(match_result, Iterable):
                assert parser.creates_children, 'Illegal return value - parser does not specify `creates_children`'
                for mainfile_key in match_result:  # type: ignore
                    assert mainfile_key and type(mainfile_key) == str, (
                        f'Child keys must be strings, got {type(mainfile_key)}')
                mainfile_keys = sorted(match_result)  # type: ignore
            else:
                mainfile_keys = None

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
            return parser, mainfile_keys

    return None, None


class ParserContext(Context):
    def __init__(self, mainfile_dir):
        self._mainfile_dir = mainfile_dir

    def raw_file(self, path, *args, **kwargs):
        return open(os.path.join(self._mainfile_dir, path), *args, **kwargs)

    def raw_path_exists(self, path: str) -> bool:
        return os.path.exists(os.path.join(self._mainfile_dir, path))


def run_parser(
        mainfile_path: str, parser: Parser, mainfile_keys: List[str] = None,
        logger=None, server_context: bool = False,
        username: str = None, password: str = None) -> List[EntryArchive]:
    '''
    Parses a file, given the path, the parser, and mainfile_keys, as returned by
    :func:`match_parser`, and returns the resulting EntryArchive objects. Parsers that have
    `create_children == False` (the most common case) will generate a list with a single entry,
    for parsers that create children the list will consist of the main entry followed by the
    child entries. The returned archive objects will have minimal metadata.
    '''
    directory = os.path.dirname(mainfile_path)
    if server_context:
        # TODO this looks totally wrong. ParserContext is not a server context at all.
        # There should be three different context. Client, Server and Parser. Currently,
        # ClientContext seems to cover both the client and the local parser use-case.
        entry_archive = EntryArchive(m_context=ParserContext(directory))
    else:
        entry_archive = EntryArchive(
            m_context=ClientContext(local_dir=directory, username=username, password=password))

    entry_archive.metadata = EntryMetadata()
    entry_archive.metadata.mainfile = mainfile_path

    entry_archives = [entry_archive]
    if mainfile_keys:
        child_archives = {}
        for mainfile_key in mainfile_keys:
            child_archive = EntryArchive()
            child_metadata = child_archive.m_create(EntryMetadata)
            child_metadata.mainfile = mainfile_path
            child_metadata.mainfile_key = mainfile_key
            child_archives[mainfile_key] = child_archive
            entry_archives.append(child_archive)
        kwargs = dict(child_archives=child_archives)
    else:
        kwargs = {}

    cwd = os.getcwd()
    try:
        mainfile_path = os.path.abspath(mainfile_path)
        os.chdir(os.path.abspath(os.path.dirname(mainfile_path)))
        parser.parse(mainfile_path, entry_archive, logger=logger, **kwargs)
    except Exception as e:
        if logger:
            logger.error('parsing was not successful', exc_info=e)
        raise e
    finally:
        os.chdir(cwd)
    for entry_archive in entry_archives:
        if entry_archive.metadata.domain is None:
            entry_archive.metadata.domain = parser.domain
    return entry_archives


# Here we resolve the path of the installation directories of various parsers.
# Note that this should be done in a way that does not yet import the modules
# themselves (this takes a while). The parser modules are imported lazily later.
prefix_electronic = f'{Path(pkgutil.get_loader("electronicparsers").path).parent.absolute()}'  # type: ignore
prefix_atomistic = f'{Path(pkgutil.get_loader("atomisticparsers").path).parent.absolute()}'  # type: ignore
prefix_workflow = f'{Path(pkgutil.get_loader("workflowparsers").path).parent.absolute()}'  # type: ignore
prefix_database = f'{Path(pkgutil.get_loader("databaseparsers").path).parent.absolute()}'  # type: ignore
prefix_eels = f'{Path(pkgutil.get_loader("eelsdbparser").path).parent.absolute()}'  # type: ignore

parsers = [
    GenerateRandomParser(),
    TemplateParser(),
    ChaosParser(),
    MatchingParserInterface(
        'electronicparsers.AbinitParser',
        metadata_path=f'{prefix_electronic}/abinit/metadata.yaml',
        mainfile_contents_re=(r'^\n*\.Version\s*[0-9.]*\s*of ABINIT\s*')
    ),
    MatchingParserInterface(
        'electronicparsers.AMSParser',
        metadata_path=f'{prefix_electronic}/ams/metadata.yaml',
        mainfile_contents_re=r'\* +\| +A M S +\| +\*'
    ),
    MatchingParserInterface(
        'electronicparsers.ATKParser',
        metadata_path=f'{prefix_electronic}/atk/metadata.yaml',
        mainfile_name_re=r'^.*\.nc', mainfile_mime_re=r'application/octet-stream'
    ),
    MatchingParserInterface(
        'electronicparsers.BigDFTParser',
        metadata_path=f'{prefix_electronic}/bigdft/metadata.yaml',
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
    MatchingParserInterface(
        'electronicparsers.CastepParser',
        metadata_path=f'{prefix_electronic}/castep/metadata.yaml',
        mainfile_contents_re=(r'\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*')
    ),
    MatchingParserInterface(
        'electronicparsers.CharmmParser',
        metadata_path=f'{prefix_electronic}/charmm/metadata.yaml',
        mainfile_contents_re=r'\s*Chemistry\s*at\s*HARvard\s*Macromolecular\s*Mechanics\s*',
        mainfile_mime_re=r'text/.*'
    ),
    MatchingParserInterface(
        'electronicparsers.CP2KParser',
        metadata_path=f'{prefix_electronic}/cp2k/metadata.yaml',
        mainfile_contents_re=(
            r'\*\*\*\* \*\*\*\* \*\*\*\*\*\*  \*\*  PROGRAM STARTED AT\s.*\n'
            r' \*\*\*\*\* \*\* \*\*\*  \*\*\* \*\*   PROGRAM STARTED ON\s*.*\n'
            r' \*\*    \*\*\*\*   \*\*\*\*\*\*    PROGRAM STARTED BY .*\n'
            r' \*\*\*\*\* \*\*    \*\* \*\* \*\*   PROGRAM PROCESS ID .*\n'
            r'  \*\*\*\* \*\*  \*\*\*\*\*\*\*  \*\*  PROGRAM STARTED IN .*\n')
    ),
    MatchingParserInterface(
        'electronicparsers.CPMDParser',
        metadata_path=f'{prefix_electronic}/cpmd/metadata.yaml',
        mainfile_contents_re=(r'\*\*\*       \*\*   \*\*\*  \*\* \*\*\*\* \*\*  \*\*   \*\*\*')
    ),
    MatchingParserInterface(
        'electronicparsers.CrystalParser',
        metadata_path=f'{prefix_electronic}/crystal/metadata.yaml',
        mainfile_contents_re=(
            fr'(\r?\n \*\s+CRYSTAL[\d]+\s+\*\r?\n \*\s*[a-zA-Z]+ : \d+[\.\d+]*)')
    ),
    MatchingParserInterface(
        'electronicparsers.Dmol3Parser',
        metadata_path=f'{prefix_electronic}/dmol3/metadata.yaml',
        mainfile_name_re=r'.*\.outmol',
        mainfile_contents_re=r'Materials Studio DMol\^3'
    ),
    MatchingParserInterface(
        'electronicparsers.ElkParser',
        metadata_path=f'{prefix_electronic}/elk/metadata.yaml',
        mainfile_contents_re=r'\| Elk version [0-9.a-zA-Z]+ started \|'
    ),
    MatchingParserInterface(
        'electronicparsers.ExcitingParser',
        metadata_path=f'{prefix_electronic}/exciting/metadata.yaml',
        mainfile_name_re=r'^.*.OUT(\.[^/]*)?$',
        mainfile_contents_re=(r'EXCITING.*started[\s\S]+?All units are atomic ')
    ),
    MatchingParserInterface(
        'electronicparsers.FHIAimsParser',
        metadata_path=f'{prefix_electronic}/fhiaims/metadata.yaml',
        mainfile_contents_re=(
            r'^(.*\n)*'
            r'?\s*Invoking FHI-aims \.\.\.')
    ),
    MatchingParserInterface(
        'electronicparsers.FleurParser',
        metadata_path=f'{prefix_electronic}/fleur/metadata.yaml',
        mainfile_contents_re=r'This output is generated by fleur.|\<fleurOutput',
        mainfile_name_re=r'.*[^/]*\.xml[^/]*',  # only the alternative mainfile name should match
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_alternative=True
    ),
    MatchingParserInterface(
        'electronicparsers.FploParser',
        metadata_path=f'{prefix_electronic}/fplo/metadata.yaml',
        mainfile_contents_re=r'\s*\|\s*FULL-POTENTIAL LOCAL-ORBITAL MINIMUM BASIS BANDSTRUCTURE CODE\s*\|\s*',
        mainfile_mime_re=r'text/.*'
    ),
    MatchingParserInterface(
        'electronicparsers.GamessParser',
        metadata_path=f'{prefix_electronic}/gamess/metadata.yaml',
        mainfile_contents_re=(
            r'\s*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s*\*\s*GAMESS VERSION =\s*(.*)\*\s*'
            r'\s*\*\s*FROM IOWA STATE UNIVERSITY\s*\*\s*')
    ),
    MatchingParserInterface(
        'electronicparsers.GaussianParser',
        metadata_path=f'{prefix_electronic}/gaussian/metadata.yaml',
        mainfile_mime_re=r'.*',
        mainfile_contents_re=(
            r'\s*Cite this work as:'
            r'\s*Gaussian [0-9]+, Revision [A-Za-z0-9\.]*,')
    ),
    MatchingParserInterface(
        'electronicparsers.GPAWParser',
        metadata_path=f'{prefix_electronic}/gpaw/metadata.yaml',
        mainfile_name_re=(r'^.*\.(gpw2|gpw)$'),
        mainfile_mime_re=r'application/(x-tar|octet-stream)'
    ),
    MatchingParserInterface(
        'electronicparsers.MolcasParser',
        metadata_path=f'{prefix_electronic}/molcas/metadata.yaml',
        mainfile_contents_re=r'M O L C A S'
    ),
    MatchingParserInterface(
        'electronicparsers.MopacParser',
        metadata_path=f'{prefix_electronic}/mopac/metadata.yaml',
        mainfile_contents_re=r'\s*\*\*\s*MOPAC\s*([0-9a-zA-Z]*)\s*\*\*\s*',
        mainfile_mime_re=r'text/.*',
    ),
    MatchingParserInterface(
        'electronicparsers.NWChemParser',
        metadata_path=f'{prefix_electronic}/nwchem/metadata.yaml',
        mainfile_contents_re=(
            r'Northwest Computational Chemistry Package \(NWChem\) (\d+\.)+\d+')
    ),
    MatchingParserInterface(
        'electronicparsers.OctopusParser',
        metadata_path=f'{prefix_electronic}/octopus/metadata.yaml',
        mainfile_contents_re=(r'\|0\) ~ \(0\) \|')
    ),
    MatchingParserInterface(
        'electronicparsers.OnetepParser',
        metadata_path=f'{prefix_electronic}/onetep/metadata.yaml',
        mainfile_contents_re=r'####### #     # ####### ####### ####### ######'
    ),
    MatchingParserInterface(
        'electronicparsers.OpenmxParser',
        metadata_path=f'{prefix_electronic}/openmx/metadata.yaml',
        mainfile_mime_re=r'(text/.*)',
        mainfile_name_re=r'.*\.out$',
        mainfile_contents_re=(r'^\*{59}\s+\*{59}\s+This calculation was performed by OpenMX'),
    ),
    MatchingParserInterface(
        'electronicparsers.OrcaParser',
        metadata_path=f'{prefix_electronic}/orca/metadata.yaml',
        mainfile_contents_re=(
            r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s+\* O   R   C   A \*\s*'
            r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s*')
    ),
    MatchingParserInterface(
        'electronicparsers.Psi4Parser',
        metadata_path=f'{prefix_electronic}/psi4/metadata.yaml',
        mainfile_contents_re=(r'Psi4: An Open-Source Ab Initio Electronic Structure Package')
    ),
    MatchingParserInterface(
        'electronicparsers.QBallParser',
        metadata_path=f'{prefix_electronic}/qball/metadata.yaml',
        mainfile_contents_re='qball',
        supported_compressions=["gz", "bz2", "xz"]
    ),
    MatchingParserInterface(
        'electronicparsers.QboxParser',
        metadata_path=f'{prefix_electronic}/qbox/metadata.yaml',
        mainfile_mime_re=r'(application/xml)|(text/.*)',
        mainfile_contents_re=(r'http://qboxcode.org')
    ),
    MatchingParserInterface(
        'electronicparsers.QuantumEspressoParser',
        metadata_path=f'{prefix_electronic}/quantumespresso/metadata.yaml',
        mainfile_contents_re=(r'(Program PWSCF.*starts)|(Current dimensions of program PWSCF are)')
    ),
    MatchingParserInterface(
        'electronicparsers.SiestaParser',
        metadata_path=f'{prefix_electronic}/siesta/metadata.yaml',
        mainfile_contents_re=(
            r'(Siesta Version: siesta-|SIESTA [0-9]\.[0-9]\.[0-9])|'
            r'(\*\s*WELCOME TO SIESTA\s*\*)')
    ),
    MatchingParserInterface(
        'electronicparsers.TurbomoleParser',
        metadata_path=f'{prefix_electronic}/turbomole/metadata.yaml',
        mainfile_contents_re=(r'Copyright \(C\) [0-9]+ TURBOMOLE GmbH, Karlsruhe')
    ),
    MatchingParserInterface(
        'electronicparsers.VASPParser',
        metadata_path=f'{prefix_electronic}/vasp/metadata.yaml',
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_name_re=r'.*[^/]*xml[^/]*',  # only the alternative mainfile name should match
        mainfile_contents_re=(
            r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
            r'?\s*<modeling>'
            r'?\s*<generator>'
            r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
            r'?|^\svasp[\.\d]+.+?(?:\(build|complex)[\s\S]+?executed on'),
        supported_compressions=['gz', 'bz2', 'xz'], mainfile_alternative=True
    ),
    MatchingParserInterface(
        'electronicparsers.Wien2kParser',
        metadata_path=f'{prefix_electronic}/wien2k/metadata.yaml',
        mainfile_name_re=r'.*\.scf$',
        mainfile_alternative=True,
        mainfile_contents_re=r'\s*---------\s*:ITE[0-9]+:\s*[0-9]+\.\s*ITERATION\s*---------'
    ),
    MatchingParserInterface(
        'electronicparsers.YamboParser',
        metadata_path=f'{prefix_electronic}/yambo/metadata.yaml',
        mainfile_contents_re=(r'Build.+\s+http://www\.yambo-code\.org')
    ),
    MatchingParserInterface(
        'electronicparsers.ABACUSParser',
        metadata_path=f'{prefix_electronic}/abacus/metadata.yaml',
        mainfile_contents_re=(r'\s*\n\s*WELCOME TO ABACUS')
    ),
    MatchingParserInterface(
        'atomisticparsers.AmberParser',
        metadata_path=f'{prefix_atomistic}/amber/metadata.yaml',
        mainfile_contents_re=r'\s*Amber\s[0-9]+\s[A-Z]+\s*[0-9]+'
    ),
    MatchingParserInterface(
        'atomisticparsers.AsapParser',
        metadata_path=f'{prefix_atomistic}/asap/metadata.yaml',
        mainfile_name_re=r'.*.traj$', mainfile_mime_re=r'application/octet-stream',
        mainfile_binary_header_re=br'AFFormatASE\-Trajectory'
    ),
    MatchingParserInterface(
        'atomisticparsers.BOPfoxParser',
        metadata_path=f'{prefix_atomistic}/bopfox/metadata.yaml',
        mainfile_contents_re=r'\-+\s+BOPfox \(v'
    ),
    MatchingParserInterface(
        'atomisticparsers.DFTBPlusParser',
        metadata_path=f'{prefix_atomistic}/dftbplus/metadata.yaml',
        mainfile_contents_re=r'\|  DFTB\+',
        mainfile_mime_re=r'text/.*'
    ),
    MatchingParserInterface(
        'atomisticparsers.DLPolyParser',
        metadata_path=f'{prefix_atomistic}/dlpoly/metadata.yaml',
        mainfile_contents_re=(r'\*\*\s+DL_POLY.+\*\*'),
    ),
    MatchingParserInterface(
        'atomisticparsers.GromacsParser',
        metadata_path=f'{prefix_atomistic}/gromacs/metadata.yaml',
        mainfile_contents_re=r'gmx mdrun, (VERSION|version)[\s\S]*Input Parameters:'
    ),
    MatchingParserInterface(
        'atomisticparsers.GromosParser',
        metadata_path=f'{prefix_atomistic}/gromos/metadata.yaml',
        mainfile_contents_re=r'Bugreports to http://www.gromos.net'
    ),
    MatchingParserInterface(
        'atomisticparsers.GulpParser',
        metadata_path=f'{prefix_atomistic}/gulp/metadata.yaml',
        mainfile_contents_re=(
            r'\s*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*'
            r'\*\*\*\*\*\*\*\*\*\*\*\*\*\s*'
            r'\s*\*\s*GENERAL UTILITY LATTICE PROGRAM\s*\*\s*'),
    ),
    MatchingParserInterface(
        'atomisticparsers.LammpsParser',
        metadata_path=f'{prefix_atomistic}/lammps/metadata.yaml',
        mainfile_contents_re=r'^LAMMPS\s+\(.+\)'
    ),
    MatchingParserInterface(
        'atomisticparsers.LibAtomsParser',
        metadata_path=f'{prefix_atomistic}/libatoms/metadata.yaml',
        mainfile_contents_re=(r'\s*<GAP_params\s'),
    ),
    MatchingParserInterface(
        'atomisticparsers.NAMDParser',
        metadata_path=f'{prefix_atomistic}/namd/metadata.yaml',
        mainfile_contents_re=r'\s*Info:\s*NAMD\s*[0-9.]+\s*for\s*',
        mainfile_mime_re=r'text/.*',
    ),
    MatchingParserInterface(
        'atomisticparsers.TinkerParser',
        metadata_path=f'{prefix_atomistic}/tinker/metadata.yaml',
        mainfile_contents_re=r'TINKER  ---  Software Tools for Molecular Design'
    ),
    MatchingParserInterface(
        'atomisticparsers.XTBParser',
        metadata_path=f'{prefix_atomistic}/xtb/metadata.yaml',
        mainfile_contents_re=r'x T B\s+\|\s+\|\s+='
    ),
    MatchingParserInterface(
        'workflowparsers.AFLOWParser',
        metadata_path=f'{prefix_workflow}/aflow/metadata.yaml',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=r'.*aflowlib\.json.*',  # only the alternative mainfile name should match
        mainfile_contents_re=(
            r"^\s*\[AFLOW\] \*+"
            r"\s*\[AFLOW\]"
            r"\s*\[AFLOW\]                     .o.        .o88o. oooo"
            r"\s*\[AFLOW\]                    .888.       888 `` `888"
            r"\s*\[AFLOW\]                   .8'888.     o888oo   888   .ooooo.  oooo oooo    ooo"
            r"\s*\[AFLOW\]                  .8' `888.     888     888  d88' `88b  `88. `88.  .8'"
            r"\s*\[AFLOW\]                 .88ooo8888.    888     888  888   888   `88..]88..8'"
            r"\s*\[AFLOW\]                .8'     `888.   888     888  888   888    `888'`888'"
            r"\s*\[AFLOW\]               o88o     o8888o o888o   o888o `Y8bod8P'     `8'  `8'  .in"
            r"|^\s*\{\"aurl\"\:\"aflowlib\.duke\.edu\:AFLOWDATA"),
        supported_compressions=['gz', 'bz2', 'xz'], mainfile_alternative=True
    ),
    MatchingParserInterface(
        'workflowparsers.ASRParser',
        metadata_path=f'{prefix_workflow}/asr/metadata.yaml',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=r'.*archive_.*\.json',
        mainfile_contents_re=(r'"name": "ASR"')
    ),
    MatchingParserInterface(
        'workflowparsers.ElasticParser',
        metadata_path=f'{prefix_workflow}/elastic/metadata.yaml',
        mainfile_contents_re=r'\s*Order of elastic constants\s*=\s*[0-9]+\s*',
        mainfile_name_re=(r'.*/INFO_ElaStic')
    ),
    MatchingParserInterface(
        'workflowparsers.FHIVibesParser',
        metadata_path=f'{prefix_workflow}/fhivibes/metadata.yaml',
        mainfile_name_re=(r'^.*\.(nc)$'),
        mainfile_mime_re=r'(application/x-hdf)',
        mainfile_binary_header_re=br'^\x89HDF',
        mainfile_contents_dict={'__has_all_keys': ['I', 'a', 'b']}
    ),
    MatchingParserInterface(
        'workflowparsers.LobsterParser',
        metadata_path=f'{prefix_workflow}/lobster/metadata.yaml',
        mainfile_name_re=r'.*lobsterout$',
        mainfile_contents_re=(r'^LOBSTER\s*v[\d\.]+.*'),
    ),
    MatchingParserInterface(
        'workflowparsers.AtomateParser',
        metadata_path=f'{prefix_workflow}/atomate/metadata.yaml',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=r'.*mp.+materials\.json',
        mainfile_contents_re=(r'"pymatgen_version":')
    ),
    MatchingParserInterface(
        'workflowparsers.PhonopyParser',
        metadata_path=f'{prefix_workflow}/phonopy/metadata.yaml',
        mainfile_name_re=(r'(.*/phonopy-FHI-aims-displacement-0*1/control.in$)|(.*/phon[^/]+yaml)')
    ),
    MatchingParserInterface(
        'eelsdbparser.EELSDBParser',
        metadata_path=f'{prefix_eels}/metadata.yaml',
        mainfile_mime_re=r'application/json',
        mainfile_contents_re=(r'https://eelsdb.eu/spectra')
    ),
    MatchingParserInterface(
        'workflowparsers.MOFStructuresParser',
        metadata_path=f'{prefix_workflow}/mofstructures/metadata.yaml',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=r'.*mof_.*\.json',
        mainfile_contents_re=r'MOF Structures'
    ),
    MatchingParserInterface(
        'workflowparsers.QuantumEspressoPhononParser',
        metadata_path=f'{prefix_workflow}/quantum_espresso_phonon/metadata.yaml',
        mainfile_contents_re=(
            r'Program PHONON.+\s*'
            r'This program is part of the open-source Quantum ESPRESSO suite')
    ),
    MatchingParserInterface(
        'workflowparsers.QuantumEspressoEPWParser',
        metadata_path=f'{prefix_workflow}/quantum_espresso_epw/metadata.yaml',
        mainfile_contents_re=(
            r'Program EPW.+\s*'
            r'This program is part of the open-source Quantum ESPRESSO suite')
    ),
    MatchingParserInterface(
        'workflowparsers.QuantumEspressoXSpectraParser',
        metadata_path=f'{prefix_workflow}/quantum_espresso_xspectra/metadata.yaml',
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_contents_re=r'\s*Program XSpectra\s*'
    ),
    MatchingParserInterface(
        'databaseparsers.OpenKIMParser',
        metadata_path=f'{prefix_database}/openkim/metadata.yaml',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_contents_re=r'openkim|OPENKIM|OpenKIM'
    ),
    MatchingParserInterface(
        'electronicparsers.Wannier90Parser',
        metadata_path=f'{prefix_electronic}/wannier90/metadata.yaml',
        mainfile_contents_re=r'\|\s*WANNIER90\s*\|'
    ),
    MatchingParserInterface(
        'electronicparsers.W2DynamicsParser',
        metadata_path=f'{prefix_electronic}/w2dynamics/metadata.yaml',
        mainfile_name_re=(r'^.*\.(h5|hdf5)$'),
        mainfile_mime_re=r'(application/x-hdf)',
        mainfile_binary_header_re=br'^\x89HDF',
        mainfile_contents_dict={'__has_all_keys': ['.axes', '.config', '.quantities']}
    ),
    MatchingParserInterface(
        'electronicparsers.SolidDMFTParser',
        metadata_path=f'{prefix_electronic}/soliddmft/metadata.yaml',
        mainfile_name_re=(r'^.*\.(h5|hdf5)$'),
        mainfile_mime_re=r'(application/x-hdf)',
        mainfile_binary_header_re=br'^\x89HDF',
        mainfile_contents_dict={'__has_all_keys': ['dft_input', 'DMFT_input', 'DMFT_results']}
    ),
    MatchingParserInterface(
        'electronicparsers.OceanParser',
        metadata_path=f'{prefix_electronic}/ocean/metadata.yaml',
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_contents_dict={'__has_all_keys': ['bse', 'structure', 'screen', 'calc']}
    ),
    MatchingParserInterface(
        'nomad.parsing.nexus.NexusParser',
        metadata_path=os.path.join(os.path.dirname(__file__), 'metadata.yaml'),
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_name_re=r'.*\.nxs',
        supported_compressions=['gz', 'bz2', 'xz']),

    TabularDataParser(),
    ELabFTWParser(),
    ArchiveParser()
]

empty_parsers = [
    EmptyParser(
        name='missing/octopus', code_name='Octopus', code_homepage='https://octopus-code.org/',
        mainfile_name_re=r'(inp)|(.*/inp)'
    ),
    EmptyParser(
        name='missing/crystal', code_name='CRYSTAL', code_homepage='https://www.crystal.unito.it/index.php',
        mainfile_name_re=r'.*\.cryst\.out'
    ),
    EmptyParser(
        name='missing/wien2k', code_name='WIEN2k', code_homepage='http://www.wien2k.at/',
        mainfile_name_re=r'.*\.scf'
    ),
    EmptyParser(
        name='missing/fhi-aims', code_name='FHI-aims', code_homepage='https://aimsclub.fhi-berlin.mpg.de/',
        mainfile_name_re=r'.*\.fhiaims'
    )
]

if config.process.use_empty_parsers:
    # There are some entries with PIDs that have mainfiles which do not match what
    # the actual parsers expect. We use the EmptyParser to produce placeholder entries
    # to keep the PIDs. These parsers will not match for new, non migrated data.
    parsers.extend(empty_parsers)

parsers.append(BrokenParser())

''' Instantiation and constructor based config of all parsers. '''

parser_dict: Dict[str, Parser] = {parser.name: parser for parser in parsers + empty_parsers}  # type: ignore
''' A dict to access parsers by name. Usually 'parsers/<...>', e.g. 'parsers/vasp'. '''

# renamed parsers
parser_dict['parsers/vaspoutcar'] = parser_dict['parsers/vasp']
parser_dict['parser/broken'] = parser_dict['parsers/broken']
parser_dict['parser/fleur'] = parser_dict['parsers/fleur']
parser_dict['parser/molcas'] = parser_dict['parsers/molcas']
parser_dict['parser/octopus'] = parser_dict['parsers/octopus']
parser_dict['parser/onetep'] = parser_dict['parsers/onetep']

# register code names as possible statistic value to the dft datamodel
code_names = []
code_metadata = {}
for parser in parsers:
    code_name = getattr(parser, 'code_name', None)
    if parser.domain == 'dft' and \
            code_name is not None and \
            code_name != 'currupted mainfile' and \
            code_name != 'Template':
        code_names.append(code_name)
        if parser.metadata:
            code_metadata[code_name] = parser.metadata.dict()
        else:
            code_metadata[code_name] = {}
code_names = sorted(set(code_names), key=lambda code_name: code_name.lower())
results.Simulation.program_name.a_elasticsearch[0].values = code_names + [config.services.unavailable_value]


def import_all_parsers():
    '''
    Imports all the parsers. This will instantiate all parser metainfo as a side
    effect.
    '''
    for parser in parsers:
        if isinstance(parser, MatchingParserInterface):
            parser.import_parser_class()
