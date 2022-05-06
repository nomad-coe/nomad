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

from nomad import config
from nomad.datamodel import EntryArchive, EntryMetadata, results
from nomad.datamodel.context import Context, ClientContext

from .parser import MissingParser, BrokenParser, Parser, ArchiveParser, MatchingParserInterface
from .artificial import EmptyParser, GenerateRandomParser, TemplateParser, ChaosParser
from .tabular import TabularDataParser


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


def run_parser(
        mainfile_path: str, parser: Parser, mainfile_keys: List[str] = None,
        logger=None, server_context: bool = False) -> List[EntryArchive]:
    '''
    Parses a file, given the path, the parser, and mainfile_keys, as returned by
    :func:`match_parser`, and returns the resulting EntryArchive objects. Parsers that have
    `create_children == False` (the most common case) will generate a list with a single entry,
    for parsers that create children the list will consist of the main entry followed by the
    child entries. The returned archive objects will have minimal metadata.
    '''
    directory = os.path.dirname(mainfile_path)
    if server_context:
        entry_archive = EntryArchive(m_context=ParserContext(directory))
    else:
        entry_archive = EntryArchive(m_context=ClientContext(local_dir=directory))

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


parsers = [
    GenerateRandomParser(),
    TemplateParser(),
    ChaosParser(),
    MatchingParserInterface(
        'electronicparsers.AbinitParser',
        name='parsers/abinit', code_name='ABINIT', code_homepage='https://www.abinit.org/',
        mainfile_contents_re=(r'^\n*\.Version\s*[0-9.]*\s*of ABINIT\s*')
    ),
    MatchingParserInterface(
        'electronicparsers.ATKParser',
        name='parsers/atk', code_name='AtomistixToolKit',
        code_homepage='https://www.synopsys.com/silicon/quantumatk.html',
        mainfile_name_re=r'^.*\.nc', mainfile_mime_re=r'application/octet-stream'
    ),
    MatchingParserInterface(
        'electronicparsers.BandParser',
        name='parsers/band', code_name='BAND',
        code_homepage='https://www.scm.com/product/band_periodicdft/',
        mainfile_contents_re=r' +\* +Amsterdam Density Functional +\(ADF\)'
    ),
    MatchingParserInterface(
        'electronicparsers.BigDFTParser',
        name='parsers/bigdft', code_name='BigDFT', code_homepage='http://bigdft.org/',
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
        name='parsers/castep', code_name='CASTEP', code_homepage='http://www.castep.org/',
        mainfile_contents_re=(r'\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*')
    ),
    MatchingParserInterface(
        'electronicparsers.CharmmParser',
        name='parsers/charmm', code_name='Charmm', domain='dft',
        mainfile_contents_re=r'\s*Chemistry\s*at\s*HARvard\s*Macromolecular\s*Mechanics\s*',
        mainfile_mime_re=r'text/.*'
    ),
    MatchingParserInterface(
        'electronicparsers.CP2KParser',
        name='parsers/cp2k', code_name='CP2K', code_homepage='https://www.cp2k.org/',
        mainfile_contents_re=(
            r'\*\*\*\* \*\*\*\* \*\*\*\*\*\*  \*\*  PROGRAM STARTED AT\s.*\n'
            r' \*\*\*\*\* \*\* \*\*\*  \*\*\* \*\*   PROGRAM STARTED ON\s*.*\n'
            r' \*\*    \*\*\*\*   \*\*\*\*\*\*    PROGRAM STARTED BY .*\n'
            r' \*\*\*\*\* \*\*    \*\* \*\* \*\*   PROGRAM PROCESS ID .*\n'
            r'  \*\*\*\* \*\*  \*\*\*\*\*\*\*  \*\*  PROGRAM STARTED IN .*\n')
    ),
    MatchingParserInterface(
        'electronicparsers.CPMDParser',
        name='parsers/cpmd', code_name='CPMD',
        code_homepage='https://www.lcrc.anl.gov/for-users/software/available-software/cpmd/',
        mainfile_contents_re=(r'\*\*\*       \*\*   \*\*\*  \*\* \*\*\*\* \*\*  \*\*   \*\*\*')
    ),
    MatchingParserInterface(
        'electronicparsers.CrystalParser',
        name='parsers/crystal',
        code_name='Crystal',
        code_homepage='https://www.crystal.unito.it/',
        mainfile_contents_re=(
            fr'(\r?\n \*\s+CRYSTAL[\d]+\s+\*\r?\n \*\s*[a-zA-Z]+ : \d+[\.\d+]*)')
    ),
    MatchingParserInterface(
        'electronicparsers.Dmol3Parser',
        name='parsers/dmol', code_name='DMol3',
        code_homepage='http://dmol3.web.psi.ch/dmol3.html', domain='dft',
        mainfile_name_re=r'.*\.outmol',
        mainfile_contents_re=r'Materials Studio DMol\^3'
    ),
    MatchingParserInterface(
        'electronicparsers.ElkParser',
        name='parsers/elk', code_name='elk', code_homepage='http://elk.sourceforge.net/',
        mainfile_contents_re=r'\| Elk version [0-9.a-zA-Z]+ started \|'
    ),
    MatchingParserInterface(
        'electronicparsers.ExcitingParser',
        name='parsers/exciting', code_name='exciting', code_homepage='http://exciting-code.org/',
        mainfile_name_re=r'^.*.OUT(\.[^/]*)?$', mainfile_contents_re=(r'EXCITING.*started')
    ),
    MatchingParserInterface(
        'electronicparsers.FHIAimsParser',
        name='parsers/fhi-aims', code_name='FHI-aims',
        code_homepage='https://aimsclub.fhi-berlin.mpg.de/',
        mainfile_contents_re=(
            r'^(.*\n)*'
            r'?\s*Invoking FHI-aims \.\.\.')
    ),
    MatchingParserInterface(
        'electronicparsers.FleurParser',
        name='parsers/fleur', code_name='fleur',
        code_homepage='https://www.flapw.de/', domain='dft',
        mainfile_contents_re=r'This output is generated by fleur.'
    ),
    MatchingParserInterface(
        'electronicparsers.FploParser',
        name='parsers/fplo', code_name='fplo', domain='dft',
        mainfile_contents_re=r'\s*\|\s*FULL-POTENTIAL LOCAL-ORBITAL MINIMUM BASIS BANDSTRUCTURE CODE\s*\|\s*',
        mainfile_mime_re=r'text/.*'
    ),
    MatchingParserInterface(
        'electronicparsers.GamessParser',
        name='parsers/gamess', code_name='GAMESS',
        code_homepage='https://www.msg.chem.iastate.edu/gamess/versions.html',
        mainfile_contents_re=(
            r'\s*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s*\*\s*GAMESS VERSION =\s*(.*)\*\s*'
            r'\s*\*\s*FROM IOWA STATE UNIVERSITY\s*\*\s*')
    ),
    MatchingParserInterface(
        'electronicparsers.GaussianParser',
        name='parsers/gaussian', code_name='Gaussian', code_homepage='http://gaussian.com/',
        mainfile_mime_re=r'.*', mainfile_contents_re=(
            r'\s*Cite this work as:'
            r'\s*Gaussian [0-9]+, Revision [A-Za-z0-9\.]*,')
    ),
    MatchingParserInterface(
        'electronicparsers.GPAWParser',
        name='parsers/gpaw', code_name='GPAW', code_homepage='https://wiki.fysik.dtu.dk/gpaw/',
        mainfile_name_re=(r'^.*\.(gpw2|gpw)$'),
        mainfile_mime_re=r'application/(x-tar|octet-stream)'
    ),
    MatchingParserInterface(
        'electronicparsers.MolcasParser',
        name='parsers/molcas', code_name='MOLCAS', code_homepage='http://www.molcas.org/',
        domain='dft', mainfile_contents_re=r'M O L C A S'
    ),
    MatchingParserInterface(
        'electronicparsers.NWChemParser',
        name='parsers/nwchem', code_name='NWChem', code_homepage='http://www.nwchem-sw.org/',
        mainfile_contents_re=(
            r'Northwest Computational Chemistry Package \(NWChem\) (\d+\.)+\d+')
    ),
    MatchingParserInterface(
        'electronicparsers.OctopusParser',
        name='parsers/octopus', code_name='Octopus', code_homepage='https://octopus-code.org/',
        mainfile_contents_re=(r'\|0\) ~ \(0\) \|')
    ),
    MatchingParserInterface(
        'electronicparsers.OnetepParser',
        name='parsers/onetep', code_name='ONETEP', code_homepage='https://www.onetep.org/',
        domain='dft', mainfile_contents_re=r'####### #     # ####### ####### ####### ######'
    ),
    MatchingParserInterface(
        'electronicparsers.OpenmxParser',
        name='parsers/openmx', code_name='OpenMX', code_homepage='http://www.openmx-square.org/',
        mainfile_mime_re=r'(text/.*)',
        mainfile_name_re=r'.*\.out$',
        mainfile_contents_re=(r'^\*{59}\s+\*{59}\s+This calculation was performed by OpenMX'),
    ),
    MatchingParserInterface(
        'electronicparsers.OrcaParser',
        name='parsers/orca', code_name='ORCA', code_homepage='https://orcaforum.kofo.mpg.de/',
        mainfile_contents_re=(
            r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s+\* O   R   C   A \*\s*'
            r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s*'
            r'\s*--- An Ab Initio, DFT and Semiempirical electronic structure package ---\s*')
    ),
    MatchingParserInterface(
        'electronicparsers.Psi4Parser',
        name='parsers/psi4', code_name='Psi4', code_homepage='https://psicode.org/',
        mainfile_contents_re=(r'Psi4: An Open-Source Ab Initio Electronic Structure Package')
    ),
    MatchingParserInterface(
        'electronicparsers.QBallParser',
        name='parsers/qball',
        code_name='qball',
        mainfile_contents_re='qball',
        supported_compressions=["gz", "bz2", "xz"]
    ),
    MatchingParserInterface(
        'electronicparsers.QboxParser',
        name='parsers/qbox', code_name='qbox', code_homepage='http://qboxcode.org/',
        domain='dft', mainfile_mime_re=r'(application/xml)|(text/.*)',
        mainfile_contents_re=(r'http://qboxcode.org')
    ),
    MatchingParserInterface(
        'electronicparsers.QuantumEspressoParser',
        name='parsers/quantumespresso', code_name='Quantum Espresso',
        code_homepage='https://www.quantum-espresso.org/',
        mainfile_contents_re=(r'(Program PWSCF.*starts)|(Current dimensions of program PWSCF are)')
    ),
    MatchingParserInterface(
        'electronicparsers.SiestaParser',
        name='parsers/siesta', code_name='Siesta',
        code_homepage='https://departments.icmab.es/leem/siesta/',
        mainfile_contents_re=(
            r'(Siesta Version: siesta-|SIESTA [0-9]\.[0-9]\.[0-9])|'
            r'(\*\s*WELCOME TO SIESTA\s*\*)')
    ),
    MatchingParserInterface(
        'electronicparsers.TurbomoleParser',
        name='parsers/turbomole', code_name='turbomole',
        code_homepage='https://www.turbomole.org/',
        mainfile_contents_re=(r'Copyright \(C\) [0-9]+ TURBOMOLE GmbH, Karlsruhe')
    ),
    MatchingParserInterface(
        'electronicparsers.VASPParser',
        name='parsers/vasp', code_name='VASP', code_homepage='https://www.vasp.at/',
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
        name='parsers/wien2k', code_name='WIEN2k', code_homepage='http://www.wien2k.at/',
        mainfile_name_re=r'.*\.scf$',
        mainfile_alternative=True,
        mainfile_contents_re=r'\s*---------\s*:ITE[0-9]+:\s*[0-9]+\.\s*ITERATION\s*---------'
    ),
    MatchingParserInterface(
        'electronicparsers.YamboParser',
        name='parsers/yambo', code_name='YAMBO', code_homepage='https://yambo-code.org/',
        mainfile_contents_re=(r'Build.+\s+http://www\.yambo-code\.org')
    ),
    MatchingParserInterface(
        'atomisticparsers.AmberParser',
        name='parsers/amber', code_name='Amber', domain='dft',
        mainfile_contents_re=r'\s*Amber\s[0-9]+\s[A-Z]+\s*[0-9]+'
    ),
    MatchingParserInterface(
        'atomisticparsers.AsapParser',
        name='parsers/asap', code_name='ASAP', domain='dft',
        mainfile_name_re=r'.*.traj$', mainfile_mime_re=r'application/octet-stream'
    ),
    MatchingParserInterface(
        'atomisticparsers.DFTBPlusParser',
        name='parsers/dftbplus', code_name='DFTB+', domain='dft',
        mainfile_contents_re=r'\|  DFTB\+',
        mainfile_mime_re=r'text/.*'
    ),
    MatchingParserInterface(
        'atomisticparsers.DLPolyParser',
        name='parsers/dl-poly', code_name='DL_POLY',
        code_homepage='https://www.scd.stfc.ac.uk/Pages/DL_POLY.aspx',
        mainfile_contents_re=(r'\*\* DL_POLY \*\*'),
    ),
    MatchingParserInterface(
        'atomisticparsers.GromacsParser',
        name='parsers/gromacs', code_name='Gromacs', code_homepage='http://www.gromacs.org/',
        domain='dft', mainfile_contents_re=r'gmx mdrun, (VERSION|version)'
    ),
    MatchingParserInterface(
        'atomisticparsers.GromosParser',
        name='parsers/gromos', code_name='Gromos', domain='dft',
        mainfile_contents_re=r'Bugreports to http://www.gromos.net'
    ),
    MatchingParserInterface(
        'atomisticparsers.GulpParser',
        name='parsers/gulp', code_name='gulp', code_homepage='http://gulp.curtin.edu.au/gulp/',
        mainfile_contents_re=(
            r'\s*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*'
            r'\*\*\*\*\*\*\*\*\*\*\*\*\*\s*'
            r'\s*\*\s*GENERAL UTILITY LATTICE PROGRAM\s*\*\s*'),
    ),
    MatchingParserInterface(
        'atomisticparsers.LammpsParser',
        name='parsers/lammps', code_name='LAMMPS', code_homepage='https://lammps.sandia.gov/',
        domain='dft', mainfile_contents_re=r'^LAMMPS'
    ),
    MatchingParserInterface(
        'atomisticparsers.LibAtomsParser',
        name='parsers/lib-atoms', code_name='libAtoms', code_homepage='https://libatoms.github.io/',
        mainfile_contents_re=(r'\s*<GAP_params\s'),
    ),
    MatchingParserInterface(
        'electronicparsers.MopacParser',
        name='parsers/mopac', code_name='MOPAC', domain='dft',
        mainfile_contents_re=r'\s*\*\*\s*MOPAC\s*([0-9a-zA-Z]*)\s*\*\*\s*',
        mainfile_mime_re=r'text/.*',
    ),
    MatchingParserInterface(
        'atomisticparsers.NAMDParser',
        name='parsers/namd', code_name='Namd', domain='dft',
        mainfile_contents_re=r'\s*Info:\s*NAMD\s*[0-9.]+\s*for\s*',
        mainfile_mime_re=r'text/.*',
    ),
    MatchingParserInterface(
        'atomisticparsers.OpenKIMParser',
        name='parsers/openkim', code_name='OpenKIM', domain='dft',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_contents_re=r'openkim|OPENKIM|OpenKIM'
    ),
    MatchingParserInterface(
        'atomisticparsers.TinkerParser',
        name='parsers/tinker', code_name='TINKER', domain='dft',
        mainfile_contents_re=r'TINKER  ---  Software Tools for Molecular Design'
    ),
    MatchingParserInterface(
        'workflowparsers.AFLOWParser',
        name='parsers/aflow', code_name='AFlow', code_homepage='http://www.aflowlib.org/',
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
        name='parsers/asr', code_name='ASR',
        code_homepage='https://asr.readthedocs.io/en/latest/index.html',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=r'.*archive_.*\.json',
        mainfile_contents_re=(r'"name": "ASR"')
    ),
    MatchingParserInterface(
        'workflowparsers.ElasticParser',
        name='parsers/elastic', code_name='elastic', code_homepage='http://exciting-code.org/elastic',
        mainfile_contents_re=r'\s*Order of elastic constants\s*=\s*[0-9]+\s*',
        mainfile_name_re=(r'.*/INFO_ElaStic')
    ),
    MatchingParserInterface(
        'workflowparsers.FHIVibesParser',
        name='parsers/fhi-vibes', code_name='FHI-vibes',
        code_homepage='https://vibes.fhi-berlin.mpg.de/',
        mainfile_name_re=(r'^.*\.(nc)$'), mainfile_mime_re=r'(application/x-hdf)',
        mainfile_binary_header_re=br'^\x89HDF'
    ),
    MatchingParserInterface(
        'workflowparsers.LobsterParser',
        name='parsers/lobster', code_name='LOBSTER',
        code_homepage='http://schmeling.ac.rwth-aachen.de/cohp/',
        mainfile_name_re=r'.*lobsterout$',
        mainfile_contents_re=(r'^LOBSTER\s*v[\d\.]+.*'),
    ),
    MatchingParserInterface(
        'workflowparsers.MPParser',
        name='parsers/mp', code_name='MaterialsProject',
        code_homepage='https://materialsproject.org',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=r'.*mp.+materials\.json',
        mainfile_contents_re=(r'"pymatgen_version":')
    ),
    MatchingParserInterface(
        'workflowparsers.PhonopyParser',
        name='parsers/phonopy', code_name='Phonopy', code_homepage='https://phonopy.github.io/phonopy/',
        mainfile_name_re=(r'(.*/phonopy-FHI-aims-displacement-0*1/control.in$)|(.*/phon.+yaml)')
    ),
    MatchingParserInterface(
        'eelsdbparser.EELSDBParser',
        name='parsers/eels', code_name='eels', code_homepage='https://eelsdb.eu/',
        domain='ems',
        mainfile_mime_re=r'application/json',
        mainfile_contents_re=(r'https://eelsdb.eu/spectra')
    ),
    MatchingParserInterface(
        'nexusparser.NexusParser',
        name='parsers/nexus', code_name='NEXUS', code_homepage='https://www.nexus.eu/',
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_name_re=(r'.*\.nxs'),
        supported_compressions=['gz', 'bz2', 'xz']
    ),
    TabularDataParser(),
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
for parser in parsers:
    if parser.domain == 'dft' and \
            getattr(parser, 'code_name', None) is not None and \
            getattr(parser, 'code_name') != 'currupted mainfile' and \
            getattr(parser, 'code_name') != 'Template':
        code_names.append(getattr(parser, 'code_name'))
code_names = sorted(set(code_names), key=lambda code_name: code_name.lower())
results.Simulation.program_name.a_elasticsearch[0].values = code_names + [
    config.services.unavailable_value, config.services.not_processed_value]


def import_all_parsers():
    '''
    Imports all the parsers. This will instantiate all parser metainfo as a side
    effect.
    '''
    for parser in parsers:
        if isinstance(parser, MatchingParserInterface):
            parser.import_parser_class()
