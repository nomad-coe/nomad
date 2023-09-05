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

'''
This module describes all configurable parameters for the nomad python code. The
configuration is used for all executed python code including API, worker, CLI, and other
scripts. To use the configuration in your own scripts or new modules, simply import
this module.

All parameters are structured into objects for two reasons. First, to have
categories. Second, to allow runtime manipulation that is not effected
by python import logic. The categories are choosen along infrastructure components:
``mongo``, ``elastic``, etc.

This module also provides utilities to read the configuration from environment variables
and .yaml files. This is done automatically on import. The precedence is env over .yaml
over defaults.
'''

import logging
import os
import os.path
import yaml
import warnings
from typing import List, Any, Union, get_type_hints
from pydantic import parse_obj_as, BaseModel

from .models import (
    NomadSettings, Services, Meta, Oasis, RabbitMQ, Celery, FS, Elastic, Keycloak,
    Mongo, Logstash, Logtransfer, Tests, Mail, Normalize, Resources, Client, DataCite, GitLab, Process,
    Reprocess, RFC3161Timestamp, BundleExport, BundleImport, Archive, UI
)
from .plugins import Plugins, Plugin, Parser, Schema
from .north import NORTH

warnings.filterwarnings('ignore', message='numpy.dtype size changed')
warnings.filterwarnings('ignore', message='numpy.ufunc size changed')
warnings.filterwarnings('ignore', category=DeprecationWarning)


def api_url(ssl: bool = True, api: str = 'api', api_host: str = None, api_port: int = None):
    '''
    Returns the url of the current running nomad API. This is for server-side use.
    This is not the NOMAD url to use as a client, use `nomad.config.client.url` instead.
    '''
    if api_port is None:
        api_port = services.api_port
    if api_host is None:
        api_host = services.api_host
    protocol = 'https' if services.https and ssl else 'http'
    host_and_port = api_host
    if api_port not in [80, 443]:
        host_and_port += ':' + str(api_port)
    base_path = services.api_base_path.strip('/')
    return f'{protocol}://{host_and_port}/{base_path}/{api}'


def gui_url(page: str = None):
    base = api_url(True)[:-3]
    if base.endswith('/'):
        base = base[:-1]

    if page is not None:
        return '%s/gui/%s' % (base, page)

    return '%s/gui' % base


def rabbitmq_url():
    return 'pyamqp://%s:%s@%s//' % (rabbitmq.user, rabbitmq.password, rabbitmq.host)


def north_url(ssl: bool = True):
    return api_url(ssl=ssl, api='north', api_host=north.hub_host, api_port=north.hub_port)


def hub_url():
    return f'http://{north.hub_host}:{north.hub_port}{services.api_base_path}/north/hub'


services = Services()
meta = Meta(deployment_url=api_url())
oasis = Oasis()
north = NORTH()
rabbitmq = RabbitMQ()
celery = Celery()
fs = FS()
elastic = Elastic()
keycloak = Keycloak()
mongo = Mongo()
logstash = Logstash()
logtransfer = Logtransfer()
tests = Tests()
mail = Mail()
normalize = Normalize()
resources = Resources()
client = Client()
datacite = DataCite()
gitlab = GitLab()
process = Process()
reprocess = Reprocess()
rfc3161_timestamp = RFC3161Timestamp()
bundle_export = BundleExport()
bundle_import = BundleImport()
archive = Archive()
ui = UI()
plugins = Plugins(options={
    'parsers/abinit': Parser(
        python_package='electronicparsers.abinit',
        mainfile_contents_re=r'^\n*\.Version\s*[0-9.]*\s*of ABINIT\s*'),
    'parsers/ams': Parser(
        python_package='electronicparsers.ams',
        mainfile_contents_re=r'\* +\| +A M S +\| +\*'),
    'parsers/atk': Parser(
        python_package='electronicparsers.atk',
        mainfile_name_re=r'^.*\.nc',
        mainfile_mime_re=r'application/octet-stream'),
    'parsers/bigdft': Parser(
        python_package='electronicparsers.bigdft',
        mainfile_contents_re=r'\|_____\|__:__\|__:__\|_____\|_____\|___ BBBBB          i     g         g\s*'),
    'parsers/castep': Parser(
        python_package='electronicparsers.castep',
        mainfile_contents_re=(r'\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*')),
    'parsers/charmm': Parser(
        python_package='electronicparsers.charmm',
        mainfile_contents_re=r'\s*Chemistry\s*at\s*HARvard\s*Macromolecular\s*Mechanics\s*',
        mainfile_mime_re=r'text/.*'),
    'parsers/cp2k': Parser(
        python_package='electronicparsers.cp2k',
        mainfile_contents_re=(
            r'\*\*\*\* \*\*\*\* \*\*\*\*\*\*  \*\*  PROGRAM STARTED AT\s.*\n'
            r' \*\*\*\*\* \*\* \*\*\*  \*\*\* \*\*   PROGRAM STARTED ON\s*.*\n'
            r' \*\*    \*\*\*\*   \*\*\*\*\*\*    PROGRAM STARTED BY .*\n'
            r' \*\*\*\*\* \*\*    \*\* \*\* \*\*   PROGRAM PROCESS ID .*\n'
            r'  \*\*\*\* \*\*  \*\*\*\*\*\*\*  \*\*  PROGRAM STARTED IN .*\n')),
    'parsers/cpmd': Parser(
        python_package='electronicparsers.cpmd',
        mainfile_contents_re=(r'\*\*\*       \*\*   \*\*\*  \*\* \*\*\*\* \*\*  \*\*   \*\*\*')),
    'parsers/crystal': Parser(
        python_package='electronicparsers.crystal',
        mainfile_contents_re=(
            fr'(\r?\n \*\s+CRYSTAL[\d]+\s+\*\r?\n \*\s*[a-zA-Z]+ : \d+[\.\d+]*)')),
    'parsers/dmol3': Parser(
        python_package='electronicparsers.dmol3',
        mainfile_name_re=r'.*\.outmol',
        mainfile_contents_re=r'Materials Studio DMol\^3'),
    'parsers/elk': Parser(
        python_package='electronicparsers.elk',
        mainfile_contents_re=r'\| Elk version [0-9.a-zA-Z]+ started \|'),
    'parsers/exciting': Parser(
        python_package='electronicparsers.exciting',
        mainfile_name_re=r'^.*.OUT(\.[^/]*)?$',
        mainfile_contents_re=(r'EXCITING.*started[\s\S]+?All units are atomic ')),
    'parsers/fhi-aims': Parser(
        python_package='electronicparsers.fhiaims',
        mainfile_contents_re=(
            r'^(.*\n)*'
            r'?\s*Invoking FHI-aims \.\.\.')),
    'parsers/fleur': Parser(
        python_package='electronicparsers.fleur',
        mainfile_contents_re=r'This output is generated by fleur.|\<fleurOutput',
        mainfile_name_re=r'.*[^/]*\.xml[^/]*',  # only the alternative mainfile name should match
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_alternative=True),
    'parsers/fplo': Parser(
        python_package='electronicparsers.fplo',
        mainfile_contents_re=r'\s*\|\s*FULL-POTENTIAL LOCAL-ORBITAL MINIMUM BASIS BANDSTRUCTURE CODE\s*\|\s*',
        mainfile_mime_re=r'text/.*'),
    'parsers/gamess': Parser(
        python_package='electronicparsers.gamess',
        mainfile_contents_re=(
            r'\s*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s*\*\s*GAMESS VERSION =\s*(.*)\*\s*'
            r'\s*\*\s*FROM IOWA STATE UNIVERSITY\s*\*\s*')),
    'parsers/gaussian': Parser(
        python_package='electronicparsers.gaussian',
        mainfile_mime_re=r'.*',
        mainfile_contents_re=(
            r'\s*Cite this work as:'
            r'\s*Gaussian [0-9]+, Revision [A-Za-z0-9\.]*,')),
    'parsers/gpaw': Parser(
        python_package='electronicparsers.gpaw',
        mainfile_name_re=(r'^.*\.(gpw2|gpw)$'),
        mainfile_mime_re=r'application/(x-tar|octet-stream)'),
    'parsers/molcas': Parser(
        python_package='electronicparsers.molcas',
        mainfile_contents_re=r'M O L C A S'),
    'parsers/mopac': Parser(
        python_package='electronicparsers.mopac',
        mainfile_contents_re=r'\s*\*\*\s*MOPAC\s*([0-9a-zA-Z\.]*)\s*\*\*\s*',
        mainfile_mime_re=r'text/.*',),
    'parsers/nwchem': Parser(
        python_package='electronicparsers.nwchem',
        mainfile_contents_re=(
            r'Northwest Computational Chemistry Package \(NWChem\) (\d+\.)+\d+')),
    'parsers/octopus': Parser(
        python_package='electronicparsers.octopus',
        mainfile_contents_re=(r'\|0\) ~ \(0\) \|')),
    'parsers/onetep': Parser(
        python_package='electronicparsers.onetep',
        mainfile_contents_re=r'####### #     # ####### ####### ####### ######'),
    'parsers/openmx': Parser(
        python_package='electronicparsers.openmx',
        mainfile_mime_re=r'(text/.*)',
        mainfile_name_re=r'.*\.out$',
        mainfile_contents_re=(r'^\*{59}\s+\*{59}\s+This calculation was performed by OpenMX'),),
    'parsers/orca': Parser(
        python_package='electronicparsers.orca',
        mainfile_contents_re=(
            r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s+\* O   R   C   A \*\s*'
            r'\s+\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\**\s*'
            r'\s*')),
    'parsers/psi4': Parser(
        python_package='electronicparsers.psi4',
        mainfile_contents_re=(r'Psi4: An Open-Source Ab Initio Electronic Structure Package')),
    'parsers/qball': Parser(
        python_package='electronicparsers.qball',
        mainfile_contents_re='qball',
        supported_compressions=["gz", "bz2", "xz"]),
    'parsers/qbox': Parser(
        python_package='electronicparsers.qbox',
        mainfile_mime_re=r'(application/xml)|(text/.*)',
        mainfile_contents_re=(r'http://qboxcode.org')),
    'parsers/quantumespresso': Parser(
        python_package='electronicparsers.quantumespresso',
        mainfile_contents_re=(r'(Program PWSCF.*starts)|(Current dimensions of program PWSCF are)'),
        supported_compressions=["gz", "bz2", "xz"]),
    'parsers/siesta': Parser(
        python_package='electronicparsers.siesta',
        mainfile_contents_re=(
            r'(Siesta Version: siesta-|SIESTA [0-9]\.[0-9]\.[0-9])|'
            r'(\*\s*WELCOME TO SIESTA\s*\*)')),
    'parsers/turbomole': Parser(
        python_package='electronicparsers.turbomole',
        mainfile_contents_re=(r'Copyright \(C\) [0-9]+ TURBOMOLE GmbH, Karlsruhe')),
    'parsers/vasp': Parser(
        python_package='electronicparsers.vasp',
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_name_re=r'.*[^/]*xml[^/]*',  # only the alternative mainfile name should match
        mainfile_contents_re=(
            r'^\s*<\?xml version="1\.0" encoding="ISO-8859-1"\?>\s*'
            r'?\s*<modeling>'
            r'?\s*<generator>'
            r'?\s*<i name="program" type="string">\s*vasp\s*</i>'
            r'?|^\svasp[\.\d]+.+?(?:\(build|complex)[\s\S]+?executed on'),
        supported_compressions=['gz', 'bz2', 'xz'], mainfile_alternative=True),
    'parsers/wien2k': Parser(
        python_package='electronicparsers.wien2k',
        mainfile_name_re=r'.*\.scf$',
        mainfile_alternative=True,
        mainfile_contents_re=r'\s*---------\s*:ITE[0-9]+:\s*[0-9]+\.\s*ITERATION\s*---------'),
    'parsers/yambo': Parser(
        python_package='electronicparsers.yambo',
        mainfile_contents_re=(r'Build.+\s+http://www\.yambo-code\.org')),
    'parsers/abacus': Parser(
        python_package='electronicparsers.abacus',
        mainfile_contents_re=(r'\s*\n\s*WELCOME TO ABACUS')),
    'parsers/amber': Parser(
        python_package='atomisticparsers.amber',
        mainfile_contents_re=r'\s*Amber\s[0-9]+\s[A-Z]+\s*[0-9]+'),
    'parsers/asap': Parser(
        python_package='atomisticparsers.asap',
        mainfile_name_re=r'.*.traj$', mainfile_mime_re=r'application/octet-stream',
        mainfile_binary_header_re=br'AFFormatASE\-Trajectory'),
    'parsers/bopfox': Parser(
        python_package='atomisticparsers.bopfox',
        mainfile_contents_re=r'\-+\s+BOPfox \(v'),
    'parsers/dftbplus': Parser(
        python_package='atomisticparsers.dftbplus',
        mainfile_contents_re=r'\|  DFTB\+',
        mainfile_mime_re=r'text/.*'),
    'parsers/dlpoly': Parser(
        python_package='atomisticparsers.dlpoly',
        mainfile_contents_re=(r'\*\*\s+DL_POLY.+\*\*'),),
    'parsers/gromacs': Parser(
        python_package='atomisticparsers.gromacs',
        mainfile_contents_re=r'gmx mdrun, (VERSION|version)[\s\S]*Input Parameters:'),
    'parsers/gromos': Parser(
        python_package='atomisticparsers.gromos',
        mainfile_contents_re=r'Bugreports to http://www.gromos.net'),
    'parsers/gulp': Parser(
        python_package='atomisticparsers.gulp',
        mainfile_contents_re=(
            r'\s*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*'
            r'\*\*\*\*\*\*\*\*\*\*\*\*\*\s*'
            r'\s*\*\s*GENERAL UTILITY LATTICE PROGRAM\s*\*\s*'),),
    'parsers/lammps': Parser(
        python_package='atomisticparsers.lammps',
        mainfile_contents_re=r'^LAMMPS\s+\(.+\)'),
    'parsers/libatoms': Parser(
        python_package='atomisticparsers.libatoms',
        mainfile_contents_re=(r'\s*<GAP_params\s'),),
    'parsers/namd': Parser(
        python_package='atomisticparsers.namd',
        mainfile_contents_re=r'\s*Info:\s*NAMD\s*[0-9.]+\s*for\s*',
        mainfile_mime_re=r'text/.*',),
    'parsers/tinker': Parser(
        python_package='atomisticparsers.tinker',
        mainfile_contents_re=r'TINKER  ---  Software Tools for Molecular Design'),
    'parsers/xtb': Parser(
        python_package='atomisticparsers.xtb',
        mainfile_contents_re=r'x T B\s+\|\s+\|\s+='),
    'parsers/aflow': Parser(
        python_package='workflowparsers.aflow',
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
        supported_compressions=['gz', 'bz2', 'xz'],
        mainfile_alternative=True),
    'parsers/asr': Parser(
        python_package='workflowparsers.asr',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=r'.*archive_.*\.json',
        mainfile_contents_re=(r'"name": "ASR"')),
    'parsers/elastic': Parser(
        python_package='workflowparsers.elastic',
        mainfile_contents_re=r'\s*Order of elastic constants\s*=\s*[0-9]+\s*',
        mainfile_name_re=(r'.*/INFO_ElaStic')),
    'parsers/fhivibes': Parser(
        python_package='workflowparsers.fhivibes',
        mainfile_name_re=(r'^.*\.(nc)$'),
        mainfile_mime_re=r'(application/x-hdf)',
        mainfile_binary_header_re=br'^\x89HDF',
        mainfile_contents_dict={'__has_all_keys': ['I', 'a', 'b']}),
    'parsers/lobster': Parser(
        python_package='workflowparsers.lobster',
        mainfile_name_re=r'.*lobsterout.*',
        mainfile_contents_re=(r'^LOBSTER\s*v[\d\.]+.*'),
        supported_compressions=['gz', 'bz2', 'xz'],),
    'parsers/atomate': Parser(
        python_package='workflowparsers.atomate',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_name_re=r'.*mp.+materials\.json',
        mainfile_contents_re=(r'"pymatgen_version":')),
    'parsers/phonopy': Parser(
        python_package='workflowparsers.phonopy',
        mainfile_name_re=(r'(.*/phonopy-FHI-aims-displacement-0*1/control.in$)|(.*/phon[^/]+yaml)')),
    'parsers/eelsdbparser': Parser(
        python_package='eelsdbparser',
        name='parsers/eelsdb',
        mainfile_mime_re=r'application/json',
        mainfile_contents_re=(r'https://eelsdb.eu/spectra')),
    'parsers/quantum_espresso_phonon': Parser(
        python_package='workflowparsers.quantum_espresso_phonon',
        mainfile_contents_re=(
            r'Program PHONON.+\s*'
            r'This program is part of the open-source Quantum ESPRESSO suite')),
    'parsers/quantum_espresso_epw': Parser(
        python_package='workflowparsers.quantum_espresso_epw',
        mainfile_contents_re=(
            r'Program EPW.+\s*'
            r'This program is part of the open-source Quantum ESPRESSO suite')),
    'parsers/quantum_espresso_xspectra': Parser(
        python_package='workflowparsers.quantum_espresso_xspectra',
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_contents_re=r'\s*Program XSpectra\s*'),
    'parsers/openkim': Parser(
        python_package='databaseparsers.openkim',
        mainfile_mime_re=r'(application/json)|(text/.*)',
        mainfile_contents_re=r'openkim|OPENKIM|OpenKIM'),
    'parsers/wannier90': Parser(
        python_package='electronicparsers.wannier90',
        mainfile_contents_re=r'\|\s*WANNIER90\s*\|'),
    'parsers/w2dynamics': Parser(
        python_package='electronicparsers.w2dynamics',
        mainfile_name_re=(r'^.*\.(h5|hdf5)$'),
        mainfile_mime_re=r'(application/x-hdf)',
        mainfile_binary_header_re=br'^\x89HDF',
        mainfile_contents_dict={'__has_all_keys': ['.axes', '.config', '.quantities']}),
    'parsers/soliddmft': Parser(
        python_package='electronicparsers.soliddmft',
        mainfile_name_re=(r'^.*\.(h5|hdf5)$'),
        mainfile_mime_re=r'(application/x-hdf)',
        mainfile_binary_header_re=br'^\x89HDF',
        mainfile_contents_dict={'__has_all_keys': ['dft_input', 'DMFT_input', 'DMFT_results']}),
    'parsers/ocean': Parser(
        python_package='electronicparsers.ocean',
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_contents_dict={'__has_all_keys': ['bse', 'structure', 'screen', 'calc']}),
    'parsers/edmft': Parser(
        python_package='electronicparsers.edmft',
        mainfile_name_re=(r'^.*\.(out)$'),
        mainfile_contents_re=r'\-\-\-\s*Preparing GF calculation\s*\-\-\-'),
    'parsers/nexus': Parser(
        python_package='nomad.parsing.nexus',
        mainfile_mime_re=r'(application/.*)|(text/.*)',
        mainfile_name_re=r'.*\.nxs',
        supported_compressions=['gz', 'bz2', 'xz']),
    'parsers/elabftw/elabftw': Parser(
        python_package='nomad.parsing.elabftw',
        parser_as_interface=True),
    'parsers/chemotion/chemotion': Parser(
        python_package='nomad.parsing.chemotion',
        parser_as_interface=True)
})


def normalize_loglevel(value, default_level=logging.INFO):
    plain_value = value
    if plain_value is None:
        return default_level
    else:
        try:
            return int(plain_value)
        except ValueError:
            return getattr(logging, plain_value)


_transformations = {
    'console_log_level': normalize_loglevel,
    'logstash_level': normalize_loglevel
}


# use std python logger, since logging is not configured while loading configuration
logger = logging.getLogger(__name__)


def _check_config():
    '''Used to check that the current configuration is valid. Should only be
    called once after the final config is loaded.

    Raises:
        AssertionError: if there is a contradiction or invalid values in the
            config file settings.
    '''
    # TODO more if this should be translated into pydantic validations.

    # The AFLOW symmetry information is checked once on import
    proto_symmetry_tolerance = normalize.prototype_symmetry_tolerance
    symmetry_tolerance = normalize.symmetry_tolerance
    if proto_symmetry_tolerance != symmetry_tolerance:
        raise AssertionError(
            "The AFLOW prototype information is outdated due to changed tolerance "
            "for symmetry detection. Please update the AFLOW prototype information "
            "by running the CLI command 'nomad admin ops prototype-update "
            "--matches-only'."
        )

    if normalize.springer_db_path and not os.path.exists(normalize.springer_db_path):
        normalize.springer_db_path = None

    if keycloak.public_server_url is None:
        keycloak.public_server_url = keycloak.server_url

    def get_external_path(path):
        if fs.external_working_directory and not os.path.isabs(path):
            return os.path.join(fs.external_working_directory, path)
        return path

    if fs.staging_external is None:
        fs.staging_external = get_external_path(fs.staging)

    if fs.public_external is None:
        fs.public_external = get_external_path(fs.public)

    if fs.north_home_external is None:
        fs.north_home_external = get_external_path(fs.north_home)

    ui.north.enabled = north.enabled
    ui.app_base = f'{"https" if services.https else "http"}://{services.api_host}:{services.api_port}{services.api_base_path.rstrip("/")}'
    ui.north_base = f'{"https" if services.https else "http"}://{north.hub_host}:{north.hub_port}{services.api_base_path.rstrip("/")}/north'


def _merge(a: Union[dict, BaseModel], b: Union[dict, BaseModel], path: List[str] = None) -> Union[dict, BaseModel]:
    '''
    Recursively merges b into a. Will add new key-value pairs, and will
    overwrite existing key-value pairs. Notice that this mutates the original
    dictionary/model a and if you want to return a copy, you will want to first
    (deep)copy the original value.
    '''

    def has(target, key):
        return key in target if isinstance(target, dict) else hasattr(target, key)

    def set(target, key, value):
        if isinstance(target, dict):
            target[key] = value
        else:
            setattr(target, key, value)

    def get(target, key):
        return target[key] if isinstance(target, dict) else getattr(target, key)

    if path is None: path = []
    for key in b.__dict__ if isinstance(b, BaseModel) else b:
        value = get(b, key)
        if has(a, key):
            child = get(a, key)
            if isinstance(value, (BaseModel, dict)) and isinstance(child, (BaseModel, dict)):
                _merge(child, value, path + [str(key)])
            else:
                set(a, key, value)
        else:
            set(a, key, value)
    return a


def _apply(key, value, raise_error: bool = True) -> None:
    '''
    Changes the config according to given key and value. The first part of a key
    (with ``_`` as a separator) is interpreted as a group of settings. E.g. ``fs_staging``
    leading to ``config.fs.staging``.
    '''
    full_key = key
    try:
        group_key, config_key = full_key.split('_', 1)
    except Exception:
        if raise_error:
            logger.error(f'config key does not exist: {full_key}')
        return

    current = globals()

    current_value: Any = None
    if group_key not in current:
        if key not in current:
            if raise_error:
                logger.error(f'config key does not exist: {full_key}')
            return
    else:
        current = current[group_key]
        if not isinstance(current, NomadSettings):
            if raise_error:
                logger.error(f'config key does not exist: {full_key}')
            return

        try:
            current_value = getattr(current, config_key)
        except AttributeError:
            if raise_error:
                logger.error(f'config key does not exist: {full_key}')
            return

        try:
            type_hints = get_type_hints(current.__class__)
            current_type = type_hints[config_key] if config_key in type_hints else type(current_value)
        except AttributeError:
            if raise_error:
                logger.error(f'config key does not exist: {full_key}')
            return

        key = config_key

    try:
        if current_value is not None:
            def default_transformation(value):
                return parse_obj_as(current_type, value)

            transformation = _transformations.get(full_key, default_transformation)
            value = transformation(value)

        if isinstance(value, (dict, BaseModel)):
            value = _merge(current_value, value)

        setattr(current, key, value)
        logger.info(f'set config setting {full_key}={value}')
    except Exception as e:
        logger.error(f'cannot set config setting {full_key}={value}: {e}')


def _apply_env_variables():
    kwargs = {
        key[len('NOMAD_'):].lower(): value
        for key, value in os.environ.items()
        if key.startswith('NOMAD_') and key != 'NOMAD_CONFIG'}

    for key, value in kwargs.items():
        _apply(key, value, raise_error=False)


def _apply_nomad_yaml():
    config_file = os.environ.get('NOMAD_CONFIG', 'nomad.yaml')

    if not os.path.exists(config_file):
        return

    with open(config_file, 'r') as stream:
        try:
            config_data = yaml.load(stream, Loader=getattr(yaml, 'FullLoader'))
        except yaml.YAMLError as e:
            logger.error(f'cannot read nomad config: {e}')
            return

    if not config_data:
        return

    for key, value in config_data.items():
        if isinstance(value, dict):
            group_key = key
            for key, value in value.items():
                _apply(f'{group_key}_{key}', value)
        else:
            _apply(key, value)


def load_config():
    '''
    Loads the configuration from nomad.yaml and environment.
    '''
    _apply_nomad_yaml()
    _apply_env_variables()
    _check_config()


load_config()
