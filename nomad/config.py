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
import inspect
import os.path
import yaml
import warnings
from typing import TypeVar, List, Any, cast
from pkg_resources import get_distribution, DistributionNotFound
from pydantic import BaseModel, Field

try:
    __version__ = get_distribution("nomad-lab").version
except DistributionNotFound:
    # package is not installed
    pass


warnings.filterwarnings('ignore', message='numpy.dtype size changed')
warnings.filterwarnings('ignore', message='numpy.ufunc size changed')
warnings.filterwarnings('ignore', category=DeprecationWarning)


NomadSettingsBound = TypeVar('NomadSettingsBound', bound='NomadSettings')


class NomadSettings(BaseModel):
    def customize(self: NomadSettingsBound, custom_settings: NomadSettingsBound, **kwargs) -> NomadSettingsBound:
        '''
        Returns a new config object, created by taking a copy of the current config and
        updating it with the settings defined in `custom_settings`. The `custom_settings` dict
        must not contain any new keys (keys not defined in this NomadSettings). If it does,
        an exception will be raised.
        '''

        rv = self.copy(deep=True)

        if custom_settings:
            for field_name in custom_settings.__fields__.keys():
                try:
                    setattr(rv, field_name, getattr(custom_settings, field_name))
                except Exception:
                    raise AssertionError(f'Invalid setting: {field_name}')

        for key, value in kwargs.items():
            if value is None:
                continue
            try:
                setattr(rv, key, value)
            except Exception:
                raise AssertionError(f'Invalid setting: {field_name}')

        return cast(NomadSettingsBound, rv)


class Services(NomadSettings):
    '''
    Contains basic configuration of the NOMAD services (app, worker, north).
    '''
    api_host = Field('localhost', description='''
        The external hostname that clients can use to reach this NOMAD installation.
    ''')
    api_port = Field(8000, description='''
        The port used to expose the NOMAD app and api to clients.
    ''')
    api_base_path = Field('/fairdi/nomad/latest', description='''
        The base path prefix for the NOMAD app and api.
    ''')
    api_secret = Field('defaultApiSecret', description='''
        A secret that is used to issue download and other tokens.
    ''')
    https = Field(False, description='''
        Set to `True`, if external clients are using *SSL* to connect to this installation.
        Requires to setup a reverse-proxy (e.g. the one used in the docker-compose
        based installation) that handles the *SSL* encryption.
    ''')
    https_upload = Field(False, description='''
        Set to `True`, if upload curl commands should suggest the use of SSL for file
        uploads. This can be configured independently of `https` to suggest large file
        via regular HTTP.
    ''')
    admin_user_id = Field('00000000-0000-0000-0000-000000000000', description='''
        The admin user `user_id`. All users are treated the same; there are no
        particular authorization information attached to user accounts. However, the
        API will grant the user with the given `user_id` more rights, e.g. using the
        `admin` owner setting in accessing data.
    ''')

    encyclopedia_base = Field(
        'https://nomad-lab.eu/prod/rae/encyclopedia/#', description='''
            This enables links to the given *encyclopedia* installation in the UI.
        ''')
    aitoolkit_enabled = Field(False, description='''
        If true, the UI will show a menu with links to the AI Toolkit notebooks on
        `nomad-lab.eu`.
    ''')

    console_log_level = Field(logging.WARNING, description='''
        The log level that controls console logging for all NOMAD services (app, worker, north).
        The level is given in Python `logging` log level numbers.
    ''')

    upload_limit = Field(10, description='''
        The maximum allowed unpublished uploads per user. If a user exceeds this
        amount, the user cannot add more uploads.
    ''')
    force_raw_file_decoding = Field(False, description='''
        By default, text raw-files are interpreted with utf-8 encoding. If this fails,
        the actual encoding is guessed. With this setting, we force to assume iso-8859-1
        encoding, if a file is not decodable with utf-8.
    ''')
    max_entry_download = Field(50000, description='''
        There is an inherent limit in page-based pagination with Elasticsearch. If you
        increased this limit with your Elasticsearch, you can also adopt this setting
        accordingly, changing the maximum amount of entries that can be paginated with
        page-base pagination.

        Page-after-value-based pagination is independent and can be used without limitations.
    ''')
    unavailable_value = Field('unavailable', description='''
        Value that is used in `results` section Enum fields (e.g. system type, spacegroup, etc.)
        to indicate that the value could not be determined.
    ''')


services = Services()


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


class Meta(NomadSettings):
    '''
    Metadata about the deployment and how it is presented to clients.
    '''
    version = Field(__version__, description='The NOMAD version string.')
    commit = Field('', description='The source-code commit that this installation\'s NOMAD version is build from.')
    deployment = Field(
        'devel', description='Human-friendly name of this nomad deployment.')
    deployment_url = Field(
        api_url(), description='The NOMAD deployment\'s url (api url).')
    label: str = Field(None, description='''
        An additional log-stash data key-value pair added to all logs. Can be used
        to differentiate deployments when analyzing logs.
    ''')
    service = Field('unknown nomad service', description='''
        Name for the service that is added to all logs. Depending on how NOMAD is
        installed, services get a name (app, worker, north) automatically.
    ''')

    name = Field(
        'NOMAD',
        description='Web-site title for the NOMAD UI.',
        deprecated=True)
    homepage = Field(
        'https://nomad-lab.eu', description='Provider homepage.', deprecated=True)
    source_url = Field(
        'https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR',
        description='URL of the NOMAD source-code repository.',
        deprecated=True)

    maintainer_email = Field(
        'markus.scheidgen@physik.hu-berlin.de',
        description='Email of the NOMAD deployment maintainer.')
    beta: dict = Field(None, description='''
        Additional data that describes how the deployment is labeled as a beta-version in the UI.
    ''')


meta = Meta()


class Oasis(NomadSettings):
    '''
    Settings related to the configuration of a NOMAD Oasis deployment.
    '''
    is_oasis = Field(False, description='Set to `True` to indicate that this deployment is a NOMAD Oasis.')
    allowed_users: str = Field(None, description='''
        A list of usernames or user account emails. These represent a white-list of
        allowed users. With this, users will need to login right-away and only the
        listed users might use this deployment. All API requests must have authentication
        information as well.''')
    uses_central_user_management = Field(False, description='''
        Set to True to use the central user-management. Typically the NOMAD backend is
        using the configured `keycloak` to access user data. With this, the backend will
        use the API of the central NOMAD (`central_nomad_deployment_url`) instead.
    ''')
    central_nomad_deployment_url = Field('https://nomad-lab.eu/prod/v1/api', description='''
        The URL of the API of the NOMAD deployment that is considered the *central* NOMAD.
    ''')


oasis = Oasis()


_jupyterhub_config_description = '''
    This setting is forwarded to jupyterhub; refer to the jupyterhub
    [documentation](https://jupyterhub.readthedocs.io/en/stable/api/app.html#).
'''


class North(NomadSettings):
    '''
    Settings related to the operation of the NOMAD remote tools hub service *north*.
    '''
    hub_connect_ip: str = Field(None, description='''
        Overwrites the default hostname that can be used from within a north container
        to reach the host system.

        Typically has to be set for non Linux hosts. Set this to `host.docker.internal`
        on windows/macos.
    ''')
    hub_connect_url: str = Field(None, description=_jupyterhub_config_description)
    hub_ip = Field('0.0.0.0', description=_jupyterhub_config_description)
    docker_network: str = Field(None, description=_jupyterhub_config_description)
    hub_host = Field('localhost', description='''
        The internal host name that NOMAD services use to connect to the jupyterhub API.
    ''')
    hub_port = Field(9000, description='''
        The internal port that NOMAD services use to connect to the jupyterhub API.
    ''')
    jupyterhub_crypt_key: str = Field(None, description=_jupyterhub_config_description)

    shared_fs = Field('.volumes/fs/north/shared', description='''
        Path to the shared folder on the host machine. This is mounted into spawned
        containers to be shared by all users.
    ''')
    users_fs = Field('.volumes/fs/north/users', description='''
        Path to a folder on the host machine. Sub-directories with the username are mounted
        into spawned containers to persist files per user.
    ''')
    nomad_host: str = Field(
        None, description='The NOMAD app host name that spawned containers use.')
    windows = Field(
        True, description='Enable windows OS hacks.')


north = North()


CELERY_WORKER_ROUTING = 'worker'
CELERY_QUEUE_ROUTING = 'queue'


class RabbitMQ(NomadSettings):
    '''
    Configures how NOMAD is connecting to RabbitMQ.
    '''
    host = Field('localhost', description='The name of the host that runs RabbitMQ.')
    user = Field('rabbitmq', description='The RabbitMQ user that is used to connect.')
    password = Field('rabbitmq', description='The password that is used to connect.')


rabbitmq = RabbitMQ()


def rabbitmq_url():
    return 'pyamqp://%s:%s@%s//' % (rabbitmq.user, rabbitmq.password, rabbitmq.host)


class Celery(NomadSettings):
    max_memory = 64e6  # 64 GB
    timeout = 1800  # 1/2 h
    acks_late = False
    routing = CELERY_QUEUE_ROUTING
    priorities = {
        'Upload.process_upload': 5,
        'Upload.delete_upload': 9,
        'Upload.publish_upload': 10
    }


celery = Celery()


class FS(NomadSettings):
    tmp = '.volumes/fs/tmp'
    staging = '.volumes/fs/staging'
    staging_external: str = None
    public = '.volumes/fs/public'
    public_external: str = None
    local_tmp = '/tmp'
    prefix_size = 2
    archive_version_suffix = 'v1'
    working_directory = os.getcwd()
    external_working_directory: str = None


fs = FS()


class Elastic(NomadSettings):
    host = 'localhost'
    port = 9200
    timeout = 60
    bulk_timeout = 600
    bulk_size = 1000
    entries_per_material_cap = 1000
    entries_index = 'nomad_entries_v1'
    materials_index = 'nomad_materials_v1'


elastic = Elastic()


class Keycloak(NomadSettings):
    server_url = 'https://nomad-lab.eu/fairdi/keycloak/auth/'
    public_server_url: str = None
    realm_name = 'fairdi_nomad_prod'
    username = 'admin'
    password = 'password'
    client_id = 'nomad_public'
    client_secret: str = None


keycloak = Keycloak()


class Mongo(NomadSettings):
    ''' Connection and usage settings for MongoDB.'''
    host: str = Field('localhost', description='The name of the host that runs mongodb.')
    port: int = Field(27017, description='The port to connect with mongodb.')
    db_name: str = Field('nomad_v1', description='The used mongodb database name.')


mongo = Mongo()


class Logstash(NomadSettings):
    enabled = False
    host = 'localhost'
    tcp_port = '5000'
    level: int = logging.DEBUG


logstash = Logstash()


class Tests(NomadSettings):
    default_timeout = 60


tests = Tests()


class Mail(NomadSettings):
    enabled = False
    with_login = False
    host = ''
    port = 8995
    user = ''
    password = ''
    from_address = 'support@nomad-lab.eu'
    cc_address = 'support@nomad-lab.eu'


mail = Mail()


class Normalize(NomadSettings):
    system_classification_with_clusters_threshold = Field(
        64, description='''
            The system size limit for running the dimensionality analysis. For very
            large systems the dimensionality analysis will get too expensive.
        ''')
    symmetry_tolerance = Field(
        0.1, description='''
            Symmetry tolerance controls the precision used by spglib in order to find
            symmetries. The atoms are allowed to move 1/2*symmetry_tolerance from
            their symmetry positions in order for spglib to still detect symmetries.
            The unit is angstroms. The value of 0.1 is used e.g. by Materials Project
            according to
            https://pymatgen.org/pymatgen.symmetry.analyzer.html#pymatgen.symmetry.analyzer.SpacegroupAnalyzer
        ''')
    prototype_symmetry_tolerance = Field(
        0.1, description='''
            The symmetry tolerance used in aflow prototype matching. Should only be
            changed before re-running the prototype detection.
        ''')
    max_2d_single_cell_size = Field(
        7, description='''
            Maximum number of atoms in the single cell of a 2D material for it to be
            considered valid.
        ''')
    cluster_threshold = Field(
        2.5, description='''
            The distance tolerance between atoms for grouping them into the same
            cluster. Used in detecting system type.
        ''')

    angle_rounding = Field(
        float(10.0), description='''
            Defines the "bin size" for rounding cell angles for the material hash in degree.
        ''')
    flat_dim_threshold = Field(
        0.1, description='''
            The threshold for a system to be considered "flat". Used e.g. when
            determining if a 2D structure is purely 2-dimensional to allow extra rigid
            transformations that are improper in 3D but proper in 2D.
        ''')

    k_space_precision = Field(
        150e6, description='''
            The threshold for point equality in k-space. Unit: 1/m.
        ''')
    band_structure_energy_tolerance = Field(
        8.01088e-21, description='''
            The energy threshold for how much a band can be on top or below the fermi
            level in order to still detect a gap. Unit: Joule.
        ''')
    springer_db_path = Field(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), 'normalizing/data/springer.msg'))


normalize = Normalize()


class Resources(NomadSettings):
    enabled = False
    db_name = 'nomad_v1_resources'
    max_time_in_mongo = Field(
        60 * 60 * 24 * 365., description='''
            Maxmimum time a resource is stored in mongodb before being updated.
        ''')
    download_retries = Field(
        2, description='Number of retries when downloading resources.')
    download_retry_delay = Field(
        10, description='Delay between retries in seconds')
    max_connections = Field(
        10, description='Maximum simultaneous connections used to download resources.')


resources = Resources()


class Client(NomadSettings):
    user: str = None
    password: str = None
    access_token: str = None
    url = 'http://nomad-lab.eu/prod/v1/api'


client = Client()


class DataCite(NomadSettings):
    mds_host = 'https://mds.datacite.org'
    enabled = False
    prefix = '10.17172'
    user = '*'
    password = '*'


datacite = DataCite()


class GitLab(NomadSettings):
    private_token = 'not set'


gitlab = GitLab()


class Reprocess(NomadSettings):
    '''
    Configures standard behaviour when reprocessing.
    Note, the settings only matter for published uploads and entries. For uploads in
    staging, we always reparse, add newfound entries, and delete unmatched entries.
    '''
    rematch_published = True
    reprocess_existing_entries = True
    use_original_parser = False
    add_matched_entries_to_published = True
    delete_unmatched_published_entries = False
    index_individual_entries = False


reprocess = Reprocess()


class Process(NomadSettings):
    store_package_definition_in_mongo = Field(
        False, description='Configures whether to store the corresponding package definition in mongodb.')
    add_definition_id_to_reference = Field(
        False, description='''
            Configures whether to attach definition id to `m_def`, note it is different from `m_def_id`.
            The `m_def_id` will be exported with the `with_def_id=True` via `m_to_dict`.
        ''')
    write_definition_id_to_archive = Field(False, description='Write `m_def_id` to the archive.')
    index_materials = True
    reuse_parser = True
    metadata_file_name = 'nomad'
    metadata_file_extensions = ('json', 'yaml', 'yml')
    auxfile_cutoff = 100
    parser_matching_size = 150 * 80  # 150 lines of 80 ASCII characters per line
    max_upload_size = 32 * (1024 ** 3)
    use_empty_parsers = False


process = Process()


class RFC3161Timestamp(NomadSettings):
    server = Field(
        'http://time.certum.pl/', description='The rfc3161ng timestamping host.')
    cert: str = Field(
        None, description='Path to the optional rfc3161ng timestamping server certificate.')
    hash_algorithm = Field(
        'sha256', description='Hash algorithm used by the rfc3161ng timestamping server.')
    username: str = None
    password: str = None


rfc3161_timestamp = RFC3161Timestamp()


class BundleExportSettings(NomadSettings):
    include_raw_files = True
    include_archive_files = True
    include_datasets = True


class BundleExport(NomadSettings):
    default_cli_bundle_export_path: str = './bundles'
    default_settings = BundleExportSettings()
    default_settings_cli: BundleExportSettings = Field(None, description='''
        Additional default settings, applied when exporting using the CLI (command-line interface).
        This allows to override some of the settings specified in the general default settings above.
    ''')


bundle_export = BundleExport()


class BundleImportSettings(NomadSettings):
    process_settings = Field(
        Reprocess(
            rematch_published=True,
            reprocess_existing_entries=True,
            use_original_parser=False,
            add_matched_entries_to_published=True,
            delete_unmatched_published_entries=False
        ), description='''
            It is possible to trigger processing of the raw files, but it is no longer the
            preferred way to import bundles. If used, the settings below control the reprocessing
            behaviour (see the config for `reprocess` for more info).
        '''
    )

    include_raw_files = True
    include_archive_files = False
    include_datasets = True

    include_bundle_info = Field(
        True, description='Keeps the bundle_info.json file, not necessary but nice to have.')
    keep_original_timestamps = Field(
        False, description='''
            If all timestamps (create time, publish time etc) should be imported from
            the bundle.
        ''')
    set_from_oasis = Field(
        True, description='If the from_oasis flag and oasis_deployment_url should be set.')

    delete_upload_on_fail = Field(
        False, description='If False, it is just removed from the ES index on failure.')
    delete_bundle_on_fail = Field(
        True, description='Deletes the source bundle if the import fails.')
    delete_bundle_on_success = Field(
        True, description='Deletes the source bundle if the import succeeds.')
    delete_bundle_include_parent_folder = Field(
        True, description='When deleting the bundle, also include parent folder, if empty.')

    trigger_processing = Field(
        True, description='If the upload should be processed when the import is done.')


class BundleImport(NomadSettings):

    required_nomad_version = Field(
        '1.1.2', description='Minimum  NOMAD version of bundles required for import.')

    default_cli_bundle_import_path = './bundles'

    allow_bundles_from_oasis = Field(
        False, description='If oasis admins can "push" bundles to this NOMAD deployment.')
    allow_unpublished_bundles_from_oasis = Field(
        False, description='If oasis admins can "push" bundles of unpublished uploads.')

    default_settings = BundleImportSettings()

    default_settings_cli = Field(
        BundleImportSettings(
            delete_bundle_on_fail=False,
            delete_bundle_on_success=False
        ),
        description='''
            Additional default settings, applied when importing using the CLI (command-line interface).
            This allows to override some of the settings specified in the general default settings above.
        '''
    )


bundle_import = BundleImport()


class Archive(NomadSettings):
    block_size = 256 * 1024
    read_buffer_size = Field(
        256 * 1024, description='GPFS needs at least 256K to achieve decent performance.')
    max_process_number = Field(
        20, description='Maximum number of processes can be assigned to process archive query.')
    min_entries_per_process = Field(
        20, description='Minimum number of entries per process.')


archive = Archive()


class UIConfig(NomadSettings):
    default_unit_system = 'Custom'
    entry_context: dict = {
        'overview': {
            'include': [
                'sections',
                'definitions',
                'nexus',
                'material',
                'electronic',
                'optoelectronic',
                'vibrational',
                'mechanical',
                'thermodynamic',
                'structural',
                'dynamical',
                'geometry_optimization',
                'spectroscopy',
                'references'
            ],
            'exclude': [
                'relatedResources'
            ],
            'options': {
                'sections': {'error': 'Could not render section card.'},
                'definitions': {'error': 'Could not render definitions card.'},
                'nexus': {'error': 'Could not render NeXus card.'},
                'material': {'error': 'Could not render material card.'},
                'electronic': {'error': 'Could not render electronic properties.'},
                'optoelectronic': {'error': 'Could not render optoelectronic properties.'},
                'vibrational': {'error': 'Could not render vibrational properties.'},
                'mechanical': {'error': 'Could not render mechanical properties.'},
                'thermodynamic': {'error': 'Could not render thermodynamic properties.'},
                'structural': {'error': 'Could not render structural properties.'},
                'geometry_optimization': {'error': 'Could not render geometry optimization.'},
                'spectroscopy': {'error': 'Could not render spectroscopic properties.'},
                'references': {'error': 'Could not render references card.'},
                'relatedResources': {'error': 'Could not render related resources card.'},
            }
        }
    }
    search_contexts: dict = {
        "include": ["entries", "eln", "materials", "solar_cells"],
        "exclude": [],
        "options": {
            "entries": {
                'label': "Entries",
                'path': "entries",
                'resource': 'entries',
                'breadcrumb': "Entries search",
                'description': "Search individual entries",
                'help': {
                    'title': 'Entries search',
                    'content': inspect.cleandoc(r'''
                        This page allows you to **search entries** within NOMAD. Entries represent
                        individual calculations or experiments that have bee uploaded into NOMAD.

                        The search page consists of three main elements: the filter panel, the search
                        bar, and the result list.

                        The filter panel on the left allows you to graphically explore and enter
                        different search filters. It also gives a visual indication of the currently
                        active search filters for each category. This is a good place to start exploring
                        the available search filters and their meaning.

                        The search bar allows you to specify filters by typing them in and pressing
                        enter. You can also start by simply typing keywords of interest, which will
                        toggle a list of suggestions. For numerical data you can also use range queries,
                        e.g. \`0.0 < band_gap <= 0.1\`.

                        Notice that the units used in the filter panel and in the queries can be changed
                        using the **units** button on the top right corner. When using the search bar,
                        you can also specify a unit by typing the unit abbreviations, e.g. \`band_gap >=
                        0.1 Ha\`

                        The result list on the right is automatically updated according to the filters
                        you have specified. You can browse through the results by scrolling through the
                        available items and loading more results as you go. Here you can also change the
                        sorting of the results, modify the displayed columns, access individual entries
                        or even download the raw data or the archive document by selecting individual
                        entries and pressing the download button that appears. The ellipsis button shown
                        for each entry will navigate you to that entry's page. This entry page will show
                        more metadata, raw files, the entry's archive, and processing logs.
                    '''),
                },
                'pagination': {
                    'order_by': 'upload_create_time',
                    'order': 'desc',
                    'page_size': 20,
                },
                'columns': {
                    'enable': [
                        'entry_name',
                        'results.material.chemical_formula_hill',
                        'entry_type',
                        'upload_create_time',
                        'authors'
                    ],
                    'include': [
                        'entry_name',
                        'results.material.chemical_formula_hill',
                        'entry_type',
                        'results.method.method_name',
                        'results.method.simulation.program_name',
                        'results.method.simulation.dft.basis_set_name',
                        'results.method.simulation.dft.xc_functional_type',
                        'results.material.structural_type',
                        'results.material.symmetry.crystal_system',
                        'results.material.symmetry.space_group_symbol',
                        'results.material.symmetry.space_group_number',
                        'results.eln.lab_ids',
                        'results.eln.sections',
                        'results.eln.methods',
                        'results.eln.tags',
                        'results.eln.instruments',
                        'mainfile',
                        'upload_create_time',
                        'authors',
                        'comment',
                        'references',
                        'datasets',
                        'published',
                    ],
                    'exclude': [],
                    'options': {
                        'entry_name': {'label': 'Name', 'align': 'left'},
                        'results.material.chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                        'entry_type': {'label': 'Entry type', 'align': 'left'},
                        'results.method.method_name': {'label': 'Method name'},
                        'results.method.simulation.program_name': {'label': 'Program name'},
                        'results.method.simulation.dft.basis_set_name': {'label': 'Basis set name'},
                        'results.method.simulation.dft.xc_functional_type': {'label': 'XC functional type'},
                        'results.material.structural_type': {'label': 'Structural type'},
                        'results.material.symmetry.crystal_system': {'label': 'Crystal system'},
                        'results.material.symmetry.space_group_symbol': {'label': 'Space group symbol'},
                        'results.material.symmetry.space_group_number': {'label': 'Space group number'},
                        'results.eln.lab_ids': {'label': 'Lab IDs'},
                        'results.eln.sections': {'label': 'Sections'},
                        'results.eln.methods': {'label': 'Methods'},
                        'results.eln.tags': {'label': 'Tags'},
                        'results.eln.instruments': {'label': 'Instruments'},
                        'mainfile': {'label': 'Mainfile', 'align': 'left'},
                        'upload_create_time': {'label': 'Upload time', 'align': 'left'},
                        'authors': {'label': 'Authors', 'align': 'left'},
                        'comment': {'label': 'Comment', 'align': 'left'},
                        'references': {'label': 'References', 'align': 'left'},
                        'datasets': {'label': 'Datasets', 'align': 'left'},
                        'published': {'label': 'Access'}
                    }
                },
                'rows': {
                    'actions': {
                        'enable': True
                    },
                    'details': {
                        'enable': True
                    },
                    'selection': {
                        'enable': True
                    }
                },
                'filter_menus': {
                    'include': [
                        'material',
                        'elements',
                        'symmetry',
                        'method',
                        'simulation',
                        'dft',
                        'gw',
                        'experiment',
                        'eels',
                        'properties',
                        'electronic',
                        'optoelectronic',
                        'vibrational',
                        'mechanical',
                        'spectroscopy',
                        'thermodynamic',
                        'geometry_optimization',
                        'eln',
                        'author',
                        'dataset',
                        'access',
                        'ids',
                        'processed_data_quantities',
                        'optimade',
                    ],
                    'exclude': [],
                    'options': {
                        'material': {'label': 'Material', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'elements': {'label': 'Elements / Formula', 'level': 1, 'size': 'large', 'menu_items': {}},
                        'symmetry': {'label': 'Symmetry', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'method': {'label': 'Method', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'simulation': {'label': 'Simulation', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'dft': {'label': 'DFT', 'level': 2, 'size': 'small', 'menu_items': {}},
                        'gw': {'label': 'GW', 'level': 2, 'size': 'small', 'menu_items': {}},
                        'experiment': {'label': 'Experiment', 'level': 1, 'size': 'small'},
                        'eels': {'label': 'EELS', 'level': 2, 'size': 'small', 'menu_items': {}},
                        'properties': {'label': 'Properties', 'level': 0, 'size': 'small'},
                        'electronic': {'label': 'Electronic', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'optoelectronic': {'label': 'Optoelectronic', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'vibrational': {'label': 'Vibrational', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'mechanical': {'label': 'Mechanical', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'spectroscopy': {'label': 'Spectroscopy', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'thermodynamic': {'label': 'Thermodynamic', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'geometry_optimization': {'label': 'Geometry Optimization', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'eln': {'label': 'Electronic Lab Notebook', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'author': {'label': 'Author / Origin', 'level': 0, 'size': 'medium', 'menu_items': {}},
                        'dataset': {'label': 'Dataset', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'access': {'label': 'Access', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'ids': {'label': 'IDs', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'processed_data_quantities': {'label': 'Processed Data Quantities', 'level': 0, 'size': 'medium', 'menu_items': {}},
                        'optimade': {'label': 'Optimade', 'level': 0, 'size': 'medium', 'menu_items': {}},
                    }
                }
            },
            "eln": {
                'label': "ELN",
                'path': "eln",
                'resource': 'entries',
                'breadcrumb': "ELN entries search",
                'description': "Search individual ELN entries",
                'help': {
                    'title': 'ELN entries search',
                    'content': inspect.cleandoc(r'''
                        This page allows you to specifically **search ELN entries** within NOMAD.
                        It is very similar to the *Entries search*, but with a reduced
                        filter set and specialized arrangement of default columns.
                    '''),
                },
                'pagination': {
                    'order_by': 'upload_create_time',
                    'order': 'desc',
                    'page_size': 20,
                },
                'columns': {
                    'enable': [
                        'entry_name',
                        'entry_type',
                        'upload_create_time',
                        'authors'
                    ],
                    'include': [
                        'entry_name',
                        'results.material.chemical_formula_hill',
                        'entry_type',
                        'results.method.method_name',
                        'results.method.simulation.program_name',
                        'results.method.simulation.dft.basis_set_name',
                        'results.method.simulation.dft.xc_functional_type',
                        'results.material.structural_type',
                        'results.material.symmetry.crystal_system',
                        'results.material.symmetry.space_group_symbol',
                        'results.material.symmetry.space_group_number',
                        'results.eln.lab_ids',
                        'results.eln.sections',
                        'results.eln.methods',
                        'results.eln.tags',
                        'results.eln.instruments',
                        'mainfile',
                        'upload_create_time',
                        'authors',
                        'comment',
                        'references',
                        'datasets',
                        'published',
                    ],
                    'exclude': [],
                    'options': {
                        'entry_name': {'label': 'Name', 'align': 'left'},
                        'results.material.chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                        'entry_type': {'label': 'Entry type', 'align': 'left'},
                        'results.method.method_name': {'label': 'Method name'},
                        'results.method.simulation.program_name': {'label': 'Program name'},
                        'results.method.simulation.dft.basis_set_name': {'label': 'Basis set name'},
                        'results.method.simulation.dft.xc_functional_type': {'label': 'XC functional type'},
                        'results.material.structural_type': {'label': 'Structural type'},
                        'results.material.symmetry.crystal_system': {'label': 'Crystal system'},
                        'results.material.symmetry.space_group_symbol': {'label': 'Space group symbol'},
                        'results.material.symmetry.space_group_number': {'label': 'Space group number'},
                        'results.eln.lab_ids': {'label': 'Lab IDs'},
                        'results.eln.sections': {'label': 'Sections'},
                        'results.eln.methods': {'label': 'Methods'},
                        'results.eln.tags': {'label': 'Tags'},
                        'results.eln.instruments': {'label': 'Instruments'},
                        'mainfile': {'label': 'Mainfile', 'align': 'left'},
                        'upload_create_time': {'label': 'Upload time', 'align': 'left'},
                        'authors': {'label': 'Authors', 'align': 'left'},
                        'comment': {'label': 'Comment', 'align': 'left'},
                        'references': {'label': 'References', 'align': 'left'},
                        'datasets': {'label': 'Datasets', 'align': 'left'},
                        'published': {'label': 'Access'}
                    }
                },
                'rows': {
                    'actions': {
                        'enable': True
                    },
                    'details': {
                        'enable': True
                    },
                    'selection': {
                        'enable': True
                    }
                },
                'filter_menus': {
                    'include': [
                        'material',
                        'elements',
                        'eln',
                        'custom_quantities',
                        'author',
                        'dataset',
                        'access',
                        'ids',
                        'processed_data_quantities',
                        'optimade',
                    ],
                    'exclude': [],
                    'options': {
                        'material': {'label': 'Material', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'elements': {'label': 'Elements / Formula', 'level': 1, 'size': 'large', 'menu_items': {}},
                        'eln': {'label': 'Electronic Lab Notebook', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'custom_quantities': {'label': 'Custom quantities', 'level': 0, 'size': 'large', 'menu_items': {}},
                        'author': {'label': 'Author / Origin', 'level': 0, 'size': 'medium', 'menu_items': {}},
                        'dataset': {'label': 'Dataset', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'access': {'label': 'Access', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'ids': {'label': 'IDs', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'processed_data_quantities': {'label': 'Processed Data Quantities', 'level': 0, 'size': 'medium', 'menu_items': {}},
                        'optimade': {'label': 'Optimade', 'level': 0, 'size': 'medium', 'menu_items': {}},
                    }
                }
            },
            "materials": {
                'label': "Materials",
                'path': "materials",
                'resource': 'materials',
                'breadcrumb': "Materials search",
                'description': "Search materials that are identified from the entries",
                'help': {
                    'title': 'Materials search',
                    'content': inspect.cleandoc(r'''
                        This page allows you to **search materials** within NOMAD. NOMAD can
                        automatically detect the material from individual entries and can then group the
                        data by using these detected materials. This allows you to search individual
                        materials which have properties that are aggregated from several entries.

                        The search page consists of three main elements: the filter panel, the search
                        bar, and the result list.

                        The filter panel on the left allows you to graphically explore and enter
                        different search filters. It also gives a visual indication of the currently
                        active search filters for each category. This is a good place to start exploring
                        the available search filters and their meaning.

                        The search bar allows you to specify filters by typing them in and pressing
                        enter. You can also start by simply typing keywords of interest, which will
                        toggle a list of suggestions. For numerical data you can also use range queries,
                        e.g. \`0.0 < band_gap <= 0.1\`.

                        The units used in the filter panel and in the queries can be changed
                        using the **units** button on the top right corner. When using the search bar,
                        you can also specify a unit by typing the unit abbreviations, e.g. \`band_gap >=
                        0.1 Ha\`.

                        Notice that by default the properties that you search can be combined from
                        several different entries. If instead you wish to search for a material with an
                        individual entry fullfilling your search criteria, uncheck the **combine results
                        from several entries**-checkbox.

                        The result list on the right is automatically updated according to the filters
                        you have specified. You can scroll through the available items and load more
                        results as you go. Here you can also change the sorting of the results, modify
                        the displayed columns and access individual materials. The ellipsis button shown
                        for each material will navigate you into the material overview page within the
                        NOMAD Encyclopedia. This page will show a more detailed overview for that
                        specific material.
                    '''),
                },
                'pagination': {
                    'order_by': 'chemical_formula_hill',
                    'order': 'asc'
                },
                'columns': {
                    'enable': [
                        'chemical_formula_hill',
                        'structural_type',
                        'symmetry.structure_name',
                        'symmetry.space_group_number',
                        'symmetry.crystal_system',
                    ],
                    'include': [
                        'chemical_formula_hill',
                        'structural_type',
                        'symmetry.structure_name',
                        'symmetry.crystal_system',
                        'symmetry.space_group_symbol',
                        'symmetry.space_group_number',
                        'material_id',
                    ],
                    'exclude': [],
                    'options': {
                        'chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                        'structural_type': {'label': 'Structural type'},
                        'symmetry.structure_name': {'label': 'Structure name'},
                        'symmetry.crystal_system': {'label': 'Crystal system'},
                        'symmetry.space_group_symbol': {'label': 'Space group symbol'},
                        'symmetry.space_group_number': {'label': 'Space group number'},
                        'material_id': {'label': 'Material ID'},
                    }
                },
                'rows': {
                    'actions': {
                        'enable': True
                    },
                    'details': {
                        'enable': False
                    },
                    'selection': {
                        'enable': False
                    }
                },
                'filter_menus': {
                    'include': [
                        'material',
                        'elements',
                        'symmetry',
                        'method',
                        'simulation',
                        'dft',
                        'gw',
                        'experiment',
                        'eels',
                        'properties',
                        'electronic',
                        'optoelectronic',
                        'vibrational',
                        'mechanical',
                        'spectroscopy',
                        'thermodynamic',
                        'geometry_optimization',
                        'eln',
                        'author',
                        'dataset',
                        'access',
                        'ids',
                        'processed_data_quantities',
                        'optimade',
                        'combine',
                    ],
                    'exclude': [],
                    'options': {
                        'material': {'label': 'Material', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'elements': {'label': 'Elements / Formula', 'level': 1, 'size': 'large', 'menu_items': {}},
                        'symmetry': {'label': 'Symmetry', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'method': {'label': 'Method', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'simulation': {'label': 'Simulation', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'dft': {'label': 'DFT', 'level': 2, 'size': 'small', 'menu_items': {}},
                        'gw': {'label': 'GW', 'level': 2, 'size': 'small', 'menu_items': {}},
                        'experiment': {'label': 'Experiment', 'level': 1, 'size': 'small'},
                        'eels': {'label': 'EELS', 'level': 2, 'size': 'small', 'menu_items': {}},
                        'properties': {'label': 'Properties', 'level': 0, 'size': 'small'},
                        'electronic': {'label': 'Electronic', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'optoelectronic': {'label': 'Optoelectronic', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'vibrational': {'label': 'Vibrational', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'mechanical': {'label': 'Mechanical', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'spectroscopy': {'label': 'Spectroscopy', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'thermodynamic': {'label': 'Thermodynamic', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'geometry_optimization': {'label': 'Geometry Optimization', 'level': 1, 'size': 'small', 'menu_items': {}},
                        'eln': {'label': 'Electronic Lab Notebook', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'author': {'label': 'Author / Origin', 'level': 0, 'size': 'medium', 'menu_items': {}},
                        'dataset': {'label': 'Dataset', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'access': {'label': 'Access', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'ids': {'label': 'IDs', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'processed_data_quantities': {'label': 'Processed Data Quantities', 'level': 0, 'size': 'medium', 'menu_items': {}},
                        'optimade': {'label': 'Optimade', 'level': 0, 'size': 'medium', 'menu_items': {}},
                        'combine': {
                            'actions': {
                                'include': ['combine'],
                                'options': {
                                    'combine': {
                                        'type': 'checkbox',
                                        'label': "Combine results from several entries",
                                        'quantity': 'combine'
                                    }
                                }
                            }
                        }
                    }
                }
            },
            "solar_cells": {
                'label': "Solar Cells",
                'path': "solar-cells",
                'resource': 'entries',
                'breadcrumb': "Solar cells search",
                'description': "Search solar cells in NOMAD",
                'help': {
                    'title': 'Solar cells search',
                    'content': inspect.cleandoc(r'''
                        This page allows you to **search solar cells** within NOMAD.
                        You can search for solar cells by their properties, by the chemistry of the absorber layer,
                        preparation method, or their origin.

                        The search page consists of three main elements: the filter panel, the search
                        bar, and the result list.

                        The filter panel on the left allows you to graphically explore and enter
                        different search filters. It also gives a visual indication of the currently
                        active search filters for each category. This is a good place to start exploring
                        the available search filters and their meaning. But clikcing in the "+" button
                        you could add the filter widgets to the central part of the page and combine them.
                        Try adding the periodic table widget to the main pannel and select some elements
                        contained in the absorber layer.

                        The search bar allows you to specify filters by typing them in and pressing
                        enter. You can also start by simply typing keywords of interest, which will
                        toggle a list of suggestions. For numerical data you can also use range queries,
                        e.g. \`0.0 < efficiency <= 20.1\`.

                        The units used in the filter panel and in the queries can be changed
                        using the **units** button on the top right corner. When using the search bar,
                        you can also specify a unit by typing the unit abbreviations, e.g. \`band_gap >=
                        0.1 Ha\`.

                        The result list on the right is automatically updated according to the filters
                        you have specified. You can browse through the results by scrolling through the
                        available items and loading more results as you go. Here you can also change the
                        sorting of the results, modify the displayed columns, download the raw data
                        or the archive document by selecting individual entries and pressing the download
                        cloud button that appears.

                        The ellipsis button (three dots) shown for each entry will navigate
                        you to that entry's page. This entry page will show more metadata,
                        raw files, the entry's archive, and processing logs.
                    '''),
                },
                'pagination': {
                    'order_by': 'results.properties.optoelectronic.solar_cell.efficiency',
                    'order': 'desc',
                    'page_size': 20,
                },
                'filters_locked': {
                    'results.material.functional_type': 'solar cell'
                },
                'columns': {
                    'enable': [
                        'results.material.chemical_formula_descriptive',
                        'results.properties.optoelectronic.solar_cell.efficiency',
                        'results.properties.optoelectronic.solar_cell.open_circuit_voltage',
                        'results.properties.optoelectronic.solar_cell.short_circuit_current_density',
                        'results.properties.optoelectronic.solar_cell.fill_factor',
                        'references'
                    ],
                    'include': [
                        'entry_name',
                        'results.material.chemical_formula_hill',
                        'results.material.chemical_formula_descriptive',
                        'results.properties.optoelectronic.solar_cell.efficiency',
                        'results.properties.optoelectronic.solar_cell.open_circuit_voltage',
                        'results.properties.optoelectronic.solar_cell.short_circuit_current_density',
                        'results.properties.optoelectronic.solar_cell.fill_factor',
                        # 'results.properties.optoelectronic.band_gap[0].value',
                        'results.properties.optoelectronic.solar_cell.device_stack',
                        'results.properties.optoelectronic.solar_cell.device_architecture',
                        'results.properties.optoelectronic.solar_cell.illumination_intensity',
                        'results.properties.optoelectronic.solar_cell.absorber_fabrication',
                        'entry_type',
                        'results.material.structural_type',
                        'results.eln.lab_ids',
                        'results.eln.sections',
                        'results.eln.methods',
                        'results.eln.tags',
                        'results.eln.instruments',
                        'mainfile',
                        'upload_create_time',
                        'authors',
                        'comment',
                        'references',
                        'datasets',
                        'published',
                    ],
                    'exclude': [],
                    'options': {
                        'entry_name': {'label': 'Name', 'align': 'left'},
                        'results.material.chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                        'results.material.chemical_formula_descriptive': {'label': 'Descriptive Formula', 'align': 'left'},
                        'entry_type': {'label': 'Entry type', 'align': 'left'},
                        'results.material.structural_type': {'label': 'Structural type'},
                        'results.properties.optoelectronic.solar_cell.efficiency': {
                            'label': 'Efficiency (%)',
                            'format': {
                                'decimals': 2,
                                'mode': 'standard',
                            },
                        },
                        'results.properties.optoelectronic.solar_cell.open_circuit_voltage': {
                            'label': 'Open circuit voltage',
                            'unit': 'V',
                            'format': {
                                'decimals': 3,
                                'mode': 'standard',
                            },
                        },
                        'results.properties.optoelectronic.solar_cell.short_circuit_current_density': {
                            'label': 'Short circuit current density',
                            'unit': 'A/m**2',
                            'format': {
                                'decimals': 3,
                                'mode': 'standard',
                            },
                        },
                        'results.properties.optoelectronic.solar_cell.fill_factor': {
                            'label': 'Fill factor',
                            'format': {
                                'decimals': 3,
                                'mode': 'standard',
                            },
                        },
                        'results.properties.optoelectronic.solar_cell.illumination_intensity': {
                            'label': 'Illum. intensity',
                            'unit': 'W/m**2',
                            'format': {
                                'decimals': 3,
                                'mode': 'standard',
                            },
                        },
                        # 'results.properties.optoelectronic.band_gap[0].value': {'label': 'Bandgap (eV)'},
                        'results.eln.lab_ids': {'label': 'Lab IDs'},
                        'results.eln.sections': {'label': 'Sections'},
                        'results.eln.methods': {'label': 'Methods'},
                        'results.eln.tags': {'label': 'Tags'},
                        'results.eln.instruments': {'label': 'Instruments'},
                        'mainfile': {'label': 'Mainfile', 'align': 'left'},
                        'upload_create_time': {'label': 'Upload time', 'align': 'left'},
                        'authors': {'label': 'Authors', 'align': 'left'},
                        'comment': {'label': 'Comment', 'align': 'left'},
                        'references': {'label': 'References', 'align': 'left'},
                        'datasets': {'label': 'Datasets', 'align': 'left'},
                        'published': {'label': 'Access'},
                    },
                },
                'rows': {
                    'actions': {
                        'enable': True
                    },
                    'details': {
                        'enable': True
                    },
                    'selection': {
                        'enable': True
                    }
                },
                'filter_menus': {
                    'include': [
                        'material',
                        'elements',
                        'properties',
                        'optoelectronic',
                        'eln',
                        'author',
                        'dataset',
                        'access',
                        'ids',
                        'processed_data_quantities',
                        'optimade',
                    ],
                    'exclude': [],
                    'options': {
                        'material': {'label': 'Material', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'elements': {'label': 'Elements / Formula', 'level': 0, 'size': 'large', 'menu_items': {}},
                        'properties': {'label': 'Properties', 'level': 0, 'size': 'small'},
                        'optoelectronic': {'label': 'Solar Cell Properties', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'eln': {'label': 'Electronic Lab Notebook', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'author': {'label': 'Author / Origin', 'level': 0, 'size': 'medium', 'menu_items': {}},
                        'dataset': {'label': 'Dataset', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'access': {'label': 'Access', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'ids': {'label': 'IDs', 'level': 0, 'size': 'small', 'menu_items': {}},
                        'processed_data_quantities': {'label': 'Processed Data Quantities', 'level': 0, 'size': 'medium', 'menu_items': {}},
                        'optimade': {'label': 'Optimade', 'level': 0, 'size': 'medium', 'menu_items': {}},
                    }
                }
            }
        }
    }


ui = UIConfig()


def north_url(ssl: bool = True):
    return api_url(ssl=ssl, api='north', api_host=north.hub_host, api_port=north.hub_port)


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
    """Used to check that the current configuration is valid. Should only be
    called once after the final config is loaded.

    Raises:
        AssertionError: if there is a contradiction or invalid values in the
            config file settings.
    """
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

    def set_external_path(source_obj, source_key, target_obj, target_key, overwrite=False):
        source_value = getattr(source_obj, source_key)
        target_value = getattr(target_obj, target_key)

        if target_value and not overwrite:
            return

        if not source_value:
            return

        if fs.external_working_directory and not os.path.isabs(source_value):
            target_value = os.path.join(fs.external_working_directory, source_value)
        else:
            target_value = source_value

        setattr(target_obj, target_key, target_value)

    set_external_path(fs, 'staging', fs, 'staging_external')
    set_external_path(fs, 'public', fs, 'public_external')
    set_external_path(north, 'users_fs', north, 'users_fs', overwrite=True)
    set_external_path(north, 'shared_fs', north, 'shared_fs', overwrite=True)


def _merge(a: dict, b: dict, path: List[str] = None) -> dict:
    '''
    Recursively merges b into a. Will add new key-value pairs, and will
    overwrite existing key-value pairs. Notice that this mutates the original
    dictionary a and if you want to return a copy, you will want to first
    (deep)copy the original dictionary.
    '''
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                _merge(a[key], b[key], path + [str(key)])
            else:
                a[key] = b[key]
        else:
            a[key] = b[key]
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

        key = config_key

    try:
        if current_value is not None and not isinstance(value, type(current_value)):
            value = _transformations.get(full_key, type(current_value))(value)

        if isinstance(value, dict):
            value = _merge(current_value, value)

        if isinstance(current, dict):
            current[key] = value
        else:
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
