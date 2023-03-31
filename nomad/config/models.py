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

import os
import json
from enum import Enum
import logging
import inspect
from typing import TypeVar, List, Dict, Any, Union, Optional, cast
from typing_extensions import Literal, Annotated  # type: ignore
from pydantic import BaseModel, Field, validator, Extra  # pylint: disable=unused-import
from pkg_resources import get_distribution, DistributionNotFound
try:
    __version__ = get_distribution('nomad-lab').version
except DistributionNotFound:
    # package is not installed
    pass


NomadSettingsBound = TypeVar('NomadSettingsBound', bound='NomadSettings')


class NomadSettings(BaseModel):
    def customize(
            self: NomadSettingsBound,
            custom_settings: Union[NomadSettingsBound, Dict[str, Any]]) -> NomadSettingsBound:
        '''
        Returns a new config object, created by taking a copy of the current config and
        updating it with the settings defined in `custom_settings`. The `custom_settings` can
        be a NomadSettings or a dictionary (in the latter case it must not contain any new keys
        (keys not defined in this NomadSettings). If it does, an exception will be raised.
        '''

        rv = self.copy(deep=True)

        if custom_settings:
            if isinstance(custom_settings, BaseModel):
                for field_name in custom_settings.__fields__.keys():
                    try:
                        setattr(rv, field_name, getattr(custom_settings, field_name))
                    except Exception:
                        raise AssertionError(f'Invalid setting: {field_name}')
            elif isinstance(custom_settings, dict):
                for key, value in custom_settings.items():
                    if value is None:
                        continue
                    try:
                        setattr(rv, key, value)
                    except Exception:
                        raise AssertionError(f'Invalid setting: ({key}: {value})')

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


class Meta(NomadSettings):
    '''
    Metadata about the deployment and how it is presented to clients.
    '''
    version = Field(__version__, description='The NOMAD version string.')
    commit = Field('', description='The source-code commit that this installation\'s NOMAD version is build from.')
    deployment = Field(
        'devel', description='Human-friendly name of this nomad deployment.')
    deployment_url: str = Field(description='The NOMAD deployment\'s url (api url).')
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
    beta: dict = Field({}, description='''
        Additional data that describes how the deployment is labeled as a beta-version in the UI.
    ''')


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


_jupyterhub_config_description = '''
    This setting is forwarded to jupyterhub; refer to the jupyterhub
    [documentation](https://jupyterhub.readthedocs.io/en/stable/api/app.html#).
'''


class NORTHToolMaintainer(BaseModel):
    name: str
    email: str


class NORTHTool(BaseModel):
    image: str
    description: str = None
    short_description: str = None
    cmd: str = None
    privileged: bool = None
    path_prefix: str = None
    mount_path: str = None
    icon: str = None
    file_extensions: List[str] = []
    maintainer: List[NORTHToolMaintainer] = []


class NORTH(NomadSettings):
    '''
    Settings related to the operation of the NOMAD remote tools hub service *north*.
    '''
    enabled: Optional[str] = Field(description='''
        Enables or disables the NORTH API and UI views. This is independent of
        whether you run a jupyter hub or not.
    ''')
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

    nomad_host: str = Field(
        None, description='The NOMAD app host name that spawned containers use.')
    windows = Field(
        True, description='Enable windows OS hacks.')

    tools: Union[str, Dict[str, NORTHTool]] = Field(
        'dependencies/nomad-remote-tools-hub/tools.json',
        description='The available north tools. Either the tools definitions as dict or a path to a .json file.')

    hub_service_api_token: str = Field('secret-token', description='''
        A secret token shared between NOMAD and the NORTH jupyterhub.
        This needs to be the token of an admin service.''')

    @validator('tools', pre=True, always=True)
    def load_tools(cls, v):  # pylint: disable=no-self-argument
        if not isinstance(v, str):
            return v

        # interpret v as file path
        path = v
        if not os.path.exists(path):
            # try to interprete path as relative to project root
            root_path = os.path.join(os.path.dirname(__file__), '../..')
            path = os.path.join(root_path, v)

        with open(path, 'rt') as f:
            return json.load(f)


class RabbitMQ(NomadSettings):
    '''
    Configures how NOMAD is connecting to RabbitMQ.
    '''
    host = Field('localhost', description='The name of the host that runs RabbitMQ.')
    user = Field('rabbitmq', description='The RabbitMQ user that is used to connect.')
    password = Field('rabbitmq', description='The password that is used to connect.')


CELERY_WORKER_ROUTING = 'worker'
CELERY_QUEUE_ROUTING = 'queue'


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


class FS(NomadSettings):
    tmp = '.volumes/fs/tmp'
    staging = '.volumes/fs/staging'
    staging_external: str = None
    public = '.volumes/fs/public'
    public_external: str = None
    north_home = '.volumes/fs/north/users'
    north_home_external: str = None
    local_tmp = '/tmp'
    prefix_size = 2
    archive_version_suffix = 'v1'
    working_directory = os.getcwd()
    external_working_directory: str = None


class Elastic(NomadSettings):
    host = 'localhost'
    port = 9200
    timeout = 60
    bulk_timeout = 600
    bulk_size = 1000
    entries_per_material_cap = 1000
    entries_index = 'nomad_entries_v1'
    materials_index = 'nomad_materials_v1'


class Keycloak(NomadSettings):
    server_url = 'https://nomad-lab.eu/fairdi/keycloak/auth/'
    public_server_url: str = None
    realm_name = 'fairdi_nomad_prod'
    username = 'admin'
    password = 'password'
    client_id = 'nomad_public'
    client_secret: str = None


class Mongo(NomadSettings):
    ''' Connection and usage settings for MongoDB.'''
    host: str = Field('localhost', description='The name of the host that runs mongodb.')
    port: int = Field(27017, description='The port to connect with mongodb.')
    db_name: str = Field('nomad_v1', description='The used mongodb database name.')


class Logstash(NomadSettings):
    enabled = False
    host = 'localhost'
    tcp_port = '5000'
    level: int = logging.DEBUG


class Tests(NomadSettings):
    default_timeout = 60


class Mail(NomadSettings):
    enabled = False
    with_login = False
    host = ''
    port = 8995
    user = ''
    password = ''
    from_address = 'support@nomad-lab.eu'
    cc_address: Optional[str]


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


class Client(NomadSettings):
    user: str = None
    password: str = None
    access_token: str = None
    url = 'http://nomad-lab.eu/prod/v1/api'


class DataCite(NomadSettings):
    mds_host = 'https://mds.datacite.org'
    enabled = False
    prefix = '10.17172'
    user = '*'
    password = '*'


class GitLab(NomadSettings):
    private_token = 'not set'


class Process(NomadSettings):
    store_package_definition_in_mongo = Field(
        False, description='Configures whether to store the corresponding package definition in mongodb.')
    add_definition_id_to_reference = Field(False, description='''
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
    redirect_stdouts: bool = Field(False, description='''
        True will redirect lines to stdout (e.g. print output) that occur during
        processing (e.g. created by parsers or normalizers) as log entries.
    ''')


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


class RFC3161Timestamp(NomadSettings):
    server = Field(
        'http://time.certum.pl/', description='The rfc3161ng timestamping host.')
    cert: str = Field(
        None, description='Path to the optional rfc3161ng timestamping server certificate.')
    hash_algorithm = Field(
        'sha256', description='Hash algorithm used by the rfc3161ng timestamping server.')
    username: str = None
    password: str = None


class BundleExportSettings(NomadSettings):
    include_raw_files: bool = Field(
        True, description='If the raw files should be included in the export')
    include_archive_files: bool = Field(
        True, description='If the parsed archive files should be included in the export')
    include_datasets: bool = Field(
        True, description='If the datasets should be included in the export')


class BundleExport(NomadSettings):
    ''' Controls behaviour related to exporting bundles. '''
    default_cli_bundle_export_path: str = Field(
        './bundles', description='Default path used when exporting bundles using the CLI command.')
    default_settings: BundleExportSettings = Field(
        BundleExportSettings(), description='''
            General default settings.
        ''')
    default_settings_cli: BundleExportSettings = Field(
        None, description='''
            Additional default settings, applied when exporting using the CLI. This allows
            to override some of the settings specified in the general default settings above.
        ''')


class BundleImportSettings(NomadSettings):
    include_raw_files: bool = Field(
        True, description='If the raw files should be included in the import')
    include_archive_files: bool = Field(
        True, description='If the parsed archive files should be included in the import')
    include_datasets: bool = Field(
        True, description='If the datasets should be included in the import')

    include_bundle_info: bool = Field(
        True, description='If the bundle_info.json file should be kept (not necessary but may be nice to have.')
    keep_original_timestamps: bool = Field(
        False, description='''
            If all timestamps (create time, publish time etc) should be imported from
            the bundle.
        ''')
    set_from_oasis: bool = Field(
        True, description='If the from_oasis flag and oasis_deployment_url should be set.')

    delete_upload_on_fail: bool = Field(
        False, description='If False, it is just removed from the ES index on failure.')
    delete_bundle_on_fail: bool = Field(
        True, description='Deletes the source bundle if the import fails.')
    delete_bundle_on_success: bool = Field(
        True, description='Deletes the source bundle if the import succeeds.')
    delete_bundle_include_parent_folder: bool = Field(
        True, description='When deleting the bundle, also include parent folder, if empty.')

    trigger_processing: bool = Field(
        False, description='If the upload should be processed when the import is done (not recommended).')
    process_settings: Reprocess = Field(
        Reprocess(
            rematch_published=True,
            reprocess_existing_entries=True,
            use_original_parser=False,
            add_matched_entries_to_published=True,
            delete_unmatched_published_entries=False
        ), description='''
            When trigger_processing is set to True, these settings control the reprocessing
            behaviour (see the config for `reprocess` for more info). NOTE: reprocessing is
            no longer the recommended method to import bundles.
        '''
    )


class BundleImport(NomadSettings):
    ''' Controls behaviour related to importing bundles. '''
    required_nomad_version: str = Field(
        '1.1.2', description='Minimum  NOMAD version of bundles required for import.')

    default_cli_bundle_import_path: str = Field(
        './bundles', description='Default path used when importing bundles using the CLI command.')

    allow_bundles_from_oasis: bool = Field(
        False, description='If oasis admins can "push" bundles to this NOMAD deployment.')
    allow_unpublished_bundles_from_oasis: bool = Field(
        False, description='If oasis admins can "push" bundles of unpublished uploads.')

    default_settings: BundleImportSettings = Field(
        BundleImportSettings(),
        description='''
            General default settings.
        ''')

    default_settings_cli: BundleImportSettings = Field(
        BundleImportSettings(
            delete_bundle_on_fail=False,
            delete_bundle_on_success=False
        ),
        description='''
            Additional default settings, applied when importing using the CLI. This allows
            to override some of the settings specified in the general default settings above.
        ''')


class Archive(NomadSettings):
    block_size = 256 * 1024
    read_buffer_size = Field(
        256 * 1024, description='GPFS needs at least 256K to achieve decent performance.')
    max_process_number = Field(
        20, description='Maximum number of processes can be assigned to process archive query.')
    min_entries_per_process = Field(
        20, description='Minimum number of entries per process.')


class UISetting(NomadSettings, extra=Extra.forbid):
    '''Extra fields are not allowed in the UI models'''


class OptionsBase(UISetting):
    '''The most basic model for defining the availability of different UI
    options.
    '''
    include: Optional[List[str]] = Field(description='''
        List of included options. If not explicitly defined, all of the options will
        be included by default.
    ''')
    exclude: Optional[List[str]] = Field(description='''
        List of excluded options. Has higher precedence than include.
    ''')


class Options(OptionsBase):
    '''Common configuration class used for enabling/disabling certain UI
    elements and defining the configuration of each element.
    '''
    options: Dict[str, Any] = Field(description='Contains the available options.')


class OptionsSingle(Options):
    '''Represents options where one value can be selected.'''
    selected: str = Field(description='Selected option.')


class OptionsMulti(Options):
    '''Represents options where multiple values can be selected.'''
    selected: List[str] = Field(description='Selected options.')


class UnitSystemEnum(str, Enum):
    CUSTOM = 'Custom'
    SI = 'SI'
    AU = 'AU'


class UnitSystems(UISetting):
    '''Controls the used unit system.'''
    selected: UnitSystemEnum = Field(description='Controls the default unit system.')


class Theme(UISetting):
    '''Theme and identity settings.'''
    title: str = Field(description='Site name in the browser tab.')


class NORTHUI(UISetting):
    '''NORTH (NOMAD Remote Tools Hub) UI configuration.'''
    enabled: bool = Field(True, description='''
        Whether the NORTH tools are available in the UI.
        The default value is read from the root-level NORTH configuration.
    ''')


class Card(UISetting):
    '''Definition for a card shown in the entry overview page.'''
    error: str = Field(description='The error message to show if an error is encountered within the card.')


class Cards(Options):
    '''Contains the overview page card definitions and controls their visibility.'''
    options: Dict[str, Card] = Field(description='Contains the available card options.')


class Entry(UISetting):
    '''Controls the entry visualization.'''
    cards: Cards = Field(description='Controls the cards that are displayed on the entry overview page.')


class Help(UISetting):
    '''Help dialog contents.'''
    title: str = Field(description='Title of the help dialog.')
    content: str = Field(description='Text contents of the help dialog. Supports markdown.')


class Pagination(UISetting):
    order_by: str = Field('upload_create_time', description='Field used for sorting.')
    order: str = Field('desc', description='Sorting order.')
    page_size: int = Field(20, description='Number of results on each page.')


class ModeEnum(str, Enum):
    STANDARD = 'standard'
    SCIENTIFIC = 'scientific'
    SEPARATORS = 'separators'


class Format(UISetting):
    '''Value formatting options.'''
    decimals: int = Field(3, description='Number of decimals to show for numbers.')
    mode: ModeEnum = Field('standard', description='Display mode for numbers.')


class AlignEnum(str, Enum):
    LEFT = 'left'
    RIGHT = 'right'
    CENTER = 'center'


class Column(UISetting):
    '''Option for a column show in the search results.'''
    label: Optional[str] = Field(description='Label shown in the header. Defaults to the quantity name.')
    align: AlignEnum = Field(AlignEnum.LEFT, description='Alignment in the table.')
    unit: Optional[str] = Field(description='''
        Unit to convert to when displaying. If not given will be displayed in
        using the default unit in the active unit system.
    ''')
    format: Optional[Format] = Field(description='Controls the formatting of the values.')


class Columns(OptionsMulti):
    '''
    Contains column definitions, controls their availability and specifies the default
    selection.
    '''
    options: Dict[str, Column] = Field(description='''
        All available column options. Note here that the key must correspond to a
        quantity path that exists in the metadata.
    ''')


class RowActions(UISetting):
    '''Controls the visualization of row actions that are shown at the end of each row.'''
    enabled: bool = Field(True, description='Whether to enable row actions. ')


class RowDetails(UISetting):
    '''
    Controls the visualization of row details that are shown upon pressing the row and
    contain basic details about the entry.
    '''
    enabled: bool = Field(True, description='Whether to show row details.')


class RowSelection(UISetting):
    '''
    Controls the selection of rows. If enabled, rows can be selected and additional
    actions performed on them.
    '''
    enabled: bool = Field(True, description='Whether to show the row selection.')


class Rows(UISetting):
    '''Controls the visualization of rows in the search results.'''
    actions: RowActions
    details: RowDetails
    selection: RowSelection


class FilterMenuActionEnum(str, Enum):
    CHECKBOX = 'checkbox'


class FilterMenuAction(UISetting):
    '''Contains definition for an action in the filter menu.'''
    type: FilterMenuActionEnum = Field(description='Action type.')
    label: str = Field(description='Label to show.')


class FilterMenuActionCheckbox(FilterMenuAction):
    '''Contains definition for checkbox action in the filter menu.'''
    quantity: str = Field(description='Targeted quantity')


class FilterMenuActions(Options):
    '''Contains filter menu action definitions and controls their availability.'''
    options: Dict[str, FilterMenuActionCheckbox] = Field(
        description='Contains options for filter menu actions.'
    )


class FilterMenuSizeEnum(str, Enum):
    S = 's'
    M = 'm'
    L = 'l'
    XL = 'xl'


class FilterMenu(UISetting):
    '''Defines the layout and functionality for a filter menu.'''
    label: Optional[str] = Field(description='Menu label to show in the UI.')
    level: Optional[int] = Field(0, description='Indentation level of the menu.')
    size: Optional[FilterMenuSizeEnum] = Field(FilterMenuSizeEnum.S, description='Width of the menu.')
    actions: Optional[FilterMenuActions]


class FilterMenus(Options):
    '''Contains filter menu definitions and controls their availability.'''
    options: Dict[str, FilterMenu] = Field(description='Contains the available filter menu options.')


class Filters(OptionsBase):
    '''Controls the availability of filters.'''


class Layout(UISetting):
    '''Defines widget size and grid positioning for different breakpoints.'''
    minH: int = Field(description='Minimum height in grid units.')
    minW: int = Field(description='Minimum width in grid units.')
    h: int = Field(description='Height in grid units')
    w: int = Field(description='Width in grid units.')
    x: int = Field(description='Horizontal start location in the grid.')
    y: int = Field(description='Vertical start location in the grid.')


class ScaleEnum(str, Enum):
    POW1 = 'linear'
    POW2 = '1/2'
    POW4 = '1/4'
    POW8 = '1/8'


class BreakpointEnum(str, Enum):
    SM = 'sm'
    MD = 'md'
    LG = 'lg'
    XL = 'xl'
    XXL = 'xxl'


class Widget(UISetting):
    '''Common configuration for all widgets.'''
    type: str = Field(description='Used to identify the widget type.')
    layout: Dict[BreakpointEnum, Layout] = Field(description='''
        Defines widget size and grid positioning for different breakpoints. The
        following breakpoints are supported: `sm`, `md`, `lg`, `xl` and `xxl`.
    ''')


class WidgetTerms(Widget):
    '''Terms widget configuration.'''
    type: Literal['terms'] = Field(description='Set as `terms` to get this widget type.')
    quantity: str = Field(description='Targeted quantity.')
    scale: ScaleEnum = Field(description='Statistics scaling.')
    showinput: bool = Field(True, description='Whether to show text input field.')


class WidgetHistogram(Widget):
    '''Histogram widget configuration.'''
    type: Literal['histogram'] = Field(description='Set as `histogram` to get this widget type.')
    quantity: str = Field(description='Targeted quantity.')
    scale: ScaleEnum = Field(description='Statistics scaling.')
    autorange: bool = Field(
        True,
        description='Whether to automatically set the range according to the data limits.'
    )
    showinput: bool = Field(
        True,
        description='Whether to show input text fields for minimum and maximum value.'
    )
    nbins: int = Field(description='''
        Maximum number of histogram bins. Notice that the actual number of bins
        may be smaller if there are fewer data items available.
    ''')


class WidgetPeriodicTable(Widget):
    '''Periodic table widget configuration.'''
    type: Literal['periodictable'] = Field(description='Set as `periodictable` to get this widget type.')
    quantity: str = Field(description='Targeted quantity.')
    scale: ScaleEnum = Field(description='Statistics scaling.')


class WidgetScatterPlot(Widget):
    '''Scatter plot widget configuration.'''
    type: Literal['scatterplot'] = Field(description='Set as `scatterplot` to get this widget type.')
    x: str = Field(description='X-axis quantity.')
    y: str = Field(description='Y-axis quantity.')
    color: Optional[str] = Field(description='Quantity used for coloring points.')
    size: int = Field(
        1000,
        description='''
        Maximum number of data points to fetch. Notice that the actual number may be less.
        '''
    )
    autorange: bool = Field(
        True,
        description='Whether to automatically set the range according to the data limits.'
    )


# The 'discriminated union' feature of Pydantic is used here:
# https://docs.pydantic.dev/usage/types/#discriminated-unions-aka-tagged-unions
WidgetAnnotated = Annotated[
    Union[WidgetTerms, WidgetHistogram, WidgetScatterPlot, WidgetPeriodicTable],
    Field(discriminator="type")]


class Dashboard(UISetting):
    '''Dashboard configuration.'''
    widgets: List[WidgetAnnotated] = Field(description='List of widgets contained in the dashboard.')  # type: ignore


class ResourceEnum(str, Enum):
    ENTRIES = 'entries'
    MATERIALS = 'materials'


class App(UISetting):
    '''Defines the layout and functionality for an App.'''
    label: str = Field(description='Name of the App.')
    path: str = Field(description='Path used in the browser address bar.')
    resource: ResourceEnum = Field(description='Targeted resource.')
    breadcrumb: str = Field(description='Path displayed in the breadcrumb.')
    category: str = Field(description='Category used to organize Apps in the explore menu.')
    description: str = Field(description='Short description of the App.')
    help: Help = Field(description='Help dialog contents.')
    pagination: Pagination = Field(Pagination(), description='Default result pagination.')
    columns: Columns = Field(description='Controls the columns shown in the results table.')
    rows: Rows = Field(description='Controls the display of entry rows in the results table.')
    filter_menus: FilterMenus = Field(description='Filter menus displayed on the left side of the screen.')
    filters: Optional[Filters] = Field(description='Controls the filters that are available in this app.')
    dashboard: Optional[Dashboard] = Field(description='Default dashboard layout.')
    filters_locked: Optional[dict] = Field(
        description='''
        Fixed query object that is applied for this search context. This filter
        will always be active for this context and will not be displayed to the
        user by default.
        '''
    )


class Apps(Options):
    '''Contains App definitions and controls their availability.'''
    options: Dict[str, App] = Field(description='Contains the available app options.')


class ExampleUploads(OptionsBase):
    '''Controls the availability of example uploads.'''


class UI(UISetting):
    '''Used to customize the user interface.'''
    theme: Theme = Field(
        Theme(**{
            'title': 'NOMAD'
        }),
        description='Controls the site theme and identity.'
    )
    unit_systems: UnitSystems = Field(
        UnitSystems(**{'selected': 'Custom'}),
        description='Controls the available unit systems.'
    )
    entry: Entry = Field(
        Entry(**{
            'cards': {
                'exclude': ['relatedResources'],
                'options': {
                    'sections': {'error': 'Could not render section card.'},
                    'definitions': {'error': 'Could not render definitions card.'},
                    'nexus': {'error': 'Could not render NeXus card.'},
                    'material': {'error': 'Could not render material card.'},
                    'solarcell': {'error': 'Could not render solar cell properties.'},
                    'electronic': {'error': 'Could not render electronic properties.'},
                    'vibrational': {'error': 'Could not render vibrational properties.'},
                    'mechanical': {'error': 'Could not render mechanical properties.'},
                    'thermodynamic': {'error': 'Could not render thermodynamic properties.'},
                    'structural': {'error': 'Could not render structural properties.'},
                    'dynamical': {'error': 'Could not render dynamical properties.'},
                    'geometry_optimization': {'error': 'Could not render geometry optimization.'},
                    'eels': {'error': 'Could not render EELS properties.'},
                    'workflow': {'error': 'Could not render workflow card.'},
                    'references': {'error': 'Could not render references card.'},
                    'relatedResources': {'error': 'Could not render related resources card.'},
                }
            }
        }),
        description='Controls the entry visualization.'
    )
    apps: Apps = Field(
        Apps(**{
            'options': {
                'entries': {
                    'label': 'Entries',
                    'path': 'entries',
                    'resource': 'entries',
                    'breadcrumb': 'Entries',
                    'category': 'All',
                    'description': 'Search entries across all domains',
                    'help': {
                        'title': 'Entries search',
                        'content': inspect.cleandoc(r'''
                            This page allows you to search **entries** within NOMAD.
                            Entries represent any individual data items that have
                            been uploaded to NOMAD, no matter whether they come from
                            theoretical calculations, experiments, lab notebooks or
                            any other source of data. This allows you to perform
                            cross-domain queries, but if you are interested in a
                            specific subfield, you should see if a specific
                            application exists for it in the explore menu to get
                            more details.
                        '''),
                    },
                    'columns': {
                        'selected': [
                            'entry_name',
                            'results.material.chemical_formula_hill',
                            'entry_type',
                            'upload_create_time',
                            'authors'
                        ],
                        'options': {
                            'entry_name': {'label': 'Name', 'align': 'left'},
                            'results.material.chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                            'entry_type': {'label': 'Entry type', 'align': 'left'},
                            'upload_name': {'label': 'Upload name', 'align': 'left'},
                            'upload_id': {'label': 'Upload id', 'align': 'left'},
                            'upload_create_time': {'label': 'Upload time', 'align': 'left'},
                            'authors': {'label': 'Authors', 'align': 'left'},
                            'results.method.method_name': {'label': 'Method name'},
                            'results.method.simulation.program_name': {'label': 'Program name'},
                            'results.method.simulation.dft.basis_set_name': {'label': 'Basis set name'},
                            'results.method.simulation.dft.xc_functional_type': {'label': 'XC Functional Type'},
                            'results.material.structural_type': {'label': 'Dimensionality'},
                            'results.material.symmetry.crystal_system': {'label': 'Crystal system'},
                            'results.material.symmetry.space_group_symbol': {'label': 'Space group symbol'},
                            'results.material.symmetry.space_group_number': {'label': 'Space group number'},
                            'results.eln.lab_ids': {'label': 'Lab IDs'},
                            'results.eln.sections': {'label': 'Sections'},
                            'results.eln.methods': {'label': 'Methods'},
                            'results.eln.tags': {'label': 'Tags'},
                            'results.eln.instruments': {'label': 'Instruments'},
                            'mainfile': {'label': 'Mainfile', 'align': 'left'},
                            'comment': {'label': 'Comment', 'align': 'left'},
                            'references': {'label': 'References', 'align': 'left'},
                            'datasets': {'label': 'Datasets', 'align': 'left'},
                            'published': {'label': 'Access'}
                        }
                    },
                    'rows': {
                        'actions': {'enabled': True},
                        'details': {'enabled': True},
                        'selection': {'enabled': True}
                    },
                    'filter_menus': {
                        'options': {
                            'material': {'label': 'Material', 'level': 0},
                            'elements': {'label': 'Elements / Formula', 'level': 1, 'size': 'xl'},
                            'structure': {'label': 'Structure', 'level': 1},
                            'method': {'label': 'Method', 'level': 0},
                            'dft': {'label': 'DFT', 'level': 1},
                            'gw': {'label': 'GW', 'level': 1},
                            'projection': {'label': 'Projection', 'level': 1},
                            'dmft': {'label': 'DMFT', 'level': 1},
                            'eels': {'label': 'EELS', 'level': 1},
                            'workflow': {'label': 'Workflow', 'level': 0},
                            'molecular_dynamics': {'label': 'Molecular dynamics', 'level': 1},
                            'geometry_optimization': {'label': 'Geometry Optimization', 'level': 1},
                            'properties': {'label': 'Properties', 'level': 0},
                            'electronic': {'label': 'Electronic', 'level': 1},
                            'vibrational': {'label': 'Vibrational', 'level': 1},
                            'mechanical': {'label': 'Mechanical', 'level': 1},
                            'usecases': {'label': 'Use Cases', 'level': 0},
                            'solarcell': {'label': 'Solar Cells', 'level': 1},
                            'author': {'label': 'Author / Origin / Dataset', 'level': 0, 'size': 'm'},
                            'metadata': {'label': 'Visibility / IDs / Schema', 'level': 0},
                            'optimade': {'label': 'Optimade', 'level': 0, 'size': 'm'},
                        }
                    },
                    'filters': {
                        'exclude': ['mainfile', 'entry_name', 'combine']
                    },
                },
                'calculations': {
                    'label': 'Calculations',
                    'path': 'calculations',
                    'resource': 'entries',
                    'breadcrumb': 'Calculations',
                    'category': 'Theory',
                    'description': 'Search calculations',
                    'help': {
                        'title': 'Calculations',
                        'content': inspect.cleandoc(r'''
                            This page allows you to search **calculations** within
                            NOMAD.  Calculations typically come from a specific
                            simulation software that uses an approximate model to
                            investigate and report different physical properties.
                        '''),
                    },
                    'columns': {
                        'selected': [
                            'results.material.chemical_formula_hill',
                            'results.method.simulation.program_name',
                            'results.method.method_name',
                            'upload_create_time',
                            'authors'
                        ],
                        'options': {
                            'results.material.chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                            'results.method.simulation.program_name': {'label': 'Program name'},
                            'results.method.method_name': {'label': 'Method name'},
                            'upload_create_time': {'label': 'Upload time', 'align': 'left'},
                            'authors': {'label': 'Authors', 'align': 'left'},
                            'results.method.simulation.dft.basis_set_name': {'label': 'Basis set name'},
                            'results.method.simulation.dft.xc_functional_type': {'label': 'XC Functional Type'},
                            'results.material.structural_type': {'label': 'Dimensionality'},
                            'results.material.symmetry.crystal_system': {'label': 'Crystal system'},
                            'results.material.symmetry.space_group_symbol': {'label': 'Space group symbol'},
                            'results.material.symmetry.space_group_number': {'label': 'Space group number'},
                            'results.eln.lab_ids': {'label': 'Lab IDs'},
                            'results.eln.sections': {'label': 'Sections'},
                            'results.eln.methods': {'label': 'Methods'},
                            'results.eln.tags': {'label': 'Tags'},
                            'results.eln.instruments': {'label': 'Instruments'},
                            'entry_name': {'label': 'Name', 'align': 'left'},
                            'entry_type': {'label': 'Entry type', 'align': 'left'},
                            'mainfile': {'label': 'Mainfile', 'align': 'left'},
                            'comment': {'label': 'Comment', 'align': 'left'},
                            'references': {'label': 'References', 'align': 'left'},
                            'datasets': {'label': 'Datasets', 'align': 'left'},
                            'published': {'label': 'Access'}
                        }
                    },
                    'rows': {
                        'actions': {'enabled': True},
                        'details': {'enabled': True},
                        'selection': {'enabled': True}
                    },
                    'filter_menus': {
                        'options': {
                            'material': {'label': 'Material', 'level': 0},
                            'elements': {'label': 'Elements / Formula', 'level': 1, 'size': 'xl'},
                            'structure': {'label': 'Structure', 'level': 1},
                            'method': {'label': 'Method', 'level': 0},
                            'dft': {'label': 'DFT', 'level': 1},
                            'gw': {'label': 'GW', 'level': 1},
                            'projection': {'label': 'Projection', 'level': 1},
                            'dmft': {'label': 'DMFT', 'level': 1},
                            'workflow': {'label': 'Workflow', 'level': 0},
                            'molecular_dynamics': {'label': 'Molecular dynamics', 'level': 1},
                            'geometry_optimization': {'label': 'Geometry Optimization', 'level': 1},
                            'properties': {'label': 'Properties', 'level': 0},
                            'electronic': {'label': 'Electronic', 'level': 1},
                            'vibrational': {'label': 'Vibrational', 'level': 1},
                            'mechanical': {'label': 'Mechanical', 'level': 1},
                            'author': {'label': 'Author / Origin / Dataset', 'level': 0, 'size': 'm'},
                            'metadata': {'label': 'Visibility / IDs / Schema', 'level': 0},
                            'optimade': {'label': 'Optimade', 'level': 0, 'size': 'm'},
                        }
                    },
                    'filters': {
                        'exclude': ['mainfile', 'entry_name', 'combine']
                    },
                    'filters_locked': {
                        'quantities': 'results.method.simulation.program_name',
                    },
                },
                'materials': {
                    'label': 'Materials',
                    'path': 'materials',
                    'resource': 'materials',
                    'breadcrumb': 'Materials',
                    'category': 'Theory',
                    'description': 'Search materials that are identified from calculations',
                    'help': {
                        'title': 'Materials',
                        'content': inspect.cleandoc(r'''
                            This page allows you to search **materials** within
                            NOMAD. NOMAD can often automatically detect the material
                            from individual calculations that contain the full
                            atomistic structure and can then group the data by using
                            these detected materials. This allows you to search
                            individual materials which have properties that are
                            aggregated from several entries. Following the link for
                            a specific material will take you to the corresponding
                            [NOMAD Encyclopedia](https://nomad-lab.eu/prod/rae/encyclopedia/#/search)
                            page for that material. NOMAD Encyclopedia is a service
                            that is specifically oriented towards materials property
                            exploration.

                            Notice that by default the properties that you search
                            can be combined from several different entries. If
                            instead you wish to search for a material with an
                            individual entry fullfilling your search criteria,
                            uncheck the **combine results from several
                            entries**-checkbox.
                        '''),
                    },
                    'pagination': {
                        'order_by': 'chemical_formula_hill',
                        'order': 'asc'
                    },
                    'columns': {
                        'selected': [
                            'chemical_formula_hill',
                            'structural_type',
                            'symmetry.structure_name',
                            'symmetry.space_group_number',
                            'symmetry.crystal_system',
                        ],
                        'options': {
                            'chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                            'structural_type': {'label': 'Dimensionality'},
                            'symmetry.structure_name': {'label': 'Structure name'},
                            'symmetry.space_group_number': {'label': 'Space group number'},
                            'symmetry.crystal_system': {'label': 'Crystal system'},
                            'symmetry.space_group_symbol': {'label': 'Space group symbol'},
                            'material_id': {'label': 'Material ID'},
                        }
                    },
                    'rows': {
                        'actions': {'enabled': True},
                        'details': {'enabled': False},
                        'selection': {'enabled': False}
                    },
                    'filter_menus': {
                        'options': {
                            'material': {'label': 'Material', 'level': 0},
                            'elements': {'label': 'Elements / Formula', 'level': 1, 'size': 'xl'},
                            'structure': {'label': 'Structure', 'level': 1},
                            'method': {'label': 'Method', 'level': 0},
                            'dft': {'label': 'DFT', 'level': 1},
                            'gw': {'label': 'GW', 'level': 1},
                            'projection': {'label': 'Projection', 'level': 1},
                            'dmft': {'label': 'DMFT', 'level': 1},
                            'workflow': {'label': 'Workflow', 'level': 0},
                            'molecular_dynamics': {'label': 'Molecular dynamics', 'level': 1},
                            'geometry_optimization': {'label': 'Geometry Optimization', 'level': 1},
                            'properties': {'label': 'Properties', 'level': 0},
                            'electronic': {'label': 'Electronic', 'level': 1},
                            'vibrational': {'label': 'Vibrational', 'level': 1},
                            'mechanical': {'label': 'Mechanical', 'level': 1},
                            'author': {'label': 'Author / Origin / Dataset', 'level': 0, 'size': 'm'},
                            'metadata': {'label': 'Visibility / IDs / Schema', 'level': 0},
                            'optimade': {'label': 'Optimade', 'level': 0, 'size': 'm'},
                            'combine': {
                                'actions': {
                                    'options': {
                                        'combine': {
                                            'type': 'checkbox',
                                            'label': 'Combine results from several entries',
                                            'quantity': 'combine'
                                        }
                                    }
                                }
                            }
                        }
                    },
                    'filters': {
                        'exclude': ['mainfile', 'entry_name']
                    }
                },
                'eln': {
                    'label': 'ELN',
                    'path': 'eln',
                    'resource': 'entries',
                    'breadcrumb': 'ELN',
                    'category': 'Experiment',
                    'description': 'Search electronic lab notebooks',
                    'help': {
                        'title': 'ELN search',
                        'content': inspect.cleandoc(r'''
                            This page allows you to specifically seach **electronic
                            lab notebooks (ELNs)** within NOMAD.  It is very similar
                            to the entries search, but with a reduced filter set and
                            specialized arrangement of default columns.
                        '''),
                    },
                    'columns': {
                        'selected': [
                            'entry_name',
                            'entry_type',
                            'upload_create_time',
                            'authors'
                        ],
                        'options': {
                            'entry_name': {'label': 'Name', 'align': 'left'},
                            'entry_type': {'label': 'Entry type', 'align': 'left'},
                            'upload_create_time': {'label': 'Upload time', 'align': 'left'},
                            'authors': {'label': 'Authors', 'align': 'left'},
                            'results.material.chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                            'results.method.method_name': {'label': 'Method name'},
                            'results.eln.lab_ids': {'label': 'Lab IDs'},
                            'results.eln.sections': {'label': 'Sections'},
                            'results.eln.methods': {'label': 'Methods'},
                            'results.eln.tags': {'label': 'Tags'},
                            'results.eln.instruments': {'label': 'Instruments'},
                            'mainfile': {'label': 'Mainfile', 'align': 'left'},
                            'comment': {'label': 'Comment', 'align': 'left'},
                            'references': {'label': 'References', 'align': 'left'},
                            'datasets': {'label': 'Datasets', 'align': 'left'},
                            'published': {'label': 'Access'}
                        }
                    },
                    'rows': {
                        'actions': {'enabled': True},
                        'details': {'enabled': True},
                        'selection': {'enabled': True}
                    },
                    'filter_menus': {
                        'options': {
                            'material': {'label': 'Material', 'level': 0},
                            'elements': {'label': 'Elements / Formula', 'level': 1, 'size': 'xl'},
                            'eln': {'label': 'Electronic Lab Notebook', 'level': 0},
                            'custom_quantities': {'label': 'User Defined Quantities', 'level': 0, 'size': 'l'},
                            'author': {'label': 'Author / Origin / Dataset', 'level': 0, 'size': 'm'},
                            'metadata': {'label': 'Visibility / IDs / Schema', 'level': 0},
                            'optimade': {'label': 'Optimade', 'level': 0, 'size': 'm'},
                        }
                    },
                    'filters': {
                        'exclude': ['mainfile', 'entry_name', 'combine']
                    },
                    'filters_locked': {
                        'quantities': 'data'
                    }
                },
                'eels': {
                    'label': 'EELS',
                    'path': 'eels',
                    'resource': 'entries',
                    'breadcrumb': 'EELS',
                    'category': 'Experiment',
                    'description': 'Search electron energy loss spectroscopy experiments',
                    'help': {
                        'title': 'EELS',
                        'content': inspect.cleandoc(r'''
                            This page allows you to spefically search **Electron
                            Energy Loss Spectroscopy (EELS) experiments** within
                            NOMAD. It is similar to the entries search, but with a
                            reduced filter set and specialized arrangement of
                            default columns.
                        '''),
                    },
                    'columns': {
                        'selected': [
                            'results.material.chemical_formula_hill',
                            'results.properties.spectroscopy.eels.detector_type',
                            'results.properties.spectroscopy.eels.resolution',
                            'upload_create_time',
                            'authors'
                        ],
                        'options': {
                            'results.material.chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                            'results.properties.spectroscopy.eels.detector_type': {'label': 'Detector type'},
                            'results.properties.spectroscopy.eels.resolution': {'label': 'Resolution'},
                            'upload_create_time': {'label': 'Upload time', 'align': 'left'},
                            'authors': {'label': 'Authors', 'align': 'left'},
                            'results.properties.spectroscopy.eels.min_energy': {},
                            'results.properties.spectroscopy.eels.max_energy': {},
                            'entry_name': {'label': 'Name', 'align': 'left'},
                            'entry_type': {'label': 'Entry type', 'align': 'left'},
                            'mainfile': {'label': 'Mainfile', 'align': 'left'},
                            'comment': {'label': 'Comment', 'align': 'left'},
                            'references': {'label': 'References', 'align': 'left'},
                            'datasets': {'label': 'Datasets', 'align': 'left'},
                            'published': {'label': 'Access'}
                        }
                    },
                    'rows': {
                        'actions': {'enabled': True},
                        'details': {'enabled': True},
                        'selection': {'enabled': True}
                    },
                    'filter_menus': {
                        'options': {
                            'material': {'label': 'Material', 'level': 0},
                            'elements': {'label': 'Elements / Formula', 'level': 1, 'size': 'xl'},
                            'method': {'label': 'Method', 'level': 0},
                            'eels': {'label': 'EELS', 'level': 1},
                            'author': {'label': 'Author / Origin / Dataset', 'level': 0, 'size': 'm'},
                            'metadata': {'label': 'Visibility / IDs / Schema', 'level': 0},
                            'optimade': {'label': 'Optimade', 'level': 0, 'size': 'm'},
                        }
                    },
                    'filters': {
                        'exclude': ['mainfile', 'entry_name', 'combine']
                    },
                    'filters_locked': {
                        'results.method.method_name': 'EELS'
                    }
                },
                'solarcells': {
                    'label': 'Solar Cells',
                    'path': 'solarcells',
                    'resource': 'entries',
                    'breadcrumb': 'Solar Cells',
                    'category': 'Use Cases',
                    'description': 'Search solar cells',
                    'help': {
                        'title': 'Solar cells',
                        'content': inspect.cleandoc(r'''
                            This page allows you to search **solar cell data**
                            within NOMAD. The filter menu on the left and the shown
                            default columns are specifically designed for solar cell
                            exploration. The dashboard directly shows useful
                            interactive statistics about the data.
                        '''),
                    },
                    'pagination': {
                        'order_by': 'results.properties.optoelectronic.solar_cell.efficiency',
                    },
                    'dashboard': {
                        'widgets': [
                            {
                                'type': 'periodictable',
                                'scale': 'linear',
                                'quantity': 'results.material.elements',
                                'layout': {
                                    'xxl': {'minH': 8, 'minW': 12, 'h': 8, 'w': 13, 'y': 0, 'x': 0},
                                    'xl': {'minH': 8, 'minW': 12, 'h': 8, 'w': 12, 'y': 0, 'x': 0},
                                    'lg': {'minH': 8, 'minW': 12, 'h': 8, 'w': 12, 'y': 0, 'x': 0},
                                    'md': {'minH': 8, 'minW': 12, 'h': 8, 'w': 12, 'y': 0, 'x': 0},
                                    'sm': {'minH': 8, 'minW': 12, 'h': 8, 'w': 12, 'y': 16, 'x': 0}
                                },
                            },
                            {
                                'type': 'scatterplot',
                                'autorange': True,
                                'size': 1000,
                                'color': 'results.properties.optoelectronic.solar_cell.short_circuit_current_density',
                                'y': 'results.properties.optoelectronic.solar_cell.efficiency',
                                'x': 'results.properties.optoelectronic.solar_cell.open_circuit_voltage',
                                'layout': {
                                    'xxl': {'minH': 3, 'minW': 3, 'h': 8, 'w': 12, 'y': 0, 'x': 24},
                                    'xl': {'minH': 3, 'minW': 3, 'h': 8, 'w': 9, 'y': 0, 'x': 12},
                                    'lg': {'minH': 3, 'minW': 3, 'h': 6, 'w': 12, 'y': 8, 'x': 0},
                                    'md': {'minH': 3, 'minW': 3, 'h': 6, 'w': 9, 'y': 8, 'x': 0},
                                    'sm': {'minH': 3, 'minW': 3, 'h': 5, 'w': 6, 'y': 0, 'x': 0}
                                },
                            },
                            {
                                'type': 'scatterplot',
                                'autorange': True,
                                'size': 1000,
                                'color': 'results.properties.optoelectronic.solar_cell.device_architecture',
                                'y': 'results.properties.optoelectronic.solar_cell.efficiency',
                                'x': 'results.properties.optoelectronic.solar_cell.open_circuit_voltage',
                                'layout': {
                                    'xxl': {'minH': 3, 'minW': 3, 'h': 8, 'w': 11, 'y': 0, 'x': 13},
                                    'xl': {'minH': 3, 'minW': 3, 'h': 8, 'w': 9, 'y': 0, 'x': 21},
                                    'lg': {'minH': 3, 'minW': 3, 'h': 6, 'w': 12, 'y': 14, 'x': 0},
                                    'md': {'minH': 3, 'minW': 3, 'h': 6, 'w': 9, 'y': 8, 'x': 9},
                                    'sm': {'minH': 3, 'minW': 3, 'h': 5, 'w': 6, 'y': 0, 'x': 6}
                                },
                            },
                            {
                                'type': 'terms',
                                'showinput': True,
                                'scale': 'linear',
                                'quantity': 'results.properties.optoelectronic.solar_cell.device_stack',
                                'layout': {
                                    'xxl': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 8, 'x': 14},
                                    'xl': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 8, 'x': 14},
                                    'lg': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 0, 'x': 12},
                                    'md': {'minH': 3, 'minW': 3, 'h': 4, 'w': 6, 'y': 4, 'x': 12},
                                    'sm': {'minH': 3, 'minW': 3, 'h': 6, 'w': 4, 'y': 10, 'x': 0}
                                },
                            },
                            {
                                'type': 'histogram',
                                'autorange': True,
                                'nbins': 30,
                                'scale': '1/4',
                                'quantity': 'results.properties.optoelectronic.solar_cell.illumination_intensity',
                                'layout': {
                                    'xxl': {'minH': 3, 'minW': 3, 'h': 3, 'w': 8, 'y': 8, 'x': 0},
                                    'xl': {'minH': 3, 'minW': 3, 'h': 3, 'w': 8, 'y': 11, 'x': 0},
                                    'lg': {'minH': 3, 'minW': 3, 'h': 4, 'w': 12, 'y': 12, 'x': 12},
                                    'md': {'minH': 3, 'minW': 3, 'h': 3, 'w': 8, 'y': 17, 'x': 10},
                                    'sm': {'minH': 3, 'minW': 3, 'h': 3, 'w': 8, 'y': 13, 'x': 4}
                                },
                            },
                            {
                                'type': 'terms',
                                'showinput': True,
                                'scale': 'linear',
                                'quantity': 'results.properties.optoelectronic.solar_cell.absorber_fabrication',
                                'layout': {
                                    'xxl': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 8, 'x': 8},
                                    'xl': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 8, 'x': 8},
                                    'lg': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 0, 'x': 18},
                                    'md': {'minH': 3, 'minW': 3, 'h': 4, 'w': 6, 'y': 0, 'x': 12},
                                    'sm': {'minH': 3, 'minW': 3, 'h': 5, 'w': 4, 'y': 5, 'x': 0}
                                },
                            },
                            {
                                'type': 'histogram',
                                'showinput': False,
                                'autorange': False,
                                'nbins': 30,
                                'scale': '1/4',
                                'quantity': 'results.properties.electronic.band_structure_electronic.band_gap.value',
                                'layout': {
                                    'xxl': {'minH': 3, 'minW': 8, 'h': 3, 'w': 8, 'y': 11, 'x': 0},
                                    'xl': {'minH': 3, 'minW': 8, 'h': 3, 'w': 8, 'y': 8, 'x': 0},
                                    'lg': {'minH': 3, 'minW': 8, 'h': 4, 'w': 12, 'y': 16, 'x': 12},
                                    'md': {'minH': 3, 'minW': 8, 'h': 3, 'w': 8, 'y': 14, 'x': 10},
                                    'sm': {'minH': 3, 'minW': 8, 'h': 3, 'w': 8, 'y': 10, 'x': 4}
                                },
                            },
                            {
                                'type': 'terms',
                                'showinput': True,
                                'scale': 'linear',
                                'quantity': 'results.properties.optoelectronic.solar_cell.electron_transport_layer',
                                'layout': {
                                    'xxl': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 8, 'x': 20},
                                    'xl': {'minH': 3, 'minW': 3, 'h': 6, 'w': 5, 'y': 8, 'x': 25},
                                    'lg': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 6, 'x': 18},
                                    'md': {'minH': 3, 'minW': 3, 'h': 6, 'w': 5, 'y': 14, 'x': 0},
                                    'sm': {'minH': 3, 'minW': 3, 'h': 5, 'w': 4, 'y': 5, 'x': 4}
                                },
                            },
                            {
                                'type': 'terms',
                                'showinput': True,
                                'scale': 'linear',
                                'quantity': 'results.properties.optoelectronic.solar_cell.hole_transport_layer',
                                'layout': {
                                    'xxl': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 8, 'x': 26},
                                    'xl': {'minH': 3, 'minW': 3, 'h': 6, 'w': 5, 'y': 8, 'x': 20},
                                    'lg': {'minH': 3, 'minW': 3, 'h': 6, 'w': 6, 'y': 6, 'x': 12},
                                    'md': {'minH': 3, 'minW': 3, 'h': 6, 'w': 5, 'y': 14, 'x': 5},
                                    'sm': {'minH': 3, 'minW': 3, 'h': 5, 'w': 4, 'y': 5, 'x': 8}
                                },
                            }
                        ]
                    },
                    'columns': {
                        'selected': [
                            'results.material.chemical_formula_descriptive',
                            'results.properties.optoelectronic.solar_cell.efficiency',
                            'results.properties.optoelectronic.solar_cell.open_circuit_voltage',
                            'results.properties.optoelectronic.solar_cell.short_circuit_current_density',
                            'results.properties.optoelectronic.solar_cell.fill_factor',
                            'references'
                        ],
                        'options': {
                            'results.material.chemical_formula_descriptive': {'label': 'Descriptive Formula', 'align': 'left'},
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
                            'references': {'label': 'References', 'align': 'left'},
                            'results.material.chemical_formula_hill': {'label': 'Formula', 'align': 'left'},
                            'results.material.structural_type': {'label': 'Dimensionality'},
                            'results.properties.optoelectronic.solar_cell.illumination_intensity': {
                                'label': 'Illum. intensity',
                                'unit': 'W/m**2',
                                'format': {'decimals': 3, 'mode': 'standard'},
                            },
                            'results.eln.lab_ids': {'label': 'Lab IDs'},
                            'results.eln.sections': {'label': 'Sections'},
                            'results.eln.methods': {'label': 'Methods'},
                            'results.eln.tags': {'label': 'Tags'},
                            'results.eln.instruments': {'label': 'Instruments'},
                            'entry_name': {'label': 'Name', 'align': 'left'},
                            'entry_type': {'label': 'Entry type', 'align': 'left'},
                            'mainfile': {'label': 'Mainfile', 'align': 'left'},
                            'upload_create_time': {'label': 'Upload time', 'align': 'left'},
                            'authors': {'label': 'Authors', 'align': 'left'},
                            'comment': {'label': 'Comment', 'align': 'left'},
                            'datasets': {'label': 'Datasets', 'align': 'left'},
                            'published': {'label': 'Access'},
                        }
                    },
                    'rows': {
                        'actions': {'enabled': True},
                        'details': {'enabled': True},
                        'selection': {'enabled': True}
                    },
                    'filter_menus': {
                        'options': {
                            'material': {'label': 'Absorber Material', 'level': 0},
                            'elements': {'label': 'Elements / Formula', 'level': 1, 'size': 'xl'},
                            'structure': {'label': 'Structure', 'level': 1},
                            'electronic': {'label': 'Electronic Properties', 'level': 0},
                            'solarcell': {'label': 'Solar Cell Properties', 'level': 0},
                            'eln': {'label': 'Electronic Lab Notebook', 'level': 0},
                            'custom_quantities': {'label': 'User Defined Quantities', 'level': 0, 'size': 'l'},
                            'author': {'label': 'Author / Origin / Dataset', 'level': 0, 'size': 'm'},
                            'metadata': {'label': 'Visibility / IDs / Schema', 'level': 0},
                            'optimade': {'label': 'Optimade', 'level': 0, 'size': 'm'},
                        }
                    },
                    'filters': {
                        'exclude': ['mainfile', 'entry_name', 'combine']
                    },
                    'filters_locked': {
                        'sections': 'nomad.datamodel.results.SolarCell'
                    }
                },
            }
        }),
        description='Contains the App definitions.'
    )
    north: NORTHUI = Field(
        NORTHUI(),
        description='NORTH (NOMAD Remote Tools Hub) UI configuration.'
    )
    example_uploads: ExampleUploads = Field(
        ExampleUploads(),
        description='Controls the available example uploads.'
    )
