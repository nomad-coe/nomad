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

import warnings
import os
import logging
import yaml
from typing import List, Union, Optional, Dict, Any
from pydantic import Field, root_validator, validator, parse_obj_as
from pkg_resources import get_distribution, DistributionNotFound
import pkgutil

try:
    __version__ = get_distribution('nomad-lab').version
except DistributionNotFound:
    # package is not installed
    pass

from .common import (
    ConfigBaseModel,
    Options,
)
from .plugins import Plugins
from .north import NORTH
from .ui import UI

warnings.filterwarnings('ignore', message='numpy.dtype size changed')
warnings.filterwarnings('ignore', message='numpy.ufunc size changed')
warnings.filterwarnings('ignore', category=DeprecationWarning)


def normalize_loglevel(value):
    """Used to normalize log level with pydantic validators."""
    plain_value = value
    if plain_value is None:
        return logging.INFO
    else:
        try:
            return int(plain_value)
        except ValueError:
            return getattr(logging, plain_value)


class Services(ConfigBaseModel):
    """
    Contains basic configuration of the NOMAD services (app, worker, north).
    """

    api_host = Field(
        'localhost',
        description="""
        The external hostname that clients can use to reach this NOMAD installation.
    """,
    )
    api_port = Field(
        8000,
        description="""
        The port used to expose the NOMAD app and api to clients.
    """,
    )
    api_base_path = Field(
        '/fairdi/nomad/latest',
        description="""
        The base path prefix for the NOMAD app and api.
    """,
    )
    api_secret = Field(
        'defaultApiSecret',
        description="""
        A secret that is used to issue download and other tokens.
    """,
    )
    api_timeout = Field(
        600,
        description="""
        If the NOMAD app is run with gunicorn as process manager, this timeout (in s) is passed
        and worker processes will be restarted, if they do not respond in time.
    """,
    )
    https = Field(
        False,
        description="""
        Set to `True`, if external clients are using *SSL* to connect to this installation.
        Requires to setup a reverse-proxy (e.g. the one used in the docker-compose
        based installation) that handles the *SSL* encryption.
    """,
    )
    https_upload = Field(
        False,
        description="""
        Set to `True`, if upload curl commands should suggest the use of SSL for file
        uploads. This can be configured independently of `https` to suggest large file
        via regular HTTP.
    """,
    )
    admin_user_id = Field(
        '00000000-0000-0000-0000-000000000000',
        description="""
        The admin user `user_id`. All users are treated the same; there are no
        particular authorization information attached to user accounts. However, the
        API will grant the user with the given `user_id` more rights, e.g. using the
        `admin` owner setting in accessing data.
    """,
    )

    encyclopedia_base = Field(
        'https://nomad-lab.eu/prod/rae/encyclopedia/#',
        description="""
            This enables links to the given *encyclopedia* installation in the UI.
        """,
    )
    aitoolkit_enabled = Field(
        False,
        description="""
        If true, the UI will show a menu with links to the AI Toolkit notebooks on
        `nomad-lab.eu`.
    """,
    )
    optimade_enabled = Field(
        True, description="""If true, the app will serve the optimade API."""
    )
    dcat_enabled = Field(
        True, description="""If true the app will serve the DCAT API."""
    )
    h5grove_enabled = Field(
        True, description="""If true the app will serve the h5grove API."""
    )

    console_log_level = Field(
        logging.WARNING,
        description="""
        The log level that controls console logging for all NOMAD services (app, worker, north).
        The level is given in Python `logging` log level numbers.
    """,
    )

    upload_limit = Field(
        10,
        description="""
        The maximum allowed unpublished uploads per user. If a user exceeds this
        amount, the user cannot add more uploads.
    """,
    )
    force_raw_file_decoding = Field(
        False,
        description="""
        By default, text raw-files are interpreted with utf-8 encoding. If this fails,
        the actual encoding is guessed. With this setting, we force to assume iso-8859-1
        encoding, if a file is not decodable with utf-8.
    """,
    )
    max_entry_download = Field(
        50000,
        description="""
        There is an inherent limit in page-based pagination with Elasticsearch. If you
        increased this limit with your Elasticsearch, you can also adopt this setting
        accordingly, changing the maximum amount of entries that can be paginated with
        page-base pagination.

        Page-after-value-based pagination is independent and can be used without limitations.
    """,
    )
    unavailable_value = Field(
        'unavailable',
        description="""
        Value that is used in `results` section Enum fields (e.g. system type, spacegroup, etc.)
        to indicate that the value could not be determined.
    """,
    )
    app_token_max_expires_in = Field(
        30 * 24 * 60 * 60,
        description="""
        Maximum expiration time for an app token in seconds. Requests with a higher value
        will be declined.
    """,
    )
    html_resource_http_max_age = Field(
        60,
        description="""
        Used for the max_age cache-control directive on statically served html, js, css
        resources.
    """,
    )
    image_resource_http_max_age = Field(
        30 * 24 * 60 * 60,
        description="""
        Used for the max_age cache-control directive on statically served image
        resources.
    """,
    )

    # Validators
    _console_log_level = validator('console_log_level', allow_reuse=True)(
        normalize_loglevel
    )

    def api_url(
        self,
        ssl: bool = True,
        api: str = 'api',
        api_host: str = None,
        api_port: int = None,
    ):
        """
        Returns the url of the current running nomad API. This is for server-side use.
        This is not the NOMAD url to use as a client, use `nomad.config.client.url` instead.
        """
        if api_port is None:
            api_port = self.api_port
        if api_host is None:
            api_host = self.api_host
        protocol = 'https' if self.https and ssl else 'http'
        host_and_port = api_host
        if api_port not in [80, 443]:
            host_and_port += ':' + str(api_port)
        base_path = self.api_base_path.strip('/')
        return f'{protocol}://{host_and_port}/{base_path}/{api}'


class Meta(ConfigBaseModel):
    """
    Metadata about the deployment and how it is presented to clients.
    """

    version = Field(__version__, description='The NOMAD version string.')
    commit = Field(
        '',
        description="The source-code commit that this installation's NOMAD version is build from.",
    )
    deployment = Field(
        'devel', description='Human-friendly name of this nomad deployment.'
    )
    deployment_url: str = Field(
        None,
        description="The NOMAD deployment's url. If not explicitly set, will default to the (api url) read from the configuration.",
    )
    label: str = Field(
        None,
        description="""
        An additional log-stash data key-value pair added to all logs. Can be used
        to differentiate deployments when analyzing logs.
    """,
    )
    service = Field(
        'unknown nomad service',
        description="""
        Name for the service that is added to all logs. Depending on how NOMAD is
        installed, services get a name (app, worker, north) automatically.
    """,
    )

    name = Field(
        'NOMAD', description='Web-site title for the NOMAD UI.', deprecated=True
    )
    homepage = Field(
        'https://nomad-lab.eu', description='Provider homepage.', deprecated=True
    )
    source_url = Field(
        'https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR',
        description='URL of the NOMAD source-code repository.',
        deprecated=True,
    )

    maintainer_email = Field(
        'markus.scheidgen@physik.hu-berlin.de',
        description='Email of the NOMAD deployment maintainer.',
    )
    beta: dict = Field(
        {},
        description="""
        Additional data that describes how the deployment is labeled as a beta-version in the UI.
    """,
    )


class Oasis(ConfigBaseModel):
    """
    Settings related to the configuration of a NOMAD Oasis deployment.
    """

    is_oasis = Field(
        False,
        description='Set to `True` to indicate that this deployment is a NOMAD Oasis.',
    )
    allowed_users: str = Field(
        None,
        description="""
        A list of usernames or user account emails. These represent a white-list of
        allowed users. With this, users will need to login right-away and only the
        listed users might use this deployment. All API requests must have authentication
        information as well.""",
    )
    uses_central_user_management = Field(
        False,
        description="""
        Set to True to use the central user-management. Typically the NOMAD backend is
        using the configured `keycloak` to access user data. With this, the backend will
        use the API of the central NOMAD (`central_nomad_deployment_url`) instead.
    """,
    )
    central_nomad_deployment_url = Field(
        'https://nomad-lab.eu/prod/v1/api',
        description="""
        The URL of the API of the NOMAD deployment that is considered the *central* NOMAD.
    """,
    )


class RabbitMQ(ConfigBaseModel):
    """
    Configures how NOMAD is connecting to RabbitMQ.
    """

    host = Field('localhost', description='The name of the host that runs RabbitMQ.')
    user = Field('rabbitmq', description='The RabbitMQ user that is used to connect.')
    password = Field('rabbitmq', description='The password that is used to connect.')


CELERY_WORKER_ROUTING = 'worker'
CELERY_QUEUE_ROUTING = 'queue'


class Celery(ConfigBaseModel):
    max_memory = 64e6  # 64 GB
    timeout = 1800  # 1/2 h
    acks_late = False
    routing = CELERY_QUEUE_ROUTING
    priorities = {
        'Upload.process_upload': 5,
        'Upload.delete_upload': 9,
        'Upload.publish_upload': 10,
    }


class FS(ConfigBaseModel):
    tmp = '.volumes/fs/tmp'
    staging = '.volumes/fs/staging'
    staging_external: str = None
    public = '.volumes/fs/public'
    public_external: str = None
    north_home = '.volumes/fs/north/users'
    north_home_external: str = None
    local_tmp = '/tmp'
    prefix_size = 2
    archive_version_suffix: Union[str, List[str]] = Field(
        'v1.2',
        description="""
        This allows to add an additional segment to the names of archive files and
        thereby allows different NOMAD installations to work with the same storage
        directories and raw files, but with separate archives.

        If this is a list, the first string is used. If the file with the first
        string does not exist on read, the system will look for the file with the
        next string, etc.
    """,
    )
    working_directory = os.getcwd()
    external_working_directory: str = None

    @root_validator()
    def __validate(cls, values):  # pylint: disable=no-self-argument
        def get_external_path(path):
            work_dir = values.get('external_working_directory')
            if work_dir and not os.path.isabs(path):
                return os.path.join(work_dir, path)
            return path

        if values.get('staging_external') is None:
            values['staging_external'] = get_external_path(values.get('staging'))

        if values.get('public_external') is None:
            values['public_external'] = get_external_path(values.get('public'))

        if values.get('north_home_external') is None:
            values['north_home_external'] = get_external_path(values.get('north_home'))

        return values


class Elastic(ConfigBaseModel):
    host = 'localhost'
    port = 9200
    timeout = 60
    bulk_timeout = 600
    bulk_size = 1000
    entries_per_material_cap = 1000
    entries_index = 'nomad_entries_v1'
    materials_index = 'nomad_materials_v1'
    username: Optional[str]
    password: Optional[str]


class Keycloak(ConfigBaseModel):
    server_url = 'https://nomad-lab.eu/fairdi/keycloak/auth/'
    public_server_url: str = None
    realm_name = 'fairdi_nomad_prod'
    username = 'admin'
    password = 'password'
    client_id = 'nomad_public'
    client_secret: str = None

    @root_validator()
    def __validate(cls, values):  # pylint: disable=no-self-argument
        if values.get('public_server_url') is None:
            values['public_server_url'] = values.get('server_url')
        return values


class Mongo(ConfigBaseModel):
    """Connection and usage settings for MongoDB."""

    host: str = Field(
        'localhost', description='The name of the host that runs mongodb.'
    )
    port: int = Field(27017, description='The port to connect with mongodb.')
    db_name: str = Field('nomad_v1', description='The used mongodb database name.')
    username: Optional[str]
    password: Optional[str]


class Logstash(ConfigBaseModel):
    enabled = False
    host = 'localhost'
    tcp_port = '5000'
    level: int = logging.DEBUG

    # Validators
    _level = validator('level', allow_reuse=True)(normalize_loglevel)


class Logtransfer(ConfigBaseModel):
    """Configuration of logtransfer and statistics service.

    Note that other configurations are also used within logtransfer

    * class Logstash (Configs: enabled, host, level, tcp_port) such that logs are send to the logstash proxy
    * class Oasis (Config: central_nomad_api_url) address to which the logs are sent to
    * class FS (Config: tmp) path where collected logfiles are stored until they are transferred
    """

    # for logtransfer, see nomad/logtransfer.py
    enable_logtransfer: bool = Field(
        False,
        description='If enabled this starts process that frequently generates logs with statistics.',
    )
    submit_interval: int = Field(
        60 * 60 * 24,
        description='Time interval in seconds after which logs are transferred.',
    )
    max_bytes: int = Field(
        int(1e7),
        description='The size of the logfile in bytes at which the logs are transferred.',
    )
    backup_count: int = Field(
        10,
        description='Number of logfiles stored before oldest rotated logfile is removed.',
    )
    log_filename: str = Field(
        'collectedlogs.txt',
        description='Filename of logfile (located in ".volumes/tmp/").',
    )
    raise_unexpected_exceptions: bool = Field(
        False,
        description='Whether to keep the server alive if an unexpected exception is raised. Set to True for testing.',
    )
    # for statistics (which are submitted to logstash/logtransfer), see nomad/statistics.py
    enable_statistics: bool = Field(
        True,
        description='If enabled this starts a process that frequently generates logs with statistics.',
    )
    statistics_interval: int = Field(
        60 * 60 * 24,
        description='Time interval in seconds in which statistics are logged.',
    )


class Tests(ConfigBaseModel):
    default_timeout = 60
    assume_auth_for_username: str = Field(
        None,
        description=(
            'Will assume that all API calls with no authentication have authentication for '
            'the user with the given username.'
        ),
    )


class Mail(ConfigBaseModel):
    enabled = False
    with_login = False
    host = ''
    port = 8995
    user = ''
    password = ''
    from_address = 'support@nomad-lab.eu'
    cc_address: Optional[str]


class Normalize(ConfigBaseModel):
    normalizers: Options = Field(
        Options(
            include=[
                'OptimadeNormalizer',
                'ResultsNormalizer',
                'MetainfoNormalizer',
            ],
            options=dict(
                PorosityNormalizer='nomad.normalizing.porosity.PorosityNormalizer',
                OptimadeNormalizer='nomad.normalizing.optimade.OptimadeNormalizer',
                ResultsNormalizer='nomad.normalizing.results.ResultsNormalizer',
                MetainfoNormalizer='nomad.normalizing.metainfo.MetainfoNormalizer',
            ),
        )
    )
    system_classification_with_clusters_threshold = Field(
        64,
        description="""
            The system size limit for running the dimensionality analysis. For very
            large systems the dimensionality analysis will get too expensive.
        """,
    )
    clustering_size_limit = Field(
        600,
        description="""
            The system size limit for running the system clustering. For very
            large systems the clustering will get too expensive.
        """,
    )
    symmetry_tolerance = Field(
        0.1,
        description="""
            Symmetry tolerance controls the precision used by spglib in order to
            find symmetries. The atoms are allowed to move this much from their
            symmetry positions in order for spglib to still detect symmetries.
            The unit is angstroms. The value of 0.1 is used e.g. by Materials
            Project according to
            https://pymatgen.org/pymatgen.symmetry.html#pymatgen.symmetry.analyzer.SpacegroupAnalyzer
        """,
    )
    prototype_symmetry_tolerance = Field(
        0.1,
        description="""
            The symmetry tolerance used in aflow prototype matching. Should only be
            changed before re-running the prototype detection.
        """,
    )
    max_2d_single_cell_size = Field(
        7,
        description="""
            Maximum number of atoms in the single cell of a 2D material for it to be
            considered valid.
        """,
    )
    cluster_threshold = Field(
        2.5,
        description="""
            The distance tolerance between atoms for grouping them into the same
            cluster. Used in detecting system type.
        """,
    )
    angle_rounding = Field(
        float(10.0),
        description="""
            Defines the "bin size" for rounding cell angles for the material hash in degree.
        """,
    )
    flat_dim_threshold = Field(
        0.1,
        description="""
            The threshold for a system to be considered "flat". Used e.g. when
            determining if a 2D structure is purely 2-dimensional to allow extra rigid
            transformations that are improper in 3D but proper in 2D.
        """,
    )
    k_space_precision = Field(
        150e6,
        description="""
            The threshold for point equality in k-space. Unit: 1/m.
        """,
    )
    band_structure_energy_tolerance = Field(
        8.01088e-21,
        description="""
            The energy threshold for how much a band can be on top or below the fermi
            level in order to still detect a gap. Unit: Joule.
        """,
    )
    springer_db_path: Optional[str] = Field(
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'normalizing/data/springer.msg'
        )
    )

    @root_validator()
    def __validate(cls, values):  # pylint: disable=no-self-argument
        proto_symmetry_tolerance = values.get('prototype_symmetry_tolerance')
        symmetry_tolerance = values.get('symmetry_tolerance')
        if proto_symmetry_tolerance != symmetry_tolerance:
            raise AssertionError(
                'The AFLOW prototype information is outdated due to changed tolerance '
                'for symmetry detection. Please update the AFLOW prototype information '
                "by running the CLI command 'nomad admin ops prototypes-update --matches-only'"
            )

        springer_db_path = values.get('springer_db_path')
        if springer_db_path and not os.path.exists(springer_db_path):
            values['springer_db_path'] = None

        return values


class Resources(ConfigBaseModel):
    enabled = False
    db_name = 'nomad_v1_resources'
    max_time_in_mongo = Field(
        60 * 60 * 24 * 365.0,
        description="""
            Maxmimum time a resource is stored in mongodb before being updated.
        """,
    )
    download_retries = Field(
        2, description='Number of retries when downloading resources.'
    )
    download_retry_delay = Field(10, description='Delay between retries in seconds')
    max_connections = Field(
        10, description='Maximum simultaneous connections used to download resources.'
    )


class Client(ConfigBaseModel):
    user: str = None
    password: str = None
    access_token: str = None
    url = 'http://nomad-lab.eu/prod/v1/api'


class DataCite(ConfigBaseModel):
    mds_host = 'https://mds.datacite.org'
    enabled = False
    prefix = '10.17172'
    user = '*'
    password = '*'


class GitLab(ConfigBaseModel):
    private_token = 'not set'


class Process(ConfigBaseModel):
    store_package_definition_in_mongo = Field(
        False,
        description='Configures whether to store the corresponding package definition in mongodb.',
    )
    add_definition_id_to_reference = Field(
        False,
        description="""
        Configures whether to attach definition id to `m_def`, note it is different from `m_def_id`.
        The `m_def_id` will be exported with the `with_def_id=True` via `m_to_dict`.
    """,
    )
    write_definition_id_to_archive = Field(
        False, description='Write `m_def_id` to the archive.'
    )
    index_materials = True
    reuse_parser = True
    metadata_file_name = 'nomad'
    metadata_file_extensions = ('json', 'yaml', 'yml')
    auxfile_cutoff = 100
    parser_matching_size = 150 * 80  # 150 lines of 80 ASCII characters per line
    max_upload_size = 32 * (1024**3)
    use_empty_parsers = False
    redirect_stdouts: bool = Field(
        False,
        description="""
        True will redirect lines to stdout (e.g. print output) that occur during
        processing (e.g. created by parsers or normalizers) as log entries.
    """,
    )
    rfc3161_skip_published = False  # skip published entries, regardless of timestamp


class Reprocess(ConfigBaseModel):
    """
    Configures standard behaviour when reprocessing.
    Note, the settings only matter for published uploads and entries. For uploads in
    staging, we always reparse, add newfound entries, and delete unmatched entries.
    """

    rematch_published = True
    reprocess_existing_entries = True
    use_original_parser = False
    add_matched_entries_to_published = True
    delete_unmatched_published_entries = False
    index_individual_entries = False


class RFC3161Timestamp(ConfigBaseModel):
    server = Field(
        'http://zeitstempel.dfn.de', description='The rfc3161ng timestamping host.'
    )
    cert: str = Field(
        None,
        description='Path to the optional rfc3161ng timestamping server certificate.',
    )
    hash_algorithm = Field(
        'sha256',
        description='Hash algorithm used by the rfc3161ng timestamping server.',
    )
    username: str = None
    password: str = None


class BundleExportSettings(ConfigBaseModel):
    include_raw_files: bool = Field(
        True, description='If the raw files should be included in the export'
    )
    include_archive_files: bool = Field(
        True, description='If the parsed archive files should be included in the export'
    )
    include_datasets: bool = Field(
        True, description='If the datasets should be included in the export'
    )


class BundleExport(ConfigBaseModel):
    """Controls behaviour related to exporting bundles."""

    default_cli_bundle_export_path: str = Field(
        './bundles',
        description='Default path used when exporting bundles using the CLI command.',
    )
    default_settings: BundleExportSettings = Field(
        BundleExportSettings(),
        description="""
            General default settings.
        """,
    )
    default_settings_cli: BundleExportSettings = Field(
        None,
        description="""
            Additional default settings, applied when exporting using the CLI. This allows
            to override some of the settings specified in the general default settings above.
        """,
    )


class BundleImportSettings(ConfigBaseModel):
    include_raw_files: bool = Field(
        True, description='If the raw files should be included in the import'
    )
    include_archive_files: bool = Field(
        True, description='If the parsed archive files should be included in the import'
    )
    include_datasets: bool = Field(
        True, description='If the datasets should be included in the import'
    )

    include_bundle_info: bool = Field(
        True,
        description='If the bundle_info.json file should be kept (not necessary but may be nice to have.',
    )
    keep_original_timestamps: bool = Field(
        False,
        description="""
            If all timestamps (create time, publish time etc) should be imported from
            the bundle.
        """,
    )
    set_from_oasis: bool = Field(
        True,
        description='If the from_oasis flag and oasis_deployment_url should be set.',
    )

    delete_upload_on_fail: bool = Field(
        False, description='If False, it is just removed from the ES index on failure.'
    )
    delete_bundle_on_fail: bool = Field(
        True, description='Deletes the source bundle if the import fails.'
    )
    delete_bundle_on_success: bool = Field(
        True, description='Deletes the source bundle if the import succeeds.'
    )
    delete_bundle_include_parent_folder: bool = Field(
        True,
        description='When deleting the bundle, also include parent folder, if empty.',
    )

    trigger_processing: bool = Field(
        False,
        description='If the upload should be processed when the import is done (not recommended).',
    )
    process_settings: Reprocess = Field(
        Reprocess(
            rematch_published=True,
            reprocess_existing_entries=True,
            use_original_parser=False,
            add_matched_entries_to_published=True,
            delete_unmatched_published_entries=False,
        ),
        description="""
            When trigger_processing is set to True, these settings control the reprocessing
            behaviour (see the config for `reprocess` for more info). NOTE: reprocessing is
            no longer the recommended method to import bundles.
        """,
    )


class BundleImport(ConfigBaseModel):
    """Controls behaviour related to importing bundles."""

    required_nomad_version: str = Field(
        '1.1.2', description='Minimum  NOMAD version of bundles required for import.'
    )

    default_cli_bundle_import_path: str = Field(
        './bundles',
        description='Default path used when importing bundles using the CLI command.',
    )

    allow_bundles_from_oasis: bool = Field(
        False,
        description='If oasis admins can "push" bundles to this NOMAD deployment.',
    )
    allow_unpublished_bundles_from_oasis: bool = Field(
        False, description='If oasis admins can "push" bundles of unpublished uploads.'
    )

    default_settings: BundleImportSettings = Field(
        BundleImportSettings(),
        description="""
            General default settings.
        """,
    )

    default_settings_cli: BundleImportSettings = Field(
        BundleImportSettings(
            delete_bundle_on_fail=False, delete_bundle_on_success=False
        ),
        description="""
            Additional default settings, applied when importing using the CLI. This allows
            to override some of the settings specified in the general default settings above.
        """,
    )


class Archive(ConfigBaseModel):
    block_size = Field(
        1 * 2**20,
        description='In case of using blocked TOC, this is the size of each block.',
    )
    read_buffer_size = Field(
        1 * 2**20,
        description='GPFS needs at least 256K to achieve decent performance.',
    )
    copy_chunk_size = Field(
        16 * 2**20,
        description="""
        The chunk size of every read of binary data.
        It is used to copy data from one file to another.
        A small value will result in more syscalls, a large value will result in higher peak memory usage.
        """,
    )
    toc_depth = Field(10, description='Depths of table of contents in the archive.')
    use_new_writer = True  # todo: to be removed
    small_obj_optimization_threshold = Field(
        1 * 2**20,
        description="""
        For any child of lists/dicts whose encoded size is smaller than this value, no TOC will be generated.""",
    )
    fast_loading = Field(
        True,
        description="""
        When enabled, this flag determines whether to read the whole dict/list at once
        when a certain mount of children has been visited.
        This reduces the number of syscalls although data may be repeatedly read.
        Otherwise, always read children one by one. This may slow down the loading as more syscalls are needed.
        """,
    )
    fast_loading_threshold = Field(
        0.6,
        description="""
        If the fraction of children that have been visited is less than this threshold, fast loading will be used.
        """,
    )
    trivial_size = Field(
        20,
        description="""
        To identify numerical lists.
        """,
    )


class Config(ConfigBaseModel):
    """Model for the NOMAD configuration."""

    services: Services = Services()
    meta: Meta = Meta()
    oasis: Oasis = Oasis()
    north: NORTH = NORTH()
    rabbitmq: RabbitMQ = RabbitMQ()
    celery: Celery = Celery()
    fs: FS = FS()
    elastic: Elastic = Elastic()
    keycloak: Keycloak = Keycloak()
    mongo: Mongo = Mongo()
    logstash: Logstash = Logstash()
    logtransfer: Logtransfer = Logtransfer()
    tests: Tests = Tests()
    mail: Mail = Mail()
    normalize: Normalize = Normalize()
    resources: Resources = Resources()
    client: Client = Client()
    datacite: DataCite = DataCite()
    gitlab: GitLab = GitLab()
    process: Process = Process()
    reprocess: Reprocess = Reprocess()
    rfc3161_timestamp: RFC3161Timestamp = RFC3161Timestamp()
    bundle_export: BundleExport = BundleExport()
    bundle_import: BundleImport = BundleImport()
    archive: Archive = Archive()
    ui: UI = UI()
    plugins: Optional[Plugins]

    def api_url(
        self,
        ssl: bool = True,
        api: str = 'api',
        api_host: str = None,
        api_port: int = None,
    ):
        """
        Returns the url of the current running nomad API. This is for server-side use.
        This is not the NOMAD url to use as a client, use `nomad.config.client.url` instead.
        """
        return self.services.api_url(ssl, api, api_host, api_port)

    def gui_url(self, page: str = None):
        base = self.api_url(True)[:-3]
        if base.endswith('/'):
            base = base[:-1]

        if page is not None:
            return '%s/gui/%s' % (base, page)

        return '%s/gui' % base

    def rabbitmq_url(self):
        return 'pyamqp://%s:%s@%s//' % (
            self.rabbitmq.user,
            self.rabbitmq.password,
            self.rabbitmq.host,
        )

    def north_url(self, ssl: bool = True):
        return self.api_url(
            ssl=ssl,
            api='north',
            api_host=self.north.hub_host,
            api_port=self.north.hub_port,
        )

    def hub_url(self):
        return f'http://{self.north.hub_host}:{self.north.hub_port}{self.services.api_base_path}/north/hub'

    @root_validator()
    def __validate(cls, values):  # pylint: disable=no-self-argument
        services = values.get('services')
        deployment_url = values.get('meta').deployment_url
        if not deployment_url:
            values.get('meta').deployment_url = services.api_url()
        north = values.get('north')
        ui = values.get('ui')
        if ui:
            if north:
                values['ui'].north.enabled = north.enabled
            if services:
                values[
                    'ui'
                ].app_base = f'{"https" if services.https else "http"}://{services.api_host}:{services.api_port}{values["services"].api_base_path.rstrip("/")}'
            if services and north:
                values[
                    'ui'
                ].north_base = f'{"https" if services.https else "http"}://{north.hub_host}:{north.hub_port}{services.api_base_path.rstrip("/")}/north'

        return values

    def load_plugins(self):
        """Used to lazy-load the plugins. We cannot instantiate the plugins
        during the initialization of the nomad.config package, because it may
        trigger circular dependency errors for plugins using the config
        (pkgutil.get_loader will run code in the package root __init__). Instead
        this function should be called to instantiate the plugins before the
        nomad application is started.

        TODO: Once we migrate to Pydantic v2, we should add the computed_field +
        cached_property decorator to the 'plugins' field instead of using this
        function.
        """
        if self.plugins is None:
            from nomad.config import _plugins
            from nomad.config.models.plugins import Plugin

            def load_plugin(values: Dict[str, Any]):
                """Loads plugin metadata from nomad_plugin.yaml"""
                python_package = values.get('python_package')
                if not python_package:
                    raise ValueError('Python plugins must provide a python_package.')

                package_path = values.get('package_path')
                if package_path is None:
                    try:
                        # We try to deduce the package path from the top-level package
                        package_path_segments = python_package.split('.')
                        root_package = package_path_segments[0]
                        package_dirs = package_path_segments[1:]
                        package_path = os.path.join(
                            os.path.dirname(
                                pkgutil.get_loader(root_package).get_filename()  # type: ignore
                            ),
                            *package_dirs,
                        )
                        if not os.path.isdir(package_path):
                            # We could not find it this way. Let's try to official way
                            package_path = os.path.dirname(
                                pkgutil.get_loader(python_package).get_filename()  # type: ignore
                            )
                    except Exception as e:
                        raise ValueError(
                            f'The python package {python_package} cannot be loaded.', e
                        )
                    values['package_path'] = package_path

                metadata_path = os.path.join(package_path, 'nomad_plugin.yaml')
                if os.path.exists(metadata_path):
                    try:
                        with open(metadata_path, 'r', encoding='UTF-8') as f:
                            metadata = yaml.load(f, Loader=yaml.SafeLoader)
                    except Exception as e:
                        raise ValueError(
                            f'Cannot load plugin metadata file {metadata_path}.', e
                        )

                    for key, value in metadata.items():
                        if key not in values:
                            values[key] = value

                return values

            for key, plugin in _plugins['options'].items():
                _plugins['options'][key] = parse_obj_as(Plugin, load_plugin(plugin))
            self.plugins = Plugins.parse_obj(_plugins)
