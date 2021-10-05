
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

import click
import sys

from nomad import config
from nomad.cli.cli import cli


@cli.group(help='''The nomad admin commands to do nasty stuff directly on the databases.
                     Remember: With great power comes great responsibility!''')
@click.pass_context
def admin(ctx):
    pass


@admin.command(help='Reset/remove all databases.')
@click.option('--remove', is_flag=True, help='Do not just reset all dbs, but also remove them.')
@click.option('--i-am-really-sure', is_flag=True, help='Must be set for the command to to anything.')
def reset(remove, i_am_really_sure):
    if not i_am_really_sure:
        print('You do not seem to be really sure about what you are doing.')
        sys.exit(1)

    from nomad import infrastructure

    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    infrastructure.reset(remove)


@admin.command(help='Reset all "stuck" in processing uploads and calc in low level mongodb operations.')
@click.option('--zero-complete-time', is_flag=True, help='Sets the complete time to epoch zero.')
def reset_processing(zero_complete_time):
    from datetime import datetime

    from nomad import infrastructure, processing as proc

    infrastructure.setup_mongo()

    def reset_collection(cls):
        in_processing = cls.objects(process_status__in=proc.ProcessStatus.STATUSES_PROCESSING)
        print('%d %s processes need to be reset due to incomplete process' % (in_processing.count(), cls.__name__))
        in_processing.update(
            process_status=proc.ProcessStatus.READY,
            current_process=None,
            current_process_step=None,
            worker_hostname=None,
            celery_task_id=None,
            errors=[], warnings=[],
            complete_time=datetime.fromtimestamp(0) if zero_complete_time else datetime.now())

    reset_collection(proc.Calc)
    reset_collection(proc.Upload)


@admin.command(help='Check and lift embargo of data with expired embargo period.')
@click.option('--dry', is_flag=True, help='Do not lift the embargo, just show what needs to be done.')
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
def lift_embargo(dry, parallel):
    from datetime import datetime
    from dateutil.relativedelta import relativedelta

    from nomad import infrastructure, processing as proc
    from nomad.search import quantity_values

    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    query = dict(with_embargo=True, published=True)

    for upload_id in quantity_values('upload_id', query=query, owner='all'):
        upload = proc.Upload.get(upload_id)
        embargo_length = upload.embargo_length

        if upload.publish_time + relativedelta(months=embargo_length) < datetime.now():
            print('need to lift the embargo of %s (publish_time=%s, embargo=%d)' % (
                upload.upload_id, upload.publish_time, embargo_length))

            if not dry:
                upload.set_upload_metadata({'embargo_length': 0})
    return


@admin.group(help='Generate scripts and commands for nomad operation.')
def ops():
    pass


# @ops.group(help='Tools for managing the DOS similarity data.')
# def similarity():
#     pass


@ops.command(help=('Dump the mongo (calculation metadata) db.'))
@click.option('--restore', is_flag=True, help='Do not dump, but restore.')
def dump(restore: bool):
    from datetime import datetime

    date_str = datetime.utcnow().strftime('%Y_%m_%d')
    print('mongodump --host {} --port {} --db {} -o /backup/fairdi/mongo/{}'.format(
        config.mongo.host, config.mongo.port, config.mongo.db_name, date_str))


@ops.command(help=('Restore the mongo (calculation metadata) db.'))
@click.argument('PATH_TO_DUMP', type=str, nargs=1)
def restore(path_to_dump):
    print('mongorestore --host {} --port {} --db {} {}'.format(
        config.mongo.host, config.mongo.port, config.mongo.db_name, path_to_dump))


@ops.command(help=('Generate an nginx.conf to serve the GUI and proxy pass to API container.'))
@click.option('--prefix', type=str, default=config.services.api_base_path, help='Alter the url path prefix.')
@click.option('--host', type=str, default=config.services.api_host, help='Alter the NOMAD app host.')
@click.option('--port', type=str, default=config.services.api_port, help='Alter the NOMAD port host.')
@click.option('--server/--no-server', default=True, help='Control writing of the outer server {} block. '
              'Useful when conf file is included within another nginx.conf.')
def nginx_conf(prefix, host, port, server):
    prefix = prefix.rstrip('/')
    prefix = '/%s' % prefix.lstrip('/')

    if server:
        print('''server {
    listen        80;
    server_name   www.example.com;
    proxy_set_header Host $host;
        ''')

    print('''
    location / {{
        proxy_pass http://{1}:{2};
    }}

    location ~ {0}\\/?(gui)?$ {{
        rewrite ^ {0}/gui/ permanent;
    }}

    location {0}/gui/ {{
        proxy_intercept_errors on;
        error_page 404 = @redirect_to_index;
        proxy_pass http://{1}:{2};
    }}

    location @redirect_to_index {{
        rewrite ^ {0}/gui/index.html break;
        proxy_pass http://{1}:{2};
    }}

    location ~ \\/gui\\/(service-worker\\.js|meta\\.json)$ {{
        add_header Last-Modified $date_gmt;
        add_header Cache-Control 'no-store, no-cache, must-revalidate, proxy-revalidate, max-age=0';
        if_modified_since off;
        expires off;
        etag off;
        proxy_pass http://{1}:{2};
    }}

    location ~ \\/api\\/uploads\\/?$ {{
        client_max_body_size 35g;
        proxy_request_buffering off;
        proxy_pass http://{1}:{2};
    }}

    location ~ \\/api\\/(raw|archive) {{
        proxy_buffering off;
        proxy_pass http://{1}:{2};
    }}

    location ~ \\/api\\/mirror {{
        proxy_buffering off;
        proxy_read_timeout 600;
        proxy_pass http://{1}:{2};
    }}

    location ~ \\/encyclopedia\\/ {{
        proxy_intercept_errors on;
        error_page 404 = @redirect_to_encyclopedia_index;
        proxy_pass http://{1}:{2};
    }}

    location @redirect_to_encyclopedia_index {{
        rewrite ^ {0}/encyclopedia/index.html break;
        proxy_pass http://{1}:{2};
    }}
'''.format(prefix, host, port))
    if server:
        print('}')


@ops.command(help='Updates the AFLOW prototype information using the latest online version and writes the results to a python module in the given FILEPATH.')
@click.argument('FILEPATH', nargs=1, type=str)
@click.option('--matches-only', is_flag=True, help='Only update the match information that depends on the symmetry analysis settings. Will not perform and online update.')
@click.pass_context
def prototypes_update(ctx, filepath, matches_only):
    from nomad.cli.aflow import update_prototypes
    update_prototypes(ctx, filepath, matches_only)


@ops.command(help='Updates the springer database in nomad.config.normalize.springer_db_path.')
@click.option('--max-n-query', default=10, type=int, help='Number of unsuccessful springer request before returning an error. Default is 10.')
@click.option('--retry-time', default=120, type=int, help='Time in seconds to retry after unsuccessful request. Default is 120.')
def springer_update(max_n_query, retry_time):
    from nomad.cli.admin import springer
    springer.update_springer(max_n_query, retry_time)


# @similarity.command(help='Updates the msgpack file containing the similarity information.')
# @click.option('--dir', "-d", "input_dir", type=str, help='Path of the folder containing the raw similarity information files')
# @click.option('--out', "-o", type=str, help='Path of the output msgpack file.')
# @click.option('--verbose', is_flag=True, help='Enable verbose output.')
# def update(input_dir, out, verbose):
#     from nomad.cli.admin import similarity
#     similarity.update(input_dir, out, verbose)


# @similarity.command(help='Ingests the given similarity information from an msgpack file into MongoDB.')
# @click.option('--in', "-i", "input_path", type=str, help='Path of the ingested msgpack file.')
# @click.option('--batch_size', type=int, default=10000, help='Batch size for MongoDB bulk ingestion.')
# @click.option('--verbose', is_flag=True, help='Enable verbose output.')
# def ingest(input_path, batch_size, verbose):
#     from nomad.cli.admin import similarity
#     similarity.ingest(input_path, batch_size, verbose)


@ops.command(help='Configures the GUIs based on NOMAD config.')
def gui_config():
    import os
    import os.path
    from nomad import config
    import glob
    import shutil

    gui_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../app/flask/static/gui'))
    run_gui_folder = os.path.join(gui_folder, '../.gui_configured')

    # copy
    shutil.rmtree(run_gui_folder, ignore_errors=True)
    shutil.copytree(gui_folder, run_gui_folder)

    # setup the env
    env_js_file = os.path.join(run_gui_folder, 'env.js')
    if not os.path.exists(env_js_file):
        with open(env_js_file, 'wt') as f:
            f.write(('''
window.nomadEnv = {
    'appBase': '%s',
    'keycloakBase': 'https://nomad-lab.eu/fairdi/keycloak/auth/',
    'keycloakRealm': '%s',
    'keycloakClientId': '%s',
    'debug': false,
    'matomoEnabled': false,
    'encyclopediaEnabled': true,
    'oasis': %s
};''' % (
                config.services.api_base_path,
                config.keycloak.realm_name,
                config.keycloak.client_id,
                'true' if config.keycloak.oasis else 'false'
            )))

    # replace base path in all GUI files
    source_file_globs = [
        '**/*.json',
        '**/*.html',
        '**/*.js',
        '**/*.js.map',
        '**/*.css']
    for source_file_glob in source_file_globs:
        source_files = glob.glob(os.path.join(run_gui_folder, source_file_glob), recursive=True)
        for source_file in source_files:
            with open(source_file, 'rt') as f:
                file_data = f.read()
            file_data = file_data.replace('/fairdi/nomad/latest', config.services.api_base_path)
            with open(source_file, 'wt') as f:
                f.write(file_data)

    gui_folder = os.path.abspath(os.path.join(
        os.path.dirname(__file__), '../../app/flask/static/encyclopedia'))

    # setup the env
    conf_js_file = os.path.join(gui_folder, 'conf.js')
    if not os.path.exists(conf_js_file):
        with open(conf_js_file, 'wt') as f:
            f.write(('''
window.nomadEnv = {
    apiRoot: "%s/api/encyclopedia/",
    keycloakBase: "%s",
    keycloakRealm: "%s",
    keycloakClientId: "%s"
};''' % (
                config.services.api_base_path,
                config.keycloak.server_url,
                config.keycloak.realm_name,
                config.keycloak.client_id
            )))
