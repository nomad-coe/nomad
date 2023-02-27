
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

from nomad import config
from nomad.cli.cli import cli
from nomad.cli.dev import get_gui_config


@cli.group(help='''The nomad admin commands to do nasty stuff directly on the databases.
                     Remember: With great power comes great responsibility!''')
@click.pass_context
def admin(ctx):
    pass


@admin.command(help='Reset/remove all databases.')
@click.option('--remove', is_flag=True, help='Do not just reset all dbs, but also remove them.')
@click.option('--i-am-really-sure', is_flag=True, help='Must be set for the command to to anything.')
def reset(remove, i_am_really_sure):
    import sys
    if not i_am_really_sure:
        print('You do not seem to be really sure about what you are doing.')
        sys.exit(1)

    from nomad import infrastructure

    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    infrastructure.reset(remove)

    from nomad.app.resources.main import remove_mongo

    remove_mongo()


@admin.command(help='Reset all uploads and entries "stuck" in processing using level mongodb operations.')
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
            worker_hostname=None,
            celery_task_id=None,
            errors=[], warnings=[],
            complete_time=datetime.fromtimestamp(0) if zero_complete_time else datetime.now())

    reset_collection(proc.Entry)
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
                upload.edit_upload_metadata(
                    edit_request_json=dict(metadata={'embargo_length': 0}),
                    user_id=config.services.admin_user_id)
    return


@admin.group(help='Generate scripts and commands for nomad operation.')
def ops():
    pass


# @ops.group(help='Tools for managing the DOS similarity data.')
# def similarity():
#     pass


@ops.command(help=('Dump the mongo db.'))
@click.option('--restore', is_flag=True, help='Do not dump, but restore.')
def dump(restore: bool):
    from datetime import datetime

    date_str = datetime.utcnow().strftime('%Y_%m_%d')
    print('mongodump --host {} --port {} --db {} -o /backup/fairdi/mongo/{}'.format(
        config.mongo.host, config.mongo.port, config.mongo.db_name, date_str))


@ops.command(help=('Restore the mongo db.'))
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

    location ~ /api/v1/uploads(/?$|.*/raw|.*/bundle?$) {{
        client_max_body_size 35g;
        proxy_request_buffering off;
        proxy_pass http://{1}:{2};
    }}

    location ~ /api/v1/.*/download {{
        proxy_buffering off;
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


@ops.command(help='Configures the GUI for production based on NOMAD config.')
def gui_config():
    import os.path
    from nomad import config
    import glob
    import shutil

    gui_folder = os.path.abspath(
        os.path.join(os.path.dirname(__file__), '../../app/static/gui'))
    run_gui_folder = os.path.join(gui_folder, '../.gui_configured')

    # copy
    shutil.rmtree(run_gui_folder, ignore_errors=True)
    shutil.copytree(gui_folder, run_gui_folder)

    # setup the env
    env_js_file = os.path.join(run_gui_folder, 'env.js')
    if not os.path.exists(env_js_file):
        with open(env_js_file, 'wt') as f:
            conf = get_gui_config(proxy=True)
            f.write(conf)

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


@admin.group(help='Commands for upgrading to a newer NOMAD version')
def upgrade():
    pass


@upgrade.command(
    help='''Converts (upgrades) records from one mongodb and migrates to another.
            Note, it is strongly recommended to run this command with loglevel verbose, i.e.

                nomad -v upgrade migrate-mongo ...

         ''')
@click.option(
    '--host', type=str, default=config.mongo.host,
    help='The mongodb host. By default same as the configured db.')
@click.option(
    '--port', type=int, default=config.mongo.port,
    help='The mongodb port. By default same as the configured db.')
@click.option(
    '--src-db-name', type=str, required=True,
    help='The name of the source database.')
@click.option(
    '--dst-db-name', type=str, default=config.mongo.db_name,
    help='The name of the destination database. By default same as the configured db.')
@click.option(
    '--upload-query', type=str,
    help='An mongo upload query. All uploads matching the query will be included in the migration.')
@click.option(
    '--entry-query', type=str,
    help='An mongo entry query. All uploads with an entry matching the query will be included in the migration.')
@click.option(
    '--ids-from-file', type=str,
    help='''Reads upload IDs from the specified file. Cannot be used together with the
            --upload-query or --entry-query options.
            This can for example be used to retry just the uploads that has previously failed
            (as these ids can be exported to file using --failed-ids-to-file). You can specify both
            --ids-from-file and --failed-ids-to-file at the same time with the same file name.''')
@click.option(
    '--failed-ids-to-file', type=str,
    help='''Write the IDs of failed and skipped uploads to the specified file.
            This can for example be used to subsequently retry just the uploads that failed
            (as these ids can be loaded from file using --ids-from-file). You can specify both
            --ids-from-file and --failed-ids-to-file at the same time with the same file name.''')
@click.option(
    '--upload-update', type=str,
    help='json with updates to apply to all converted uploads')
@click.option(
    '--entry-update', type=str,
    help='json with updates to apply to all converted entries')
@click.option(
    '--overwrite', type=click.Choice(['always', 'if-newer', 'never'], case_sensitive=False), default='never',
    help='''If an upload already exists in the destination db, this option determines whether
            it and its child records should be overwritten with the data from the source db.
            Possible values are "always", "if-newer", "never". Selecting "always" always overwrites,
            "never" never overwrites, and "if-newer" overwrites if the upload either doesn't exist
            in the destination, or it exists but its complete_time (i.e. last time it was
            processed) is older than in the source db.''')
@click.option(
    '--fix-problems', is_flag=True,
    help='''If a minor, fixable problem is encountered, fixes it automaticall; otherwise fail.''')
@click.option(
    '--dry', is_flag=True,
    help='Dry run (not writing anything to the destination database).')
def migrate_mongo(
        host, port, src_db_name, dst_db_name, upload_query, entry_query,
        ids_from_file, failed_ids_to_file, upload_update, entry_update, overwrite, fix_problems, dry):
    import json
    import sys

    from pymongo.database import Database
    from nomad import infrastructure
    from nomad.cli.admin import migrate

    config.mongo.host = host
    config.mongo.port = port
    config.mongo.db_name = dst_db_name
    infrastructure.setup_mongo()

    db_src: Database = infrastructure.mongo_client.get_database(src_db_name)
    db_dst: Database = infrastructure.mongo_client.get_database(dst_db_name)

    if not dry:
        migrate.create_collections_if_needed(db_dst)

    upload_ids = None
    if upload_query and entry_query:
        print('Cannot specify both upload-query and entry-query')
        return -1
    if ids_from_file:
        if upload_query or entry_query:
            print('Cannot specify a query when using --ids-from-file.')
            return -1
        try:
            with open(ids_from_file, 'r') as f:
                upload_ids = [line.strip() for line in f.readlines() if line.strip()]
        except FileNotFoundError:
            print(f'Could not open file {ids_from_file}', file=sys.stderr)
            return -1
    elif upload_query:
        upload_query = json.loads(upload_query)
    elif entry_query:
        entry_query = json.loads(entry_query)

    if upload_update:
        upload_update = json.loads(upload_update)
    if entry_update:
        entry_update = json.loads(entry_update)

    if entry_query:
        print('Querying entries...')
        src_entry_collection = db_src.calc if 'calc' in db_src.collection_names() else db_src.entry
        upload_ids = list(src_entry_collection.distinct('upload_id', entry_query))
    if upload_ids:
        upload_query = {'_id': {'$in': upload_ids}}
    print('Querying uploads...')
    uploads = db_src.upload.find(upload_query)

    migrate.migrate_mongo_uploads(
        db_src, db_dst, uploads, failed_ids_to_file, upload_update, entry_update, overwrite,
        fix_problems, dry)


@admin.command(
    help='''
        Rewrites the existing dataset URLs in existing DOI records with freshly generated dataset URLs.
        This is useful, if the URL layout has changed.''')
@click.argument('DOIs', nargs=-1)
@click.option('--dry', is_flag=True, help='Just test if DOI exists and print is current URL.')
@click.option('--save-existing-records', help='A filename to store the existing DOI records in.')
def rewrite_doi_urls(dois, dry, save_existing_records):
    import json
    import requests

    from nomad.doi import edit_doi_url, _create_dataset_url

    existing_records = []

    if len(dois) == 0:
        from nomad import infrastructure
        from nomad.datamodel import Dataset
        infrastructure.setup_mongo()

        datasets = Dataset.m_def.a_mongo.objects(doi__exists=True)
        dois = [dataset.doi for dataset in datasets]

    try:
        for doi in dois:
            # TODO remove this
            if doi == '10.17172/NOMAD/2016.10.14-1':
                continue

            # check if doi exits
            response = requests.get(f'https://api.datacite.org/dois/{doi}')
            if response.status_code == 404:
                print(f'Cannot rewrite {doi}. DOI does not exist.')
                continue

            data = response.json()
            existing_records.append(data)

            if data['data']['attributes']['url'] == _create_dataset_url(doi):
                print(f'Already up-to-date {doi}')
                continue

            print(f'Updating {doi} ...')
            if not dry:
                edit_doi_url(doi)
    finally:
        if save_existing_records:
            with open(save_existing_records, 'wt') as f:
                json.dump(existing_records, f, indent=2)
