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

import click
import datetime
import elasticsearch.helpers

from nomad import processing as proc, search, datamodel, infrastructure, utils, config

from nomad.cli.cli import cli


@cli.group(help='''The nomad admin commands to do nasty stuff directly on the databases.
                     Remember: With great power comes great responsibility!''')
@click.pass_context
def admin(ctx):
    pass


@admin.command(help='(Re-)index all calcs.')
@click.option('--threads', type=int, default=1, help='Number of threads to use.')
@click.option('--dry', is_flag=True, help='Do not index, just compute entries.')
def index(threads, dry):
    infrastructure.setup_logging()
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    all_calcs = proc.Calc.objects().count()
    print('indexing %d ...' % all_calcs)

    def elastic_updates():
        with utils.ETA(all_calcs, '   index %10d or %10d calcs, ETA %s') as eta:
            for calc in proc.Calc.objects():
                eta.add()
                entry = None
                entry = search.Entry.from_calc_with_metadata(
                    datamodel.CalcWithMetadata(**calc.metadata))
                entry = entry.to_dict(include_meta=True)
                entry['_op_type'] = 'index'
                yield entry

    if dry:
        for _ in elastic_updates():
            pass
    if threads > 1:
        print('  use %d threads' % threads)
        for _ in elasticsearch.helpers.parallel_bulk(
                infrastructure.elastic_client, elastic_updates(), chunk_size=500,
                thread_count=threads):
            pass
    else:
        elasticsearch.helpers.bulk(
            infrastructure.elastic_client, elastic_updates())
    search.refresh()

    print('')
    print('indexing completed')


@admin.group(help='Generate scripts and commands for nomad operation.')
def ops():
    pass


@ops.command(help=('Dump the mongo (calculation metadata) db.'))
@click.option('--restore', is_flag=True, help='Do not dump, but restore.')
def dump(restore: bool):
    date_str = datetime.datetime.utcnow().strftime('%Y_%m_%d')
    print('mongodump --host {} --port {} --db {} -o /backup/fairdi/mongo/{}'.format(
        config.mongo.host, config.mongo.port, config.mongo.db_name, date_str))


@ops.command(help=('Restore the mongo (calculation metadata) db.'))
@click.argument('PATH_TO_DUMP', type=str, nargs=1)
def restore(path_to_dump):
    print('mongorestore --host {} --port {} --db {} {}'.format(
        config.mongo.host, config.mongo.port, config.mongo.db_name, path_to_dump))


@ops.command(help=('Generate an nginx.conf to serve the GUI and proxy pass to API container.'))
@click.option('--prefix', type=str, default='/example_nomad', help='Url path prefix. Default is /example_nomd, can be empty str.')
def nginx_conf(prefix):
    prefix = prefix.rstrip('/')
    prefix = '/%s' % prefix.lstrip('/')

    print('''\
server {{
    listen        80;
    server_name   www.example.com;

    location /{0} {{
        return 301 /example-nomad/gui;
    }}

    location {1}/gui {{
        root      /app/;
        rewrite ^{1}/gui/(.*)$ /nomad/$1 break;
        try_files $uri {1}/gui/index.html;
    }}

    location {1}/gui/service-worker.js {{
        add_header Last-Modified $date_gmt;
        add_header Cache-Control 'no-store, no-cache, must-revalidate, proxy-revalidate, max-age=0';
        if_modified_since off;
        expires off;
        etag off;
        root      /app/;
        rewrite ^{1}/gui/service-worker.js /nomad/service-worker.js break;
    }}

    location {1}/api {{
        proxy_set_header Host $host;
        proxy_pass_request_headers on;
        proxy_pass http://api:8000;
    }}

    location {1}/api/uploads {{
        client_max_body_size 35g;
        proxy_request_buffering off;
        proxy_set_header Host $host;
        proxy_pass_request_headers on;
        proxy_pass http://api:8000;
    }}

    location {1}/api/raw {{
        proxy_buffering off;
        proxy_set_header Host $host;
        proxy_pass_request_headers on;
        proxy_pass http://api:8000;
    }}
}}
    '''.format(prefix.lstrip('/'), prefix))


@ops.command(help=('Generate a proxy pass config for apache2 reverse proxy servers.'))
@click.option('--prefix', type=str, default='uploads', help='The path prefix under which everything is proxy passed.')
@click.option('--host', type=str, default='130.183.207.104', help='The host to proxy to.')
@click.option('--port', type=str, default='30001', help='The port to proxy to.')
def apache_conf(prefix, host, port):
    print('''\
ProxyPass "/{0}" "http://{1}:{2}/{0}"
ProxyPassReverse "/{0}" "http://{1}:{2}/{0}"
<Proxy http://{1}:{2}/{0}>
    ProxyPreserveHost On
    Order deny,allow
    Allow from all
</Proxy>'''.format(prefix, host, port))
