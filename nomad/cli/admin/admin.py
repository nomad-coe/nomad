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
import sys
import io
import re
import uuid
import json
from io import StringIO

import numpy as np
import requests
from ase import Atoms
import ase.io
from bs4 import BeautifulSoup
from matid import SymmetryAnalyzer

from nomad import processing as proc, search, datamodel, infrastructure, utils, config
from nomad.normalizing.structure import get_normalized_wyckoff
from nomad.cli.cli import cli
from nomad import config


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

    infrastructure.setup_logging()
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    infrastructure.reset(remove)


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

    location {1}/gui/meta.json {{
        add_header Last-Modified $date_gmt;
        add_header Cache-Control 'no-store, no-cache, must-revalidate, proxy-revalidate, max-age=0';
        if_modified_since off;
        expires off;
        etag off;
        root      /app/;
        rewrite ^{1}/gui/meta.json /nomad/meta.json break;
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
@click.option('--prefix', type=str, default='app', help='The path prefix under which everything is proxy passed.')
@click.option('--host', type=str, default='130.183.207.104', help='The host to proxy to.')
@click.option('--port', type=str, default='30001', help='The port to proxy to.')
def apache_conf(prefix, host, port):
    print('''\
ProxyPass "/{0}" "http://{1}:{2}/{0}"
ProxyPassReverse "/{0}" "http://{1}:{2}/{0}"
<Proxy http://{1}:{2}/{0}>
    ProxyPreserveHost On
    <IfModule !mod_access_compat.c>
         Require all granted
     </IfModule>
     <IfModule mod_access_compat.c>
         Order allow,deny
         Allow from all
     </IfModule>
</Proxy>

RequestHeader set "X-Forwarded-Proto" expr=%{{REQUEST_SCHEME}}
RequestHeader set "X-Forwarded-SSL" expr=%{{HTTPS}}

ProxyPass /fairdi/keycloak http://{1}:8002/fairdi/keycloak
ProxyPassReverse /fairdi/keycloak http://{1}:8002/fairdi/keycloak
<Proxy http://{1}:8002/app>
     ProxyPreserveHost On
     <IfModule !mod_access_compat.c>
         Require all granted
     </IfModule>
     <IfModule mod_access_compat.c>
         Order allow,deny
         Allow from all
     </IfModule>
</Proxy>

RewriteEngine on
RewriteCond %{QUERY_STRING} ^pid=([^&]+)$
RewriteRule ^/NomadRepository-1.1/views/calculation.zul$ /{0}/gui/entry/pid/%1? [R=301]

AllowEncodedSlashes On
'''.format(prefix, host, port))


def write_prototype_data_file(aflow_prototypes: dict, filepath) -> None:
    """Writes the prototype data file in a compressed format to a python
    module.

    Args:
        aflow_prototypes
    """
    class NoIndent(object):
        def __init__(self, value):
            self.value = value

    class NoIndentEncoder(json.JSONEncoder):
        """A custom JSON encoder that can pretty-print objects wrapped in the
        NoIndent class.
        """
        def __init__(self, *args, **kwargs):
            super(NoIndentEncoder, self).__init__(*args, **kwargs)
            self.kwargs = dict(kwargs)
            del self.kwargs['indent']
            self._replacement_map = {}

        def default(self, o):  # pylint: disable=E0202
            if isinstance(o, NoIndent):
                key = uuid.uuid4().hex
                self._replacement_map[key] = json.dumps(o.value, **self.kwargs)
                return "@@%s@@" % (key,)
            else:
                return super(NoIndentEncoder, self).default(o)

        def encode(self, o):
            result = super(NoIndentEncoder, self).encode(o)
            for k, v in self._replacement_map.items():
                result = result.replace('"@@%s@@"' % (k,), v)
            return result

    prototype_dict = aflow_prototypes["prototypes_by_spacegroup"]
    for prototypes in prototype_dict.values():
        for prototype in prototypes:
            # Save the information back in a prettified form
            prototype["atom_positions"] = NoIndent(prototype["atom_positions"])
            prototype["atom_labels"] = NoIndent(prototype["atom_labels"])
            prototype["lattice_vectors"] = NoIndent(prototype["lattice_vectors"])
            try:
                prototype["normalized_wyckoff_matid"] = NoIndent(prototype["normalized_wyckoff_matid"])
            except KeyError:
                pass

    # Save the updated data
    with io.open(filepath, "w", encoding="utf8") as f:
        json_dump = json.dumps(aflow_prototypes, ensure_ascii=False, indent=4, sort_keys=True, cls=NoIndentEncoder)
        json_dump = re.sub(r"\"(-?\d+(?:[\.,]\d+)?)\"", r'\1', json_dump)  # Removes quotes around numbers
        f.write("# -*- coding: utf-8 -*-\naflow_prototypes = {}\n".format(json_dump))


@ops.command(help='Updates the AFLOW prototype information using the latest online version and writes the results to a python module in the given FILEPATH.')
@click.argument('FILEPATH', nargs=1, type=str)
@click.option('--matches-only', is_flag=True, help='Only update the match information that depends on the symmetry analysis settings. Will not perform and online update.')
@click.pass_context
def prototypes_update(ctx, filepath, matches_only):

    if matches_only:
        from nomad.normalizing.data.aflow_prototypes import aflow_prototypes
    else:
        # The basic AFLOW prototype data is available in a Javascript file. Here we
        # retrieve it and read only the prototype list from it.
        prototypes_file_url = 'http://aflowlib.org/CrystalDatabase/js/table_sort.js'
        r = requests.get(prototypes_file_url, allow_redirects=True)
        datastring = r.content.decode("utf-8")
        datastring = datastring.split('];')[0]
        datastring = datastring.split('= [')[1]
        data = json.loads('[' + datastring + ']')

        newdictarray = []
        n_prototypes = 0
        n_missing = 0
        for protodict in data:
            n_prototypes += 1
            newdict = {}

            # Make prototype plaintext
            prototype = BeautifulSoup(protodict["Prototype"], "html5lib").getText()

            # Add to new dictionary
            newdict['Notes'] = protodict['Notes']
            newdict['Prototype'] = prototype
            newdict['Space Group Symbol'] = protodict['Space Group Symbol']
            newdict['Space Group Number'] = protodict['Space Group Number']
            newdict['Pearsons Symbol'] = protodict['Pearson Symbol']
            newdict['Strukturbericht Designation'] = protodict['Strukturbericht Designation']
            newdict['aflow_prototype_id'] = protodict['AFLOW Prototype']
            newdict['aflow_prototype_url'] = 'http://www.aflowlib.org/CrystalDatabase/' + protodict['href'][2:]

            # Download cif or poscar if possible make ASE Atoms object if possible
            # to obtain labels, positions, cell
            cifurl = 'http://www.aflowlib.org/CrystalDatabase/CIF/' + protodict['href'][2:-5] + '.cif'
            r = requests.get(cifurl, allow_redirects=True)
            cif_str = r.content.decode("utf-8")
            cif_file = StringIO()
            cif_file.write(cif_str)
            cif_file.seek(0)
            try:
                atoms = ase.io.read(cif_file, format='cif')
            except Exception:
                print("Error in getting prototype structure from CIF: {}", format(cifurl))
                # Then try to get structure from POSCAR
                try:
                    poscarurl = 'http://www.aflowlib.org/CrystalDatabase/POSCAR/' + protodict['href'][2:-5] + '.poscar'
                    r = requests.get(poscarurl, allow_redirects=True)
                    poscar_str = r.content.decode("utf-8")
                    poscar_file = StringIO()
                    poscar_file.write(poscar_str)
                    poscar_file.seek(0)
                    atoms = ase.io.read(poscar_file, format='vasp')
                except Exception:
                    print("Error in getting prototype structure from POSCAR: {}".format(poscarurl))
                    print("Could not read prototype structure from CIF or POSCAR file for prototype: {}, {}, ".format(prototype, newdict['aflow_prototype_url']))
                    n_missing += 1
                    continue

            atom_positions = atoms.get_positions()
            atom_labels = atoms.get_chemical_symbols()
            cell = atoms.get_cell()

            newdict['lattice_vectors'] = cell.tolist()
            newdict['atom_positions'] = atom_positions.tolist()
            newdict['atom_labels'] = atom_labels
            newdictarray.append(newdict)

            print("Processed: {}".format(len(newdictarray)))

        # Sort prototype dictionaries by spacegroup and make dictionary
        structure_types_by_spacegroup = {}
        for i_sg in range(1, 231):
            protos_sg = []
            for newdict in newdictarray:
                if newdict['Space Group Number'] == i_sg:
                    protos_sg.append(newdict)
            structure_types_by_spacegroup[i_sg] = protos_sg

        # Wrap in a dictionary that can hold other data, e.g. the symmemtry tolerance parameter.
        aflow_prototypes = {
            "prototypes_by_spacegroup": structure_types_by_spacegroup
        }
        print(
            "Extracted latest AFLOW prototypes online. Total number of "
            "successfully fetched prototypes: {}, missing: {}"
            .format(n_prototypes, n_missing)
        )

    # Update matches
    n_prototypes = 0
    n_failed = 0
    n_unmatched = 0
    prototype_dict = aflow_prototypes["prototypes_by_spacegroup"]

    for aflow_spg_number, prototypes in prototype_dict.items():
        n_prototypes += len(prototypes)
        for prototype in prototypes:

            # Read prototype structure
            pos = np.array(prototype["atom_positions"])
            labels = prototype["atom_labels"]
            cell = np.array(prototype["lattice_vectors"])
            atoms = Atoms(
                symbols=labels,
                positions=pos,
                cell=cell,
                pbc=True
            )

            # Try to first see if the space group can be matched with the one in AFLOW
            tolerance = config.normalize.symmetry_tolerance
            try:
                symm = SymmetryAnalyzer(atoms, tolerance)
                spg_number = symm.get_space_group_number()
                wyckoff_matid = symm.get_wyckoff_letters_conventional()
                norm_system = symm.get_conventional_system()
            except Exception:
                n_failed += 1
            else:
                # If the space group is matched, add the MatID normalized Wyckoff
                # letters to the data.
                if spg_number == aflow_spg_number:
                    atomic_numbers = norm_system.get_atomic_numbers()
                    normalized_wyckoff_matid = get_normalized_wyckoff(atomic_numbers, wyckoff_matid)
                    prototype["normalized_wyckoff_matid"] = normalized_wyckoff_matid
                else:
                    n_unmatched += 1
    print(
        "Updated matches in AFLOW prototype library. Total number of "
        "prototypes: {}, unmatched: {}, failed: {}"
        .format(n_prototypes, n_unmatched, n_failed)
    )

    aflow_prototypes["matid_symmetry_tolerance"] = tolerance

    # Write data file to the specified path
    write_prototype_data_file(aflow_prototypes, filepath)
