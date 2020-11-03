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
import typing
import click
import datetime
import elasticsearch_dsl
import elasticsearch
import sys
import threading

from nomad import processing as proc, search, datamodel, infrastructure, utils, config
from nomad.cli.cli import cli


def __run_parallel(
        uploads, parallel: int, callable, label: str):
    if isinstance(uploads, (tuple, list)):
        uploads_count = len(uploads)

    else:
        uploads_count = uploads.count()
        uploads = list(uploads)  # copy the whole mongo query set to avoid cursor timeouts

    cv = threading.Condition()
    threads: typing.List[threading.Thread] = []

    state = dict(
        completed_count=0,
        skipped_count=0,
        available_threads_count=parallel)

    logger = utils.get_logger(__name__)

    print('%d uploads selected, %s ...' % (uploads_count, label))

    def process_upload(upload: proc.Upload):
        logger.info('%s started' % label, upload_id=upload.upload_id)

        completed = False
        if callable(upload, logger):
            completed = True

        with cv:
            state['completed_count'] += 1 if completed else 0
            state['skipped_count'] += 1 if not completed else 0
            state['available_threads_count'] += 1

            print(
                '   %s %s and skipped %s of %s uploads' %
                (label, state['completed_count'], state['skipped_count'], uploads_count))

            cv.notify()

    for upload in uploads:
        with cv:
            cv.wait_for(lambda: state['available_threads_count'] > 0)
            state['available_threads_count'] -= 1
            thread = threading.Thread(target=lambda: process_upload(upload))
            threads.append(thread)
            thread.start()

    for thread in threads:
        thread.join()


def __run_processing(
        uploads, parallel: int, process, label: str, reprocess_running: bool = False):

    def run_process(upload, logger):
        if upload.process_running and not reprocess_running:
            logger.warn(
                'cannot trigger %s, since the upload is already/still processing' % label,
                current_process=upload.current_process,
                current_task=upload.current_task, upload_id=upload.upload_id)
            return False
        else:
            upload.reset()
            process(upload)
            upload.block_until_complete(interval=.5)

            if upload.tasks_status == proc.FAILURE:
                logger.info('%s with failure' % label, upload_id=upload.upload_id)

            logger.info('%s complete' % label, upload_id=upload.upload_id)
            return True

    __run_parallel(uploads, parallel=parallel, callable=run_process, label=label)


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

    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    infrastructure.reset(remove)


@admin.command(help='Reset all "stuck" in processing uploads and calc in low level mongodb operations.')
@click.option('--zero-complete-time', is_flag=True, help='Sets the complete time to epoch zero.')
def reset_processing(zero_complete_time):
    infrastructure.setup_mongo()

    def reset_collection(cls):
        in_processing = cls.objects(process_status__in=[proc.PROCESS_RUNNING, proc.base.PROCESS_CALLED])
        print('%d %s processes need to be reset due to incomplete process' % (in_processing.count(), cls.__name__))
        in_processing.update(
            process_status=None,
            current_process=None,
            worker_hostname=None,
            celery_task_id=None,
            errors=[], warnings=[],
            complete_time=datetime.datetime.fromtimestamp(0) if zero_complete_time else datetime.datetime.now(),
            current_task=None,
            tasks_status=proc.base.CREATED)

        in_tasks = cls.objects(tasks_status__in=[proc.PENDING, proc.RUNNING])
        print('%d %s processes need to be reset due to incomplete tasks' % (in_tasks.count(), cls.__name__))
        in_tasks.update(
            current_task=None,
            tasks_status=proc.base.CREATED,
            errors=[], warnings=[],
            complete_time=datetime.datetime.fromtimestamp(0) if zero_complete_time else datetime.datetime.now())

    reset_collection(proc.Calc)
    reset_collection(proc.Upload)


@admin.command(help='Check and lift embargo of data with expired embargo period.')
@click.option('--dry', is_flag=True, help='Do not lift the embargo, just show what needs to be done.')
@click.option('--parallel', default=1, type=int, help='Use the given amount of parallel processes. Default is 1.')
def lift_embargo(dry, parallel):
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    request = search.SearchRequest()
    request.q = elasticsearch_dsl.Q('term', with_embargo=True) & elasticsearch_dsl.Q('term', published=True)
    request.quantity('upload_id', 1000)
    result = request.execute()

    uploads_to_repack = []
    for upload_id in result['quantities']['upload_id']['values']:
        upload = proc.Upload.get(upload_id)
        embargo_length = upload.embargo_length
        if embargo_length is None:
            embargo_length = 36
            upload.embargo_length = 36

        if upload.upload_time + datetime.timedelta(days=int(embargo_length * 365 / 12)) < datetime.datetime.now():
            print('need to lift the embargo of %s (upload_time=%s, embargo=%d)' % (
                upload.upload_id, upload.upload_time, embargo_length))

            if not dry:
                proc.Calc._get_collection().update_many(
                    {'upload_id': upload_id},
                    {'$set': {'metadata.with_embargo': False}})
                uploads_to_repack.append(upload)
                upload.save()

                with upload.entries_metadata() as entries:
                    search.index_all(entries)

    if not dry:
        __run_processing(uploads_to_repack, parallel, lambda upload: upload.re_pack(), 're-packing')


@admin.command(help='(Re-)index all calcs.')
@click.option('--threads', type=int, default=1, help='Number of threads to use.')
@click.option('--dry', is_flag=True, help='Do not index, just compute entries.')
def index(threads, dry):
    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    all_calcs = proc.Calc.objects().count()
    print('indexing %d ...' % all_calcs)

    def elastic_updates():
        with utils.ETA(all_calcs, '   index %10d or %10d calcs, ETA %s') as eta:
            for calc in proc.Calc.objects():
                eta.add()
                entry_metadata = datamodel.EntryMetadata.m_from_dict(calc.metadata)
                entry = entry_metadata.a_elastic.create_index_entry().to_dict(include_meta=True)
                entry['_op_type'] = 'index'
                yield entry

    if dry:
        for _ in elastic_updates():
            pass
    else:
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


@admin.command()
@click.option('--threads', type=int, default=1, help='Number of threads to use.')
@click.option('--code', multiple=True, type=str, help='Index only calculcations of given codes.')
@click.option('--dry', is_flag=True, help='Do not index, just compute entries.')
@click.option('--in-place', is_flag=True, default=False, help='Perform indexing in the current elastic search index. Meant only for small reindex operations.')
@click.option('-n', type=int, default=None, help='Number of calculations to process. Leave undefined to process all calculations.')
@click.option('--source',
              type=click.Choice(['mongo', 'es'], case_sensitive=True))
def index_materials(threads, code, dry, in_place, n, source):
    """(Re-)index all materials.

    This command will completely rebuild the materials index. The index is
    built from the material metainfo stored in MongoDB. The materials index can
    be used normally during the reindexing.
    """
    from nomad.datamodel.material import Material, Calculation
    from nomad.datamodel.encyclopedia import EncyclopediaMetadata
    from nomad.search import material_document
    from nomad.datamodel.material import Material, Calculation, Method, Properties, IdealizedStructure, Energies, Workflow, Bulk

    chunk_size = 500
    infrastructure.setup_mongo()
    client = infrastructure.setup_elastic()

    # In order to do the reindexing with zero downtime, two different indices
    # are rotated and an alias is used
    old_index_name = list(client.indices.get(config.elastic.materials_index_name).keys())[0]
    if in_place:
        target_index_name = old_index_name
    else:
        if old_index_name == config.elastic.materials_index_name + "_a":
            target_index_name = config.elastic.materials_index_name + "_b"
        elif old_index_name == config.elastic.materials_index_name + "_b":
            target_index_name = config.elastic.materials_index_name + "_a"
        else:
            raise ValueError(
                "Unrecognized index name accociated with the alias {}"
                .format(config.elastic.materials_index_name)
            )

    if source == "mongo":
        all_calcs = proc.Calc.objects().count()
        print('indexing materials from %d calculations ...' % all_calcs)

        # Bulk update
        def elastic_updates():
            with utils.ETA(all_calcs, '   index %10d of %10d calcs, ETA %s') as eta:
                mongo_db = infrastructure.mongo_client[config.mongo.db_name]
                mongo_collection = mongo_db['archive']
                i_calc = 0
                for mongo_archive in mongo_collection.find():
                    i_calc += 1
                    if n is not None:
                        if i_calc > n:
                            return
                    eta.add()

                    # Do not process entries that do not have the material
                    # information
                    try:
                        status = mongo_archive["section_metadata"]["encyclopedia"]["status"]
                        if status != EncyclopediaMetadata.status.type.success:
                            raise AttributeError
                    except (KeyError, AttributeError, IndexError):
                        continue

                    # Create material information
                    metadata = mongo_archive["section_metadata"]
                    encyclopedia = EncyclopediaMetadata.m_from_dict(metadata["encyclopedia"])
                    dft = metadata["dft"]
                    material: Material = Material()
                    material.material_id = encyclopedia.material.material_id
                    material.material_type = encyclopedia.material.material_type
                    material.material_name = encyclopedia.material.material_name
                    material.material_classification = encyclopedia.material.material_classification
                    material.formula = encyclopedia.material.formula
                    material.formula_reduced = encyclopedia.material.formula_reduced
                    material.species_and_counts = encyclopedia.material.species_and_counts
                    material.species = encyclopedia.material.species
                    enc_bulk = encyclopedia.material.bulk
                    if enc_bulk:
                        bulk = Bulk.m_from_dict(enc_bulk.m_to_dict())
                        material.m_add_sub_section(Material.bulk, bulk)

                    # Create calculation info for this entry
                    calc = Calculation()
                    calc.calc_id = metadata["calc_id"]
                    calc.upload_id = metadata["upload_id"]
                    mongo_calc = proc.Calc.get(calc.calc_id)
                    calc.published = mongo_calc["metadata"]["published"]
                    calc.with_embargo = mongo_calc["metadata"]["with_embargo"]
                    calc.owners = [mongo_calc["metadata"]["uploader"]] + mongo_calc["metadata"]["shared_with"]
                    enc_idealized_structure = encyclopedia.material.idealized_structure
                    idealized_structure = IdealizedStructure()
                    cell_volume = enc_idealized_structure.cell_volume
                    if cell_volume is not None:
                        idealized_structure.cell_volume = cell_volume
                    idealized_structure.lattice_parameters = enc_idealized_structure.lattice_parameters
                    calc.m_add_sub_section(Calculation.idealized_structure, idealized_structure)
                    enc_method = encyclopedia.method
                    method = Method.m_from_dict(enc_method.m_to_dict())
                    method.program_name = dft["code_name"]
                    method.program_version = dft["code_version"]
                    method.basis_set = dft["basis_set"]
                    calc.m_add_sub_section(Calculation.method, method)
                    enc_props = encyclopedia.properties

                    # Properties may not exist at all
                    if enc_props is not None:
                        properties = Properties()

                        # Energies may not be present in all calculations
                        try:
                            energies = Energies.m_from_dict(enc_props.energies.m_to_dict())
                            properties.m_add_sub_section(Properties.energies, energies)
                        except AttributeError:
                            pass

                        properties.has_electronic_dos = enc_props.electronic_dos is not None
                        properties.has_electronic_band_structure = enc_props.electronic_band_structure is not None
                        properties.has_thermodynamical_properties = enc_props.thermodynamical_properties is not None
                        atomic_density = enc_props.atomic_density
                        if atomic_density is not None:
                            properties.atomic_density = atomic_density
                        mass_density = enc_props.mass_density
                        if mass_density is not None:
                            properties.mass_density = mass_density
                        band_gap = enc_props.band_gap
                        if band_gap is not None:
                            properties.band_gap = band_gap
                        band_gap_direct = enc_props.band_gap_direct
                        if band_gap_direct is not None:
                            properties.band_gap_direct = band_gap_direct
                        calc.m_add_sub_section(Calculation.properties, properties)

                    workflow = Workflow()
                    workflow.workflow_type = encyclopedia.calculation.calculation_type
                    calc.m_add_sub_section(Calculation.workflow, workflow)
                    material.m_add_sub_section(Material.calculations, calc)

                    # Update entry that inserts the full material info if entry
                    # does not exists, otherwise only adds the calculation into the
                    # nested subdocument
                    entry = {}
                    entry['_op_type'] = 'update'
                    entry['_index'] = target_index_name
                    entry['_id'] = material.material_id
                    entry['_type'] = 'doc'
                    entry['_source'] = {
                        "upsert": material.m_to_dict(include_defaults=False, partial="es"),
                        "doc_as_upsert": False,
                        "script": {
                            "source": "ctx._source.calculations.add(params.calc)",
                            "params": {
                                "calc": calc.m_to_dict(include_defaults=False, partial="es")
                            },
                        }
                    }
                    yield entry
    elif source == "es":
        s = elasticsearch_dsl.Search(index=config.elastic.index_name)
        filters = [elasticsearch_dsl.Q("term", encyclopedia__status="success")]
        if code:
            filters.append(elasticsearch_dsl.Q("terms", dft__code_name=code))
        query = elasticsearch_dsl.Q(
            "bool",
            filter=filters,
        )
        s = s.query(query)
        s = s.extra(**{
            "size": 0,
        })
        all_calcs = s.execute().hits.total
        print('indexing materials from %d calculations ...' % all_calcs)

        def elastic_updates():
            with utils.ETA(all_calcs, '   index %10d of %10d calcs, ETA %s', chunk_size) as eta:

                s = elasticsearch_dsl.Search(index=config.elastic.index_name)
                filters = [elasticsearch_dsl.Q("term", encyclopedia__status="success")]
                if code:
                    filters.append(elasticsearch_dsl.Q("terms", dft__code_name=code))
                query = elasticsearch_dsl.Q(
                    "bool",
                    filter=filters,
                )
                s = s.query(query)
                s = s.extra(**{
                    "size": chunk_size,
                })
                i_calc = 0
                for hit in s.scan():
                    i_calc += 1
                    if n is not None:
                        if i_calc > n:
                            return
                    eta.add()

                    material: Material = Material()
                    calc = Calculation()

                    # Check that all required information exists. If not, the
                    # calculation is skipped.
                    try:
                        material.material_id = hit.encyclopedia.material.material_id
                        material.material_type = hit.encyclopedia.material.material_type
                        material.formula = hit.encyclopedia.material.formula
                        material.formula_reduced = hit.encyclopedia.material.formula_reduced
                        material.species_and_counts = hit.encyclopedia.material.species_and_counts
                        material.species = hit.encyclopedia.material.species
                        calc.calc_id = hit.calc_id
                        calc.upload_id = hit.upload_id
                        calc.published = hit.published
                        calc.with_embargo = hit.with_embargo
                        calc.owners = [x.user_id for x in hit.owners]
                        idealized_structure = IdealizedStructure.m_from_dict(hit.encyclopedia.material.idealized_structure.to_dict())
                        calc.m_add_sub_section(Calculation.idealized_structure, idealized_structure)

                        method = Method.m_from_dict(hit.encyclopedia.method.to_dict())
                        method.program_name = hit.dft.code_name
                        method.program_version = hit.dft.code_version
                        method.basis_set = hit.dft.basis_set
                        calc.m_add_sub_section(Calculation.method, method)

                        workflow = Workflow()
                        workflow.workflow_type = hit.encyclopedia.calculation.calculation_type
                        calc.m_add_sub_section(Calculation.workflow, workflow)
                    except AttributeError:
                        continue

                    # Not all materials have a name
                    try:
                        material.material_name = hit.encyclopedia.material.material_name
                    except AttributeError:
                        pass

                    # Not all materials have a bulk section
                    try:
                        bulk = Bulk.m_from_dict(hit.encyclopedia.material.bulk)
                        material.m_add_sub_section(Material.bulk, bulk)
                    except AttributeError:
                        pass

                    # Properties may not exist at all
                    try:
                        enc_properties = hit.encyclopedia.properties
                    except AttributeError:
                        pass
                    else:
                        properties = Properties()

                        # Energies may not be present in all calculations
                        try:
                            energies = Energies.m_from_dict(enc_properties.energies.to_dict())
                            properties.m_add_sub_section(Properties.energies, energies)
                        except AttributeError:
                            pass

                        # Gather the boolean flags that indicate the presence of
                        # certain properties
                        try:
                            properties.has_electronic_dos = enc_properties.electronic_dos is not None
                        except AttributeError:
                            properties.has_electronic_dos = False
                        try:
                            properties.has_electronic_band_structure = enc_properties.electronic_band_structure is not None
                        except AttributeError:
                            properties.has_electronic_band_structure = False
                        try:
                            properties.has_thermodynamical_properties = enc_properties.thermodynamical_properties is not None
                        except AttributeError:
                            properties.has_thermodynamical_properties = False

                        # Not all materials have an atomic density
                        try:
                            properties.atomic_density = enc_properties.atomic_density
                        except AttributeError:
                            pass

                        # Not all materials have a mass density
                        try:
                            properties.mass_density = enc_properties.mass_density
                        except AttributeError:
                            pass

                        # Not all materials have band gaps
                        try:
                            properties.band_gap = enc_properties.band_gap
                        except AttributeError:
                            pass

                        # Not all materials have band gap type
                        try:
                            properties.band_gap_direct = enc_properties.band_gap_direct
                        except AttributeError:
                            pass

                        calc.m_add_sub_section(Calculation.properties, properties)

                    material.m_add_sub_section(Material.calculations, calc)

                    # Update entry that inserts the full material info if entry
                    # does not exists, otherwise only adds the calculation into
                    # the nested subdocument
                    entry = {}
                    entry['_op_type'] = 'update'
                    entry['_index'] = target_index_name
                    entry['_id'] = material.material_id
                    entry['_type'] = 'doc'
                    entry['_source'] = {
                        "upsert": material.m_to_dict(include_defaults=False, partial="es"),
                        "doc_as_upsert": False,
                        "script": {
                            "params": {
                                "calc": calc.m_to_dict(include_defaults=False, partial="es")
                            },
                        }
                    }
                    if in_place:
                        entry['_source']["script"]["source"] = "ctx._source.calculations.removeIf(x -> x.calc_id == params.calc.calc_id); ctx._source.calculations.add(params.calc)"
                    else:
                        entry['_source']["script"]["source"] = "ctx._source.calculations.add(params.calc)"
                    yield entry

    if dry:
        for _ in elastic_updates():
            pass
    else:
        # Create new index into which the data will be inserted. The old index will
        # keep working while the new index is being built
        material_document.init(index=target_index_name)

        if threads > 1:
            print('  use %d threads' % threads)
            for _ in elasticsearch.helpers.parallel_bulk(
                    infrastructure.elastic_client, elastic_updates(), chunk_size=chunk_size,
                    thread_count=threads):
                pass
        else:
            elasticsearch.helpers.bulk(
                infrastructure.elastic_client, elastic_updates())
            search.refresh()

        # Changes materials index alias to point to the new index and remove the
        # old index.
        if not in_place:
            new_index = elasticsearch_dsl.Index(target_index_name)
            new_index.put_alias(name=config.elastic.materials_index_name)
            old_index = elasticsearch_dsl.Index(old_index_name)
            old_index.delete()

    print('')
    print('indexing completed')


@admin.group(help='Generate scripts and commands for nomad operation.')
def ops():
    pass


@ops.group(help='Tools for managing the DOS similarity data.')
def similarity():
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
    proxy_set_header Host $host;

    location / {{
        proxy_pass http://app:8000;
    }}

    location ~ {1}\\/?(gui)?$ {{
        rewrite ^ {1}/gui/ permanent;
    }}

    location {1}/gui/ {{
        proxy_intercept_errors on;
        error_page 404 = @redirect_to_index;
        proxy_pass http://app:8000;
    }}

    location @redirect_to_index {{
        rewrite ^ {1}/gui/index.html break;
        proxy_pass http://app:8000;
    }}

    location ~ \\/gui\\/(service-worker\\.js|meta\\.json)$ {{
        add_header Last-Modified $date_gmt;
        add_header Cache-Control 'no-store, no-cache, must-revalidate, proxy-revalidate, max-age=0';
        if_modified_since off;
        expires off;
        etag off;
        proxy_pass http://app:8000;
    }}

    location ~ \\/api\\/uploads\\/?$ {{
        client_max_body_size 35g;
        proxy_request_buffering off;
        proxy_pass http://app:8000;
    }}

    location ~ \\/api\\/(raw|archive) {{
        proxy_buffering off;
        proxy_pass http://app:8000;
    }}

    location ~ \\/api\\/mirror {{
        proxy_buffering off;
        proxy_read_timeout 600;
        proxy_pass http://app:8000;
    }}
}}'''.format(prefix))


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
'''.format(prefix, host, port))  # type: ignore


@ops.command(help='Updates the AFLOW prototype information using the latest online version and writes the results to a python module in the given FILEPATH.')
@click.argument('FILEPATH', nargs=1, type=str)
@click.option('--matches-only', is_flag=True, help='Only update the match information that depends on the symmetry analysis settings. Will not perform and online update.')
@click.pass_context
def prototypes_update(ctx, filepath, matches_only):
    from nomad.cli.admin import prototypes
    prototypes.update_prototypes(ctx, filepath, matches_only)


@ops.command(help='Updates the springer database in nomad.config.normalize.springer_db_path.')
@click.option('--max-n-query', default=10, type=int, help='Number of unsuccessful springer request before returning an error. Default is 10.')
@click.option('--retry-time', default=120, type=int, help='Time in seconds to retry after unsuccessful request. Default is 120.')
def springer_update(max_n_query, retry_time):
    from nomad.cli.admin import springer
    springer.update_springer(max_n_query, retry_time)


@similarity.command(help='Updates the msgpack file containing the similarity information.')
@click.option('--dir', "-d", "input_dir", type=str, help='Path of the folder containing the raw similarity information files')
@click.option('--out', "-o", type=str, help='Path of the output msgpack file.')
@click.option('--verbose', is_flag=True, help='Enable verbose output.')
def update(input_dir, out, verbose):
    from nomad.cli.admin import similarity
    similarity.update(input_dir, out, verbose)


@similarity.command(help='Ingests the given similarity information from an msgpack file into MongoDB.')
@click.option('--in', "-i", "input_path", type=str, help='Path of the ingested msgpack file.')
@click.option('--batch_size', type=int, default=10000, help='Batch size for MongoDB bulk ingestion.')
@click.option('--verbose', is_flag=True, help='Enable verbose output.')
def ingest(input_path, batch_size, verbose):
    from nomad.cli.admin import similarity
    similarity.ingest(input_path, batch_size, verbose)
