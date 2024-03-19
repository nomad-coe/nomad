import math
import os
from datetime import datetime, timezone
from typing import List, Tuple

import pytest

from nomad import bundles, datamodel, processing, utils
from nomad.archive import (
    read_archive,
    to_json,
    write_archive,
    write_partial_archive_to_mongo,
)
from nomad.config import config
from nomad.datamodel import EntryArchive, OptimadeEntry, User
from nomad.datamodel.datamodel import SearchableQuantity
from nomad.metainfo.elasticsearch_extension import schema_separator
from nomad.processing import ProcessStatus
from nomad.utils.exampledata import ExampleData
from tests.config import python_schema_name, yaml_schema_name, yaml_schema_root
from tests.normalizing.conftest import run_normalize
from tests.parsing import test_parsing
from tests.processing import test_data as test_processing
from tests.test_files import empty_file, example_file_vasp_with_binary
from tests.utils import (
    build_url,
    create_template_upload_file,
    set_upload_entry_metadata,
)


@pytest.fixture(scope='session')
def example_mainfile() -> Tuple[str, str]:
    return ('parsers/template', 'tests/data/templates/template.json')


@pytest.fixture(scope='function', params=['empty_file', 'example_file'])
def example_upload(request, tmp) -> str:
    if request.param == 'empty_file':
        return create_template_upload_file(tmp, mainfiles=[], auxfiles=0)

    return create_template_upload_file(
        tmp, mainfiles=['tests/data/proc/templates/template.json']
    )


@pytest.fixture(scope='function')
def non_empty_example_upload(tmp):
    return create_template_upload_file(
        tmp, mainfiles=['tests/data/proc/templates/template.json']
    )


@pytest.fixture(scope='session')
def non_empty_example_upload_vasp_with_binary():
    return example_file_vasp_with_binary


@pytest.fixture(scope='session')
def empty_upload():
    return empty_file


@pytest.fixture(scope='module')
def example_user_metadata(other_test_user, test_user) -> dict:
    return {
        'comment': 'test comment',
        'references': ['http://external.ref/one', 'http://external.ref/two'],
        'entry_coauthors': [other_test_user.user_id],
        '_pid': '256',
        'external_id': 'external_test_id',
    }


@pytest.fixture(scope='module')
def internal_example_user_metadata(example_user_metadata) -> dict:
    return {
        key[1:] if key[0] == '_' else key: value
        for key, value in example_user_metadata.items()
    }


@pytest.fixture(scope='session')
def parsed(example_mainfile: Tuple[str, str]) -> EntryArchive:
    """Provides a parsed entry in the form of an EntryArchive."""
    parser, mainfile = example_mainfile
    return test_parsing.run_singular_parser(parser, mainfile)


@pytest.fixture(scope='session')
def parsed_ems() -> EntryArchive:
    """Provides a parsed experiment in the form of a EntryArchive."""
    return test_parsing.run_singular_parser(
        'parsers/eels', 'tests/data/parsers/eels.json'
    )


@pytest.fixture(scope='session')
def normalized(parsed: EntryArchive) -> EntryArchive:
    """Provides a normalized entry in the form of a EntryArchive."""
    return run_normalize(parsed)


@pytest.fixture(scope='function')
def uploaded(example_upload: str, raw_files_function) -> Tuple[str, str]:
    """
    Provides a uploaded with uploaded example file and gives the upload_id.
    Clears files after test.
    """
    example_upload_id = os.path.basename(example_upload).replace('.zip', '')
    return example_upload_id, example_upload


@pytest.fixture(scope='function')
def non_empty_uploaded(
    non_empty_example_upload: str, raw_files_function
) -> Tuple[str, str]:
    example_upload_id = os.path.basename(non_empty_example_upload).replace('.zip', '')
    return example_upload_id, non_empty_example_upload


@pytest.fixture(scope='function')
def oasis_publishable_upload(
    api_v1,
    proc_infra,
    non_empty_processed: processing.Upload,
    internal_example_user_metadata,
    monkeypatch,
    test_user,
):
    """
    Creates a published upload which can be used with Upload.publish_externally. Some monkeypatching
    is done which replaces IDs when importing.
    """
    # Create a published upload
    set_upload_entry_metadata(non_empty_processed, internal_example_user_metadata)
    non_empty_processed.publish_upload()
    non_empty_processed.block_until_complete(interval=0.01)

    suffix = '_2'  # Will be added to all IDs in the mirrored upload
    upload_id = non_empty_processed.upload_id

    # Do some tricks to add suffix to the ID fields
    old_bundle_importer_open = bundles.BundleImporter.open

    def new_bundle_importer_open(self, *args, **kwargs):
        old_bundle_importer_open(self, *args, **kwargs)
        # Change the id's in the bundle_info dict behind the scenes when loading the bundle
        bundle_info = self.bundle_info
        bundle_info['upload_id'] += suffix
        bundle_info['upload']['_id'] += suffix
        for entry_dict in bundle_info['entries']:
            entry_dict['_id'] = utils.generate_entry_id(
                upload_id + suffix,
                entry_dict['mainfile'],
                entry_dict.get('mainfile_key'),
            )
            entry_dict['upload_id'] += suffix

    old_bundle_import_files = bundles.BundleImporter._import_files

    def new_bundle_import_files(self, *args, **kwargs):
        old_bundle_import_files(self, *args, **kwargs)
        # Overwrite the archive files with files containing the updated IDs
        archive_path = self.upload_files.os_path
        for file_name in os.listdir(archive_path):
            if file_name.endswith('.msg'):
                full_path = os.path.join(archive_path, file_name)
                new_data = []
                with read_archive(full_path) as data:
                    for entry_id in data.keys():
                        archive_dict = to_json(data[entry_id])
                        section_metadata = archive_dict['metadata']
                        section_metadata['upload_id'] += suffix
                        new_entry_id = utils.generate_entry_id(
                            section_metadata['upload_id'],
                            section_metadata['mainfile'],
                            section_metadata.get('mainfile_key'),
                        )
                        section_metadata['entry_id'] = new_entry_id
                        new_data.append((new_entry_id, archive_dict))
                write_archive(full_path, len(new_data), new_data)

    monkeypatch.setattr('nomad.bundles.BundleImporter.open', new_bundle_importer_open)
    monkeypatch.setattr(
        'nomad.bundles.BundleImporter._import_files', new_bundle_import_files
    )

    # Further monkey patching
    def new_post(url, data, params={}, **kwargs):
        return api_v1.post(
            build_url(url.lstrip('/api/v1/'), params), data=data.read(), **kwargs
        )

    monkeypatch.setattr('requests.post', new_post)
    monkeypatch.setattr('nomad.config.oasis.is_oasis', True)
    monkeypatch.setattr('nomad.config.keycloak.username', test_user.username)

    monkeypatch.setattr('nomad.config.oasis.central_nomad_deployment_url', '/api')

    # create a dataset to also test this aspect of oasis uploads
    entry = non_empty_processed.successful_entries[0]
    datamodel.Dataset(
        dataset_id='dataset_id', dataset_name='dataset_name', user_id=test_user.user_id
    ).a_mongo.save()
    entry.datasets = ['dataset_id']
    entry.save()
    return non_empty_processed.upload_id, suffix


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def processed(
    uploaded: Tuple[str, str], test_user: User, proc_infra, mails
) -> processing.Upload:
    """
    Provides a processed upload. Upload was uploaded with test_user.
    """
    return test_processing.run_processing(uploaded, test_user)


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def processeds(
    non_empty_example_upload: str, test_user: User, proc_infra
) -> List[processing.Upload]:
    result: List[processing.Upload] = []
    for i in range(2):
        upload_id = '%s_%d' % (
            os.path.basename(non_empty_example_upload).replace('.zip', ''),
            i,
        )
        result.append(
            test_processing.run_processing(
                (upload_id, non_empty_example_upload), test_user
            )
        )

    return result


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def non_empty_processed(
    non_empty_uploaded: Tuple[str, str], test_user: User, proc_infra
) -> processing.Upload:
    """
    Provides a processed upload. Upload was uploaded with test_user.
    """
    return test_processing.run_processing(non_empty_uploaded, test_user)


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def published(
    non_empty_processed: processing.Upload, internal_example_user_metadata
) -> processing.Upload:
    """
    Provides a processed published upload. Upload was uploaded with test_user and is embargoed.
    """
    set_upload_entry_metadata(non_empty_processed, internal_example_user_metadata)
    non_empty_processed.publish_upload(embargo_length=12)
    try:
        non_empty_processed.block_until_complete(interval=0.01)
    except Exception:
        pass

    return non_empty_processed


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.fixture(scope='function')
def published_wo_user_metadata(
    non_empty_processed: processing.Upload,
) -> processing.Upload:
    """
    Provides a processed upload. Upload was uploaded with test_user.
    """
    non_empty_processed.publish_upload()
    try:
        non_empty_processed.block_until_complete(interval=0.01)
    except Exception:
        pass

    return non_empty_processed


@pytest.fixture(scope='module')
def example_data(
    elastic_module,
    raw_files_module,
    mongo_module,
    test_user,
    other_test_user,
    normalized,
):
    """
    Provides a couple of uploads and entries including metadata, raw-data, and
    archive files.

    id_embargo:
        1 entry, 1 material, published with embargo
    id_embargo_w_coauthor:
        1 entry, 1 material, published with embargo and coauthor
    id_embargo_w_reviewer:
        1 entry, 1 material, published with embargo and reviewer
    id_unpublished:
        1 entry, 1 material, unpublished
    id_unpublished_w_coauthor:
        1 entry, 1 material, unpublished with coauthor
    id_unpublished_w_reviewer:
        1 entry, 1 material, unpublished with reviewer
    id_published:
        23 entries, 6 materials published without embargo
        partial archive exists only for id_01
        raw files and archive file for id_02 are missing
        id_10, id_11 reside in the same directory
    id_child_entries:
        1 parent entry and 2 child entries from one mainfile, 1 material, unpublished
    id_processing:
        unpublished upload without any entries, in status processing
    id_empty:
        unpublished upload without any entries
    """
    data = ExampleData(main_author=test_user)

    # 6 uploads with different combinations of main_type and sub_type
    for main_type in ('embargo', 'unpublished'):
        for sub_type in ('', 'w_coauthor', 'w_reviewer'):
            upload_id = 'id_' + main_type + ('_' if sub_type else '') + sub_type
            if main_type == 'embargo':
                published = True
                embargo_length = 12
                upload_name = 'name_' + upload_id[3:]
            else:
                published = False
                embargo_length = 0
                upload_name = None
            entry_id = upload_id + '_1'
            coauthors = [other_test_user.user_id] if sub_type == 'w_coauthor' else None
            reviewers = [other_test_user.user_id] if sub_type == 'w_reviewer' else None
            data.create_upload(
                upload_id=upload_id,
                upload_name=upload_name,
                coauthors=coauthors,
                reviewers=reviewers,
                published=published,
                embargo_length=embargo_length,
            )
            data.create_entry(
                upload_id=upload_id,
                entry_id=entry_id,
                material_id=upload_id,
                mainfile=f'test_content/{entry_id}/mainfile.json',
            )

    # one upload with 23 entries, published, no embargo
    data.create_upload(
        upload_id='id_published', upload_name='name_published', published=True
    )
    for i in range(1, 24):
        entry_id = 'id_%02d' % i
        material_id = 'id_%02d' % (int(math.floor(i / 4)) + 1)
        mainfile = 'test_content/subdir/test_entry_%02d/mainfile.json' % i
        kwargs = dict(
            optimade=OptimadeEntry(nelements=2, elements=['H', 'O']),
        )
        if i == 11:
            mainfile = 'test_content/subdir/test_entry_10/mainfile_11.json'
        if i == 1:
            kwargs['pid'] = '123'
        data.create_entry(
            upload_id='id_published',
            entry_id=entry_id,
            material_id=material_id,
            mainfile=mainfile,
            **kwargs,
        )

        if i == 1:
            archive = data.archives[entry_id]
            write_partial_archive_to_mongo(archive)

    # 3 entries from one mainfile, 1 material, unpublished
    upload_id = 'id_child_entries'
    data.create_upload(
        upload_id=upload_id, upload_name='name_child_entries', published=False
    )
    for mainfile_key in (None, 'child1', 'child2'):
        data.create_entry(
            upload_id=upload_id,
            entry_id=upload_id + '_' + (mainfile_key or 'main'),
            material_id=upload_id,
            mainfile=f'test_content/mainfile_w_children.json',
            mainfile_key=mainfile_key,
        )

    # one upload, no entries, still processing
    data.create_upload(
        upload_id='id_processing', published=False, process_status=ProcessStatus.RUNNING
    )

    # one upload, no entries, unpublished
    data.create_upload(upload_id='id_empty', published=False)

    data.save(with_files=False)
    del data.archives['id_02']
    data.save(with_files=True, with_es=False, with_mongo=False)

    # yield

    # # The data is deleted after fixture goes out of scope
    # data.delete()


@pytest.fixture(scope='function')
def example_data_schema_python(
    elastic_module, raw_files_module, mongo_module, test_user, normalized
):
    """
    Contains entries that store data using a python schema.
    """
    data = ExampleData(main_author=test_user)
    upload_id = 'id_plugin_schema_published'
    date_value = datetime.now(timezone.utc)

    data.create_upload(upload_id=upload_id, upload_name=upload_id, published=True)
    for i in range(0, 15):
        data.create_entry(
            upload_id=upload_id,
            entry_id=f'test_entry_{i}',
            mainfile=f'test_content/test.archive{i}.json',
            search_quantities=[
                SearchableQuantity(
                    id=f'data.name{schema_separator}{python_schema_name}',
                    definition='nomadschemaexample.schema.MySchema.name',
                    path_archive='data.name',
                    str_value=f'test{i}',
                ),
                SearchableQuantity(
                    id=f'data.valid{schema_separator}{python_schema_name}',
                    definition='nomadschemaexample.schema.MySchema.valid',
                    path_archive='data.valid',
                    bool_value=i % 2 == 0,
                ),
                SearchableQuantity(
                    id=f'data.message{schema_separator}{python_schema_name}',
                    definition='nomadschemaexample.schema.MySchema.message',
                    path_archive='data.message',
                    str_value='A' if i % 2 == 0 else 'B',
                ),
                SearchableQuantity(
                    id=f'data.frequency{schema_separator}{python_schema_name}',
                    definition='nomadschemaexample.schema.MySchema.frequency',
                    path_archive='data.frequency',
                    float_value=i + 0.5,
                ),
                SearchableQuantity(
                    id=f'data.count{schema_separator}{python_schema_name}',
                    definition='nomadschemaexample.schema.MySchema.count',
                    path_archive='data.count',
                    int_value=i,
                ),
                SearchableQuantity(
                    id=f'data.timestamp{schema_separator}{python_schema_name}',
                    definition='nomadschemaexample.schema.MySchema.timestamp',
                    path_archive='data.timestamp',
                    datetime_value=date_value,
                ),
                SearchableQuantity(
                    id=f'data.child.name{schema_separator}{python_schema_name}',
                    definition='nomadschemaexample.schema.MySection.name',
                    path_archive='data.child.name',
                    str_value=f'test_child{i}',
                ),
                SearchableQuantity(
                    id=f'data.child_repeating.name{schema_separator}{python_schema_name}',
                    definition='nomadschemaexample.schema.MySection.name',
                    path_archive='data.child_repeating.0.name',
                    str_value=f'test_child_repeating{i}',
                ),
            ],
        )
    data.save(with_files=False)

    yield

    # The data is deleted after fixture goes out of scope
    data.delete()


@pytest.fixture(scope='function')
def example_data_nexus(
    elastic_module, raw_files_module, mongo_module, test_user, normalized
):
    """
    Contains entries that store data using a python schema.
    """
    data = ExampleData(main_author=test_user)
    upload_id = 'id_nexus_published'

    data.create_upload(upload_id=upload_id, upload_name=upload_id, published=True)
    data.create_entry(
        upload_id=upload_id,
        entry_id=f'test_entry_nexus',
        mainfile=f'test_content/test.archive.json',
        search_quantities=[
            SearchableQuantity(
                id=f'nexus.NXiv_temp.ENTRY.DATA.temperature__field',
                definition='nexus.NXiv_temp.ENTRY.DATA.temperature__field',
                path_archive='nexus.NXiv_temp.ENTRY.0.DATA.0.temperature__field',
                float_value=273.15,
            ),
            SearchableQuantity(
                id=f'nexus.NXiv_temp.ENTRY.definition__field',
                definition='nexus.NXiv_temp.NXentry.definition__field',
                path_archive='nexus.NXiv_temp.ENTRY.0.definition__field',
                str_value='NXiv_temp',
            ),
        ],
    )
    data.save(with_files=False)

    yield

    # The data is deleted after fixture goes out of scope
    data.delete()


@pytest.fixture(scope='function')
def example_data_schema_yaml(
    elastic_module,
    raw_files_function,
    no_warn,
    raw_files_module,
    mongo_module,
    test_user,
    normalized,
):
    """
    Contains entries that store data using a python schema.
    """
    data = ExampleData(main_author=test_user)
    upload_id = 'id_plugin_schema_published'
    date_value = datetime.now(timezone.utc)

    data.create_upload(upload_id=upload_id, upload_name=upload_id, published=True)
    for i in range(0, 15):
        data.create_entry(
            upload_id=upload_id,
            entry_id=f'test_entry_{i}',
            mainfile=f'test_content/test.archive{i}.json',
            search_quantities=[
                SearchableQuantity(
                    id=f'data.name{schema_separator}{yaml_schema_name}',
                    definition=f'{yaml_schema_root}/quantities/0',
                    path_archive='data.name',
                    str_value=f'test{i}',
                ),
                SearchableQuantity(
                    id=f'data.valid{schema_separator}{yaml_schema_name}',
                    definition=f'{yaml_schema_root}/quantities/1',
                    path_archive='data.valid',
                    bool_value=i % 2 == 0,
                ),
                SearchableQuantity(
                    id=f'data.message{schema_separator}{yaml_schema_name}',
                    definition=f'{yaml_schema_root}/quantities/1',
                    path_archive='data.message',
                    str_value='A' if i % 2 == 0 else 'B',
                ),
                SearchableQuantity(
                    id=f'data.frequency{schema_separator}{yaml_schema_name}',
                    definition=f'{yaml_schema_root}/quantities/5',
                    path_archive='data.frequency',
                    float_value=i + 0.5,
                ),
                SearchableQuantity(
                    id=f'data.count{schema_separator}{yaml_schema_name}',
                    definition=f'{yaml_schema_root}/quantities/4',
                    path_archive='data.count',
                    int_value=i,
                ),
                SearchableQuantity(
                    id=f'data.timestamp{schema_separator}{yaml_schema_name}',
                    definition=f'{yaml_schema_root}/quantities/1',
                    path_archive='data.timestamp',
                    datetime_value=date_value,
                ),
                SearchableQuantity(
                    id=f'data.child.name{schema_separator}{yaml_schema_name}',
                    definition=f'{yaml_schema_root}/section_definitions/0/quantities/0',
                    path_archive='data.child.name',
                    str_value=f'test_child{i}',
                ),
                SearchableQuantity(
                    id=f'data.child_repeating.name{schema_separator}{yaml_schema_name}',
                    definition=f'{yaml_schema_root}/section_definitions/0/quantities/0',
                    path_archive='data.child_repeating.0.name',
                    str_value=f'test_child_repeating{i}',
                ),
            ],
        )
    data.save(with_files=False)

    yield

    # The data is deleted after fixture goes out of scope
    data.delete()


@pytest.fixture(scope='function')
def example_data_writeable(mongo_function, test_user, normalized):
    data = ExampleData(main_author=test_user)

    # one upload with one entry, published
    data.create_upload(upload_id='id_published_w', published=True, embargo_length=12)
    data.create_entry(
        upload_id='id_published_w',
        entry_id='id_published_w_entry',
        mainfile='test_content/test_embargo_entry/mainfile.json',
    )

    # one upload with one entry, unpublished
    data.create_upload(upload_id='id_unpublished_w', published=False, embargo_length=12)
    data.create_entry(
        upload_id='id_unpublished_w',
        entry_id='id_unpublished_w_entry',
        mainfile='test_content/test_embargo_entry/mainfile.json',
    )

    # one upload, no entries, running a blocking processing
    data.create_upload(
        upload_id='id_processing_w',
        published=False,
        process_status=ProcessStatus.RUNNING,
        current_process='publish_upload',
    )

    # one upload, no entries, unpublished
    data.create_upload(upload_id='id_empty_w', published=False)

    data.save()

    yield

    data.delete()


@pytest.fixture(scope='function')
def example_datasets(mongo_function, test_user, other_test_user):
    dataset_specs = (
        ('test_dataset_1', test_user, None),
        ('test_dataset_2', test_user, 'test_doi_2'),
        ('test_dataset_3', other_test_user, None),
    )
    datasets = []
    for dataset_name, user, doi in dataset_specs:
        now = datetime.utcnow()
        dataset = datamodel.Dataset(
            dataset_id=utils.create_uuid(),
            dataset_name=dataset_name,
            doi=doi,
            user_id=user.user_id,
            dataset_create_time=now,
            dataset_modified_time=now,
            dataset_type='owned',
        )
        dataset.a_mongo.create()
        datasets.append(dataset)

    yield datasets

    while datasets:
        datasets.pop().a_mongo.delete()
