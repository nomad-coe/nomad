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

from typing import Generator, Tuple
import pytest
from datetime import datetime
import os.path
import json
import re
import shutil

from nomad import utils, infrastructure, config, datamodel
from nomad.files import UploadFiles, StagingUploadFiles, PublicUploadFiles
from nomad.processing import Upload, Calc
from nomad.processing.base import task as task_decorator, FAILURE, SUCCESS

from tests.test_search import assert_search_upload
from tests.test_files import assert_upload_files


def test_send_mail(mails, monkeypatch):
    infrastructure.send_mail('test name', 'test@email.de', 'test message', 'subject')

    for message in mails.messages:
        assert re.search(r'test message', message.data.decode('utf-8')) is not None


@pytest.fixture(scope='function', autouse=True)
def mongo_forall(mongo):
    pass


@pytest.fixture
def uploaded_id_with_warning(raw_files) -> Generator[Tuple[str, str], None, None]:
    example_file = 'tests/data/proc/examples_with_warning_template.zip'
    example_upload_id = os.path.basename(example_file).replace('.zip', '')

    yield example_upload_id, example_file


def run_processing(uploaded: Tuple[str, str], test_user) -> Upload:
    uploaded_id, uploaded_path = uploaded
    upload = Upload.create(
        upload_id=uploaded_id, user=test_user, upload_path=uploaded_path)
    upload.upload_time = datetime.utcnow()

    assert upload.tasks_status == 'RUNNING'
    assert upload.current_task == 'uploading'

    upload.process_upload()  # pylint: disable=E1101
    upload.block_until_complete(interval=.01)

    return upload


def assert_processing(upload: Upload, published: bool = False):
    assert not upload.tasks_running
    assert upload.current_task == 'cleanup'
    assert upload.upload_id is not None
    assert len(upload.errors) == 0
    assert upload.tasks_status == SUCCESS
    assert upload.joined

    upload_files = UploadFiles.get(upload.upload_id, is_authorized=lambda: True)
    if published:
        assert isinstance(upload_files, PublicUploadFiles)
    else:
        assert isinstance(upload_files, StagingUploadFiles)

    for calc in Calc.objects(upload_id=upload.upload_id):
        assert calc.parser is not None
        assert calc.mainfile is not None
        assert calc.tasks_status == SUCCESS
        assert calc.metadata['published'] == published

        with upload_files.archive_file(calc.calc_id) as archive_json:
            archive = json.load(archive_json)
        assert 'section_run' in archive
        assert 'section_entry_info' in archive

        with upload_files.archive_log_file(calc.calc_id, 'rt') as f:
            assert 'a test' in f.read()
        assert len(calc.errors) == 0

        with upload_files.raw_file(calc.mainfile) as f:
            f.read()

        for path in calc.metadata['files']:
            with upload_files.raw_file(path) as f:
                f.read()

        assert upload.get_calc(calc.calc_id).metadata is not None


def test_processing(processed, no_warn, mails, monkeypatch):
    assert_processing(processed)

    assert len(mails.messages) == 1
    assert re.search(r'Processing completed', mails.messages[0].data.decode('utf-8')) is not None


def test_processing_with_large_dir(test_user, proc_infra):
    upload_path = 'tests/data/proc/examples_large_dir.zip'
    upload_id = upload_path[:-4]
    upload = run_processing((upload_id, upload_path), test_user)
    for calc in upload.calcs:
        assert len(calc.warnings) == 1


def test_publish(non_empty_processed: Upload, no_warn, example_user_metadata, monkeypatch):
    processed = non_empty_processed
    processed.compress_and_set_metadata(example_user_metadata)

    additional_keys = ['with_embargo']

    processed.publish_upload()
    try:
        processed.block_until_complete(interval=.01)
    except Exception:
        pass

    upload = processed.to_upload_with_metadata(example_user_metadata)

    assert_upload_files(upload, PublicUploadFiles, published=True)
    assert_search_upload(upload, additional_keys, published=True)

    assert_processing(Upload.get(upload.upload_id, include_published=True), published=True)


def test_republish(non_empty_processed: Upload, no_warn, example_user_metadata, monkeypatch):
    processed = non_empty_processed
    processed.compress_and_set_metadata(example_user_metadata)

    additional_keys = ['with_embargo']

    processed.publish_upload()
    processed.block_until_complete(interval=.01)
    assert Upload.get('examples_template') is not None

    processed.publish_upload()
    processed.block_until_complete(interval=.01)

    upload = processed.to_upload_with_metadata(example_user_metadata)

    assert_upload_files(upload, PublicUploadFiles, published=True)
    assert_search_upload(upload, additional_keys, published=True)


def test_publish_failed(
        non_empty_uploaded: Tuple[str, str], example_user_metadata, test_user,
        monkeypatch, proc_infra):

    mock_failure(Calc, 'parsing', monkeypatch)

    processed = run_processing(non_empty_uploaded, test_user)
    processed.compress_and_set_metadata(example_user_metadata)

    additional_keys = ['with_embargo']

    processed.publish_upload()
    try:
        processed.block_until_complete(interval=.01)
    except Exception:
        pass

    upload = processed.to_upload_with_metadata(example_user_metadata)

    assert_search_upload(upload, additional_keys, published=True, processed=False)


@pytest.mark.timeout(config.tests.default_timeout)
def test_processing_with_warning(proc_infra, test_user, with_warn):
    example_file = 'tests/data/proc/examples_with_warning_template.zip'
    example_upload_id = os.path.basename(example_file).replace('.zip', '')

    upload = run_processing((example_upload_id, example_file), test_user)
    assert_processing(upload)


@pytest.mark.timeout(config.tests.default_timeout)
def test_process_non_existing(proc_infra, test_user, with_error):
    upload = run_processing(('__does_not_exist', '__does_not_exist'), test_user)

    assert not upload.tasks_running
    assert upload.current_task == 'extracting'
    assert upload.tasks_status == FAILURE
    assert len(upload.errors) > 0


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.mark.parametrize('with_failure', [None, 'before', 'after'])
def test_re_processing(published: Upload, example_user_metadata, monkeypatch, with_failure):
    if with_failure == 'before':
        calc = published.all_calcs(0, 1).first()
        calc.tasks_status = FAILURE
        calc.errors = ['example error']
        calc.save()
        assert published.failed_calcs > 0

    assert published.published
    assert published.upload_files.to_staging_upload_files() is None

    with_failure = with_failure == 'after'

    old_upload_time = published.last_update
    first_calc = published.all_calcs(0, 1).first()
    old_calc_time = first_calc.metadata['last_processing']

    if not with_failure:
        with published.upload_files.archive_log_file(first_calc.calc_id) as f:
            old_log_line = f.readline()
    old_archive_files = list(
        archive_file
        for archive_file in os.listdir(published.upload_files.os_path)
        if 'archive' in archive_file)
    for archive_file in old_archive_files:
        with open(published.upload_files.join_file(archive_file).os_path, 'wt') as f:
            f.write('')

    if with_failure:
        raw_files = 'tests/data/proc/examples_template_unparsable.zip'
    else:
        raw_files = 'tests/data/proc/examples_template_different_atoms.zip'
    shutil.copyfile(
        raw_files, published.upload_files.join_file('raw-restricted.plain.zip').os_path)

    upload = published.to_upload_with_metadata(example_user_metadata)

    # reprocess
    monkeypatch.setattr('nomad.config.version', 're_process_test_version')
    monkeypatch.setattr('nomad.config.commit', 're_process_test_commit')
    published.reset()
    published.re_process_upload()
    try:
        published.block_until_complete(interval=.01)
    except Exception:
        pass

    published.reload()
    first_calc.reload()

    # assert new process time
    assert published.last_update > old_upload_time
    assert first_calc.metadata['last_processing'] > old_calc_time

    # assert new process version
    assert first_calc.metadata['nomad_version'] == 're_process_test_version'
    assert first_calc.metadata['nomad_commit'] == 're_process_test_commit'

    # assert changed archive files
    if not with_failure:
        for archive_file in old_archive_files:
            assert os.path.getsize(published.upload_files.join_file(archive_file).os_path) > 0

    # assert changed archive log files
    if not with_failure:
        with published.upload_files.archive_log_file(first_calc.calc_id) as f:
            new_log_line = f.readline()
        assert old_log_line != new_log_line

    # assert maintained user metadata (mongo+es)
    assert_upload_files(upload, PublicUploadFiles, published=True)
    assert_search_upload(upload, published=True)
    if not with_failure:
        assert_processing(Upload.get(upload.upload_id, include_published=True), published=True)

    # assert changed calc metadata (mongo)
    if not with_failure:
        assert first_calc.metadata['atoms'][0] == 'H'
    else:
        assert first_calc.metadata['atoms'][0] == 'Br'


@pytest.mark.timeout(config.tests.default_timeout)
@pytest.mark.parametrize('with_failure', [None, 'before', 'after'])
def test_re_pack(published: Upload, example_user_metadata, monkeypatch, with_failure):
    upload_id = published.upload_id
    calc = Calc.objects(upload_id=upload_id).first()
    assert calc.metadata['with_embargo']
    calc.metadata['with_embargo'] = False
    calc.save()

    published.re_pack()
    try:
        published.block_until_complete(interval=.01)
    except Exception:
        pass

    upload_files = PublicUploadFiles(upload_id)
    for raw_file in upload_files.raw_file_manifest():
        with upload_files.raw_file(raw_file) as f:
            f.read()
    for calc in Calc.objects(upload_id=upload_id):
        with upload_files.archive_file(calc.calc_id) as f:
            f.read()


def mock_failure(cls, task, monkeypatch):
    def mock(self):
        raise Exception('fail for test')

    mock.__name__ = task
    mock = task_decorator(mock)

    monkeypatch.setattr('nomad.processing.data.%s.%s' % (cls.__name__, task), mock)


@pytest.mark.parametrize('task', ['extracting', 'parse_all', 'cleanup', 'parsing'])
@pytest.mark.timeout(config.tests.default_timeout)
def test_task_failure(monkeypatch, uploaded, task, proc_infra, test_user, with_error):
    # mock the task method to through exceptions
    if hasattr(Upload, task):
        cls = Upload
    elif hasattr(Calc, task):
        cls = Calc
    else:
        assert False

    mock_failure(cls, task, monkeypatch)

    # run the test
    upload = run_processing(uploaded, test_user)

    assert not upload.tasks_running

    if task != 'parsing':
        assert upload.tasks_status == FAILURE
        assert upload.current_task == task
        assert len(upload.errors) > 0
    else:
        # there is an empty example with no calcs, even if past parsing_all task
        utils.get_logger(__name__).error('fake')
        if upload.total_calcs > 0:  # pylint: disable=E1101
            assert upload.tasks_status == SUCCESS
            assert upload.current_task == 'cleanup'
            assert len(upload.errors) == 0
            for calc in upload.all_calcs(0, 100):  # pylint: disable=E1101
                assert calc.tasks_status == FAILURE
                assert calc.current_task == 'parsing'
                assert len(calc.errors) > 0

# TODO timeout
# consume_ram, segfault, and exit are not testable with the celery test worker
@pytest.mark.parametrize('failure', ['exception'])
def test_malicious_parser_task_failure(proc_infra, failure, test_user):
    example_file = 'tests/data/proc/chaos_%s.zip' % failure
    example_upload_id = os.path.basename(example_file).replace('.zip', '')

    upload = run_processing((example_upload_id, example_file), test_user)

    assert not upload.tasks_running
    assert upload.current_task == 'cleanup'
    assert len(upload.errors) == 0
    assert upload.tasks_status == SUCCESS

    calcs = Calc.objects(upload_id=upload.upload_id)
    assert calcs.count() == 1
    calc = next(calcs)
    assert not calc.tasks_running
    assert calc.tasks_status == FAILURE
    assert len(calc.errors) == 1


def test_ems_data(proc_infra, test_user, monkeypatch):
    monkeypatch.setattr('nomad.config.domain', 'EMS')
    monkeypatch.setattr('nomad.datamodel.Domain.instance', datamodel.Domain.instances['EMS'])
    monkeypatch.setattr('nomad.datamodel.CalcWithMetadata', datamodel.Domain.instance.domain_entry_class)

    upload = run_processing(('test_ems_upload', 'tests/data/proc/example_ems.zip'), test_user)

    additional_keys = ['method', 'experiment_location', 'experiment_time', 'formula', 'chemical']
    assert upload.total_calcs == 1
    assert len(upload.calcs) == 1

    upload_with_metadata = upload.to_upload_with_metadata()
    assert_upload_files(upload_with_metadata, StagingUploadFiles, published=False)
    assert_search_upload(upload_with_metadata, additional_keys, published=False)
