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

import pytest
import json
import shutil
import os
import os.path
from zipfile import ZipFile

from nomad.files import Objects, ObjectFile, ArchiveFile, UploadFile, ArchiveLogFile, \
    BaggedDataContainer, ZippedDataContainer
from nomad import config

# example_file uses an artificial parser for faster test execution, can also be
# changed to examples_vasp.zip for using vasp parser
example_file = 'tests/data/proc/examples_template.zip'
empty_file = 'tests/data/proc/empty.zip'

example_bucket = 'test_bucket'
example_data = dict(test_key='test_value')


@pytest.fixture(scope='function')
def clear_files():
    """ Utility fixture that removes all files from files and tmp after test. """
    try:
        yield
    finally:
        try:
            shutil.rmtree(config.fs.objects)
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(config.fs.tmp)
        except FileNotFoundError:
            pass


class TestObjects:
    @pytest.fixture()
    def existing_example_file(self, clear_files):
        with ObjectFile(example_bucket, 'example_file', ext='json').open(mode='wt') as out:
            json.dump(example_data, out)

        yield 'example_file', 'json'

    def test_size(self, existing_example_file):
        name, ext = existing_example_file
        assert ObjectFile(example_bucket, name, ext).size > 0

    def test_exists(self, existing_example_file):
        name, ext = existing_example_file
        assert ObjectFile(example_bucket, name, ext).exists()

    def test_not_exists(self):
        assert not ObjectFile(example_bucket, 'does_not_exist').exists()

    def test_open(self, existing_example_file):
        name, ext = existing_example_file

        assert ObjectFile(example_bucket, name, ext).exists()
        with ObjectFile(example_bucket, name, ext=ext).open() as f:
            json.load(f)

    def test_delete(self, existing_example_file):
        name, ext = existing_example_file
        ObjectFile(example_bucket, name, ext).delete()
        assert not ObjectFile(example_bucket, name, ext).exists()

    def test_delete_all(self, existing_example_file):
        name, ext = existing_example_file
        Objects.delete_all(example_bucket)
        assert not ObjectFile(example_bucket, name, ext).exists()


class TestBaggedDataContainer:

    @pytest.fixture(scope='function')
    def example_directory(self, clear_files):
        directory = os.path.join(config.fs.tmp, 'test_container')
        os.makedirs(directory, exist_ok=True)

        with ZipFile(example_file) as zip_file:
            zip_file.extractall(directory)

        yield directory

    @pytest.fixture(scope='function')
    def example_container(self, example_directory):
        yield BaggedDataContainer.create(example_directory)

    def assert_container(self, container):
        assert container.manifest is not None
        assert len(container.manifest) == 5
        assert container.hash is not None
        assert container.metadata is not None

    def test_make(self, example_container):
        self.assert_container(example_container)

    def test_metadata(self, example_directory, example_container):
        example_container.metadata['test'] = dict(k1='v1', k2=True, k3=0)
        example_container.save_metadata()

        example_container = BaggedDataContainer(example_directory)
        self.assert_container(example_container)
        assert example_container.metadata['test']['k1'] == 'v1'
        assert example_container.metadata['test']['k2']
        assert example_container.metadata['test']['k3'] == 0

    def test_file(self, example_container):
        file = example_container.get_file('examples_template/template.json')
        assert file is not None
        with file.open('r') as f:
            assert json.load(f)


class TestZippedDataContainer(TestBaggedDataContainer):
    @pytest.fixture(scope='function')
    def example_container(self, example_directory):
        BaggedDataContainer.create(example_directory)
        return ZippedDataContainer.create(example_directory)

    def test_metadata(self, example_directory, example_container):
        pass


@pytest.fixture(scope='function', params=[False, True])
def archive_config(monkeypatch, request):
    new_config = config.FilesConfig(
        config.files.uploads_bucket,
        config.files.raw_bucket,
        config.files.archive_bucket,
        request.param)
    monkeypatch.setattr(config, 'files', new_config)
    yield


@pytest.fixture(scope='function')
def archive(clear_files, archive_config):
    archive = ArchiveFile('__test_upload_hash/__test_calc_hash')
    with archive.write_archive_json() as out:
        json.dump(example_data, out)
    yield archive


class TestArchiveFile:

    def test_archive(self, archive: ArchiveFile, no_warn):
        assert archive.exists()

        with archive.read_archive_json() as file:
            result = json.load(file)

        assert 'test_key' in result
        assert result['test_key'] == 'test_value'

    def test_delete_archive(self, archive: ArchiveFile, no_warn):
        archive.delete()
        assert not archive.exists()

    def test_delete_archives(self, archive: ArchiveFile, no_warn):
        ArchiveFile.delete_archives(archive.object_id.split('/')[0])
        assert not archive.exists()


class TestUploadFile:

    @pytest.fixture()
    def upload_same_file(self, clear_files):
        upload = UploadFile('__test_upload_id2')
        shutil.copyfile(example_file, upload.os_path)
        yield upload

    @pytest.fixture()
    def upload(self, clear_files):
        upload = UploadFile('__test_upload_id')
        shutil.copyfile(example_file, upload.os_path)
        yield upload

    def assert_upload(self, upload: UploadFile):
        assert upload.exists()

        assert len(upload.filelist) == 5
        has_json = False
        for filename in upload.filelist:
            the_file = upload.get_file(filename)
            assert the_file.exists()
            assert the_file.size >= 0
            if the_file.path.endswith('.json'):
                has_json = True
                assert the_file.size > 0
                with the_file.open() as f:
                    f.read()
                break
        assert has_json

    def test_upload_extracted(self, upload: UploadFile):
        with upload:
            self.assert_upload(upload)

    def test_persist(self, upload: UploadFile):
        with upload:
            zipped_container = upload.persist()

        assert zipped_container.exists()
        assert zipped_container.os_path.endswith('%s.zip' % upload.upload_hash())

    def test_delete_upload(self, upload: UploadFile):
        upload.delete()
        assert not upload.exists()

    def test_hash(self, upload: UploadFile, upload_same_file: UploadFile, no_warn):
        with upload:
            hash = upload.upload_hash()
            assert hash is not None
            assert isinstance(hash, str)

        with upload_same_file:
            assert hash == upload_same_file.upload_hash()

    def test_siblings(self, upload: UploadFile, no_warn):
        with upload:
            siblings = list(upload.get_siblings('examples_template/template.json'))
            assert len(siblings) == 4
            assert all(sibling.endswith('.aux') for sibling in siblings)


class TestLocalUploadFile(TestUploadFile):
    @pytest.fixture()
    def upload_same_file(self, clear_files):
        upload = UploadFile('__test_upload_id2', local_path=example_file)
        yield upload

    @pytest.fixture()
    def upload(self, clear_files):
        upload = UploadFile('__test_upload_id', local_path=example_file)
        yield upload

    def test_delete_upload(self, upload: UploadFile):
        upload.delete()
        assert upload.exists()


@pytest.fixture(scope='function')
def archive_log(clear_files, archive_config):
    archive_log = ArchiveLogFile('__test_upload_hash/__test_calc_hash')
    with archive_log.open('wt') as f:
        f.write('This is a test')

    yield archive_log


class TestArchiveLogFile:

    def test_archive_log_file(self, archive_log):
        assert archive_log.exists()
        with archive_log.open('rt') as f:
            assert 'test' in f.read()
