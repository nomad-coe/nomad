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

from nomad.files import Objects, ObjectFile, ArchiveFile, UploadFile, ArchiveLogFile
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
            shutil.rmtree(config.fs.tmp)
        except FileNotFoundError:
            pass


class TestObjects:
    @pytest.fixture()
    def existing_example_file(self, clear_files):
        out = ObjectFile(example_bucket, 'example_file', ext='json').open(mode='wt')
        json.dump(example_data, out)
        out.close()

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
        file = ObjectFile(example_bucket, name, ext=ext).open()
        json.load(file)
        file.close()

    def test_delete(self, existing_example_file):
        name, ext = existing_example_file
        ObjectFile(example_bucket, name, ext).delete()
        assert not ObjectFile(example_bucket, name, ext).exists()

    def test_delete_all(self, existing_example_file):
        name, ext = existing_example_file
        Objects.delete_all(example_bucket)
        assert not ObjectFile(example_bucket, name, ext).exists()


@pytest.fixture(scope='function', params=[False, True])
def archive_config(monkeypatch, request):
    new_config = config.FilesConfig(
        config.files.uploads_bucket,
        config.files.repository_bucket,
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

    def test_upload(self, upload: UploadFile):
        assert upload.exists()

        with upload:
            assert len(upload.filelist) == 5
            # now just try to open the first file (not directory), without error
            for filename in upload.filelist:
                the_file = upload.get_file(filename)
                if the_file.os_path.endswith('.xml'):
                    the_file.open()
                    the_file.close()
                    break

    def test_delete_upload(self, upload: UploadFile):
        upload.delete()
        assert not upload.exists()

    def test_hash(self, upload: UploadFile, upload_same_file: UploadFile, no_warn):
        with upload:
            hash = upload.hash()
            assert hash is not None
            assert isinstance(hash, str)

        with upload_same_file:
            assert hash == upload_same_file.hash()


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
    f = archive_log.open('wt')
    f.write('This is a test')
    f.close()

    yield archive_log


class TestArchiveLogFile:

    def test_archive_log_file(self, archive_log):
        assert archive_log.exists()
        assert 'test' in archive_log.open('rt').read()
