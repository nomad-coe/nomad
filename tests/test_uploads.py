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

from typing import Generator
import os
import os.path
import shutil
import pytest

from nomad import config
from nomad.uploads import DirectoryObject, PathObject
from nomad.uploads import Metadata, MetadataTimeout, PublicMetadata, StagingMetadata
from nomad.uploads import StagingUploadFiles, PublicUploadFiles, UploadFiles

from tests.test_files import example_file, example_file_contents, example_file_mainfile


class TestObjects:

    @pytest.fixture(scope='function')
    def test_bucket(self):
        yield 'test_bucket'

        bucket = os.path.join(config.fs.objects, 'test_bucket')
        if os.path.exists(bucket):
            shutil.rmtree(os.path.join(config.fs.objects, 'test_bucket'))

    def test_file_dir_existing(self, test_bucket):
        file = PathObject(test_bucket, 'sub/test_id')
        assert not os.path.exists(os.path.dirname(file.os_path))

    def test_directory_create(self, test_bucket):
        directory = DirectoryObject(test_bucket, 'sub/test_id', create=True)
        assert directory.exists()
        assert os.path.isdir(directory.os_path)

    def test_directory_existing(self, test_bucket):
        directory = DirectoryObject(test_bucket, 'sub/test_id', create=False)
        assert not directory.exists()
        assert not os.path.exists(directory.os_path)
        assert not os.path.exists(os.path.dirname(directory.os_path))

    def test_directory_join_dir_create(self, test_bucket):
        directory = DirectoryObject(test_bucket, 'sub/parent', create=True)
        directory = directory.join_dir('test_id')
        assert directory.exists()
        assert os.path.isdir(directory.os_path)

    def test_directory_join_dir_existing(self, test_bucket):
        directory = DirectoryObject(test_bucket, 'sub/parent', create=False)
        directory = directory.join_dir('test_id')
        assert not directory.exists()
        assert not os.path.exists(directory.os_path)
        assert not os.path.exists(os.path.dirname(directory.os_path))

    def test_directory_join_dir_join_create(self, test_bucket):
        directory = DirectoryObject(test_bucket, 'sub/parent', create=False)
        directory = directory.join_dir('test_id', create=True)
        assert directory.exists()
        assert os.path.isdir(directory.os_path)

    def test_directory_join_dir_join_exist(self, test_bucket):
        directory = DirectoryObject(test_bucket, 'sub/parent', create=True)
        directory = directory.join_dir('test_id', create=False)
        assert not directory.exists()
        assert not os.path.exists(directory.os_path)
        assert os.path.exists(os.path.dirname(directory.os_path))

    def test_directory_join_file_dir_create(self, test_bucket):
        directory = DirectoryObject(test_bucket, 'sub/parent', create=True)
        file = directory.join_file('test_id')
        assert os.path.exists(directory.os_path)
        assert os.path.exists(os.path.dirname(file.os_path))

    def test_directory_join_file_dir_existing(self, test_bucket):
        directory = DirectoryObject(test_bucket, 'sub/parent', create=False)
        file = directory.join_file('test_id')
        assert not os.path.exists(directory.os_path)
        assert not os.path.exists(os.path.dirname(file.os_path))

    def test_directory_prefix(self, test_bucket):
        directory = DirectoryObject(test_bucket, 'sub/parent', create=True, prefix=True)
        assert os.path.join('sub/par/parent') in directory.os_path
        assert os.path.isdir(directory.os_path)


example_calc = {
    'hash': '0',
    'mainfile': 'examples_template/template.json',
    'data': 'value'
}
example_calc_hash = example_calc['hash']


def assert_example_calc(calc):
    assert calc is not None
    assert calc['data'] == example_calc['data']


class MetadataContract:
    @pytest.fixture(scope='function')
    def test_dir(self):
        path = os.path.join(config.fs.tmp, 'test_dir')
        os.makedirs(path)
        yield path
        shutil.rmtree(path)

    @pytest.fixture(scope='function')
    def md(self, test_dir):
        raise NotImplementedError()

    def test_open_empty(self, md):
        pass

    def test_insert(self, md: Metadata):
        md.insert(example_calc)
        assert len(md) == 1
        assert_example_calc(md.get(example_calc_hash))

    def test_insert_fail(self, md: Metadata):
        failed = False
        md.insert(example_calc)
        try:
            md.insert(example_calc)
        except Exception:
            failed = True

        assert failed
        assert len(md) == 1

    def test_update(self, md: Metadata):
        md.insert(example_calc)
        md.update(example_calc_hash, dict(data='updated'))
        assert len(md) == 1
        assert md.get(example_calc_hash)['data'] == 'updated'

    def test_update_fail(self, md: Metadata):
        failed = False
        try:
            md.update(example_calc_hash, dict(data='updated'))
        except KeyError:
            failed = True
        assert failed
        assert len(md) == 0

    def test_get(self, md: Metadata):
        md.insert(example_calc)
        assert_example_calc(md.get(example_calc_hash))

    def test_get_fail(self, md: Metadata):
        failed = False
        try:
            md.get(example_calc_hash)
        except KeyError:
            failed = True
        assert failed


class TestStagingMetadata(MetadataContract):
    @pytest.fixture(scope='function')
    def md(self, test_dir):
        with StagingMetadata(DirectoryObject(None, None, os_path=test_dir)) as md:
            yield md


class TestPublicMetadata(MetadataContract):

    @pytest.fixture(scope='function')
    def md(self, test_dir):
        with PublicMetadata(test_dir) as md:
            yield md

    def test_lock(self, test_dir):
        timeout = False
        with PublicMetadata(test_dir):
            try:
                with PublicMetadata(test_dir, lock_timeout=0.1):
                    pass
            except MetadataTimeout:
                timeout = True
        assert timeout


class UploadFilesContract:
    @pytest.fixture(scope='function')
    def test_upload_id(self) -> Generator[str, None, None]:
        yield 'test_upload'
        for bucket in [config.files.staging_bucket, config.files.public_bucket]:
            directory = DirectoryObject(bucket, 'test_upload', prefix=True)
            if directory.exists():
                directory.delete()

    @pytest.fixture(scope='function')
    def test_upload(self, test_upload_id) -> Generator[UploadFiles, None, None]:
        raise NotImplementedError()

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id) -> Generator[UploadFiles, None, None]:
        raise NotImplementedError()

    def test_create(self, empty_test_upload):
        pass

    def test_rawfile(self, test_upload):
        with test_upload.raw_file(example_file_mainfile) as f:
            assert len(f.read()) > 0

    def test_archive(self, test_upload):
        with test_upload.archive_file(example_calc_hash) as f:
            assert f.read() == b'archive'

    def test_metadata(self, test_upload):
        with test_upload.metadata as md:
            assert_example_calc(md.get(example_calc_hash))

    def test_update_metadata(self, test_upload):
        with test_upload.metadata as md:
            md.update(example_calc_hash, dict(data='updated'))

        with test_upload.metadata as md:
            assert md.get(example_calc_hash)['data'] == 'updated'


class TestStagingUploadFiles(UploadFilesContract):
    @staticmethod
    def create_test_upload(test_upload_id):
        upload = StagingUploadFiles(test_upload_id, create=True, archive_ext='txt')
        upload.add_rawfiles(example_file)
        with upload.archive_file(example_calc_hash, read=False) as f:
            f.write(b'archive')
        upload.metadata.insert(example_calc)
        return upload

    @pytest.fixture(scope='function')
    def test_upload(self, test_upload_id):
        yield TestStagingUploadFiles.create_test_upload(test_upload_id)

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id) -> Generator[UploadFiles, None, None]:
        yield StagingUploadFiles(test_upload_id, create=True)

    def test_add_rawfiles_zip(self, test_upload):
        for filepath in example_file_contents:
            with test_upload.raw_file(filepath) as f:
                content = f.read()
                if filepath == example_file_mainfile:
                    assert len(content) > 0

    def test_write_archive(self, test_upload):
        with test_upload.archive_file(example_calc_hash) as f:
            assert f.read() == b'archive'

    def test_calc_hash(self, test_upload):
        assert test_upload.calc_hash(example_file_mainfile) is not None

    def test_pack(self, test_upload):
        test_upload.pack()

    def test_all_rawfiles(self, test_upload: StagingUploadFiles):
        for filepath in test_upload.all_rawfiles:
            assert os.path.isfile(filepath)


class TestPublicUploadFiles(UploadFilesContract):

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id) -> Generator[UploadFiles, None, None]:
        yield PublicUploadFiles(test_upload_id, 'txt')

    @pytest.fixture(scope='function')
    def test_upload(self, test_upload_id):
        staging_upload = TestStagingUploadFiles.create_test_upload(test_upload_id)
        staging_upload.pack()
        yield PublicUploadFiles(test_upload_id, 'txt')
