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

import os
import os.path
import shutil
import pytest

from nomad import config
from nomad.uploads import DirectoryObject, FileObject
from nomad.uploads import Metadata, MetadataTimeout
from nomad.uploads import StagingUploadFiles, PublicUploadFiles


class TestObjects:

    @pytest.fixture(scope='function')
    def test_bucket(self):
        yield 'test_bucket'

        bucket = os.path.join(config.fs.objects, 'test_bucket')
        if os.path.exists(bucket):
            shutil.rmtree(os.path.join(config.fs.objects, 'test_bucket'))

    def test_file_create(self, test_bucket):
        file = FileObject(test_bucket, 'sub/test_id', create=True)
        assert os.path.isdir(os.path.dirname(file.os_path))
        assert not file.exists()

    def test_file_dir_existing(self, test_bucket):
        file = FileObject(test_bucket, 'sub/test_id', create=False)
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


class TestMetadata:
    @pytest.fixture(scope='function')
    def test_dir(self):
        path = os.path.join(config.fs.tmp, 'test_dir')
        os.makedirs(path)
        yield path
        shutil.rmtree(path)

    @pytest.fixture(scope='function')
    def md(self, test_dir):
        with Metadata(test_dir) as md:
            yield md

    def test_open_empty(self, test_dir):
        with Metadata(test_dir):
            pass

    def test_modify(self, test_dir):
        with Metadata(test_dir) as md:
            md.data['key'] = 'value'

        with Metadata(test_dir) as md:
            assert 'key' in md.data
            assert md.data['key'] == 'value'
            assert len(md.data) == 1

    def test_lock(self, test_dir):
        timeout = False
        with Metadata(test_dir):
            try:
                with Metadata(test_dir, lock_timeout=0.1):
                    pass
            except MetadataTimeout:
                timeout = True
        assert timeout

    def test_insert(self, md: Metadata):
        md.insert(dict(hash='0', data='test'))
        assert len(md.data) == 1
        assert '0' in md.data
        assert md.data['0']['data'] == 'test'

    def test_insert_fail(self, md: Metadata):
        failed = False
        md.insert(dict(hash='0', data='test'))
        try:
            md.insert(dict(hash='0', data='test'))
        except Exception:
            failed = True

        assert failed
        assert len(md.data) == 1

    def test_update(self, md: Metadata):
        md.insert(dict(hash='0', data='test'))
        md.update(dict(hash='0', data='updated'))
        assert len(md.data) == 1
        assert md.data['0']['data'] == 'updated'

    def test_update_fail(self, md: Metadata):
        failed = False
        try:
            md.update(dict(hash='0', data='updated'))
        except KeyError:
            failed = True
        assert failed
        assert len(md.data) == 0

    def test_get(self, md: Metadata):
        md.insert(dict(hash='0', data='test'))
        assert md.get('0') is not None
        assert 'data' in md.get('0')

    def test_get_fail(self, md: Metadata):
        failed = False
        try:
            md.get('0')
        except KeyError:
            failed = True
        assert failed


class TestStagingUploadFiles:

    @pytest.fixture(scope='function')
    def test_upload(self):
        yield 'test_upload'
        for bucket in [config.files.staging_bucket, config.files.public_bucket]:
            directory = DirectoryObject(bucket, 'test_upload')
            if directory.exists():
                directory.delete()

    def test_create(self, test_upload):
        StagingUploadFiles(test_upload, create=True)

    def test_metadata_lock(self, test_upload):
        failed = False
        with StagingUploadFiles(test_upload, create=True):
            try:
                with StagingUploadFiles(test_upload):
                    pass
            except MetadataTimeout:
                failed = True
        assert failed
