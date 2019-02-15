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

from typing import Generator, Any, Dict
import os
import os.path
import shutil
import pytest
import json

from nomad import config
from nomad.files import DirectoryObject, PathObject
from nomad.files import Metadata, PublicMetadata, StagingMetadata
from nomad.files import StagingUploadFiles, PublicUploadFiles, UploadFiles, Restricted, \
    ArchiveBasedStagingUploadFiles


# example_file uses an artificial parser for faster test execution, can also be
# changed to examples_vasp.zip for using vasp parser
example_file = 'tests/data/proc/examples_template.zip'
example_file_contents = [
    'examples_template/template.json',
    'examples_template/1.aux',
    'examples_template/2.aux',
    'examples_template/3.aux',
    'examples_template/4.aux']
example_file_mainfile = 'examples_template/template.json'
empty_file = 'tests/data/proc/empty.zip'

example_bucket = 'test_bucket'
example_data = dict(test_key='test_value')


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

    @pytest.mark.parametrize('dirpath', ['test', os.path.join('sub', 'test')])
    @pytest.mark.parametrize('create', [True, False])
    @pytest.mark.parametrize('prefix', [True, False])
    def test_directory(self, test_bucket: str, dirpath: str, create: bool, prefix: bool) -> None:
        directory = DirectoryObject(test_bucket, dirpath, create=create, prefix=prefix)
        assert directory.exists() == create
        assert os.path.isdir(directory.os_path) == create
        assert directory.os_path.endswith(os.path.join('tes' if prefix else '', 'test'))

    @pytest.mark.parametrize('dirpath', ['test', os.path.join('sub', 'test')])
    @pytest.mark.parametrize('create', [True, False])
    @pytest.mark.parametrize('join_create', [True, False])
    @pytest.mark.parametrize('prefix', [True, False])
    def test_directory_join(self, test_bucket: str, dirpath: str, create: bool, prefix: bool, join_create: bool) -> None:
        directory = DirectoryObject(test_bucket, 'parent', create=create, prefix=prefix)
        directory = directory.join_dir(dirpath, create=join_create)

        assert directory.exists() == join_create
        assert os.path.isdir(directory.os_path) == join_create
        assert dirpath.endswith(os.path.join('', 'test'))

    @pytest.mark.parametrize('filepath', ['test', 'sub/test'])
    @pytest.mark.parametrize('create', [True, False])
    def test_directory_join_file_dir_create(self, test_bucket: str, filepath: str, create: bool):
        directory = DirectoryObject(test_bucket, 'parent', create=create)
        file = directory.join_file(filepath)
        assert os.path.exists(directory.os_path) == create
        assert os.path.exists(os.path.dirname(file.os_path)) == create


example_calc: Dict[str, Any] = {
    'calc_id': '0',
    'mainfile': 'examples_template/template.json',
    'data': 'value'
}
example_calc_id = example_calc['calc_id']


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

    def test_get(self, md: Metadata):
        assert_example_calc(md.get(example_calc_id))

    def test_get_fail(self, md: Metadata):
        failed = False
        try:
            md.get('unknown')
        except KeyError:
            failed = True
        assert failed


class TestStagingMetadata(MetadataContract):
    @pytest.fixture(scope='function')
    def md(self, test_dir):
        md = StagingMetadata(DirectoryObject(None, None, os_path=test_dir))
        md.insert(example_calc)
        return md

    def test_remove(self, md: StagingMetadata):
        md.remove(example_calc)
        failed = False
        try:
            assert md.get(example_calc['calc_id'])
        except KeyError:
            failed = True
        assert failed

    def test_insert(self, md: StagingMetadata):
        md.remove(example_calc)
        md.insert(example_calc)
        assert len(md) == 1
        assert_example_calc(md.get(example_calc_id))

    def test_insert_fail(self, md: StagingMetadata):
        failed = False
        try:
            md.insert(example_calc)
        except Exception:
            failed = True

        assert failed
        assert len(md) == 1

    def test_update(self, md: StagingMetadata):
        md.update(example_calc_id, dict(data='updated'))
        assert len(md) == 1
        assert md.get(example_calc_id)['data'] == 'updated'

    def test_update_fail(self, md: StagingMetadata):
        failed = False
        try:
            md.update('unknown', dict(data='updated'))
        except KeyError:
            failed = True
        assert failed
        assert len(md) == 1


class TestPublicMetadata(MetadataContract):

    @pytest.fixture(scope='function')
    def md(self, test_dir):
        md = PublicMetadata(test_dir)
        md._create([example_calc])
        return md


class UploadFilesFixtures:

    @pytest.fixture(scope='function')
    def test_upload_id(self) -> Generator[str, None, None]:
        for bucket in [config.files.staging_bucket, config.files.public_bucket]:
            directory = DirectoryObject(bucket, 'test_upload', prefix=True)
            if directory.exists():
                directory.delete()
        yield 'test_upload'
        for bucket in [config.files.staging_bucket, config.files.public_bucket]:
            directory = DirectoryObject(bucket, 'test_upload', prefix=True)
            if directory.exists():
                directory.delete()


class UploadFilesContract(UploadFilesFixtures):

    @pytest.fixture(scope='function', params=['r'])
    def test_upload(self, request, test_upload_id) -> UploadFiles:
        raise NotImplementedError()

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id) -> Generator[UploadFiles, None, None]:
        raise NotImplementedError()

    def test_create(self, empty_test_upload):
        assert UploadFiles.get(empty_test_upload.upload_id).__class__ == empty_test_upload.__class__

    def test_rawfile(self, test_upload):
        try:
            with test_upload.raw_file(example_file_mainfile) as f:
                assert len(f.read()) > 0
            if not test_upload._is_authorized():
                assert not test_upload.metadata.get(example_calc_id).get('restricted', False)
        except Restricted:
            assert not test_upload._is_authorized()
            assert test_upload.metadata.get(example_calc_id).get('restricted', False)

    @pytest.mark.parametrize('prefix', [None, 'examples'])
    def test_raw_file_manifest(self, test_upload: StagingUploadFiles, prefix: str):
        raw_files = list(test_upload.raw_file_manifest(path_prefix=prefix))
        assert sorted(file for file in raw_files if file.startswith('examples')) == sorted(example_file_contents)

    @pytest.mark.parametrize('test_logs', [True, False])
    def test_archive(self, test_upload, test_logs: bool):
        try:
            if test_logs:
                with test_upload.archive_log_file(example_calc_id, 'rt') as f:
                    assert f.read() == 'archive'
            else:
                f = test_upload.archive_file(example_calc_id, 'rt')
                assert json.load(f) == 'archive'

            if not test_upload._is_authorized():
                assert not test_upload.metadata.get(example_calc_id).get('restricted', False)
        except Restricted:
            assert not test_upload._is_authorized()
            assert test_upload.metadata.get(example_calc_id).get('restricted', False)

    def test_metadata(self, test_upload):
        assert_example_calc(test_upload.metadata.get(example_calc_id))


def create_staging_upload(upload_id: str, calc_specs: str) -> StagingUploadFiles:
    """
    Create an upload according to given spec. Additional arguments are given to
    the StagingUploadFiles contstructor.

    Arguments:
        upload_id: The id that should be given to this test upload.
        calc_specs: A string that determines the properties of the given upload.
            With letters determining example calcs being public `p` or restricted `p`.
            The calcs will be copies of calcs in `example_file`.
            First calc is at top level, following calcs will be put under 1/, 2/, etc.
    """
    upload = StagingUploadFiles(upload_id, create=True, is_authorized=lambda: True)

    prefix = 0
    for calc_spec in calc_specs:
        upload.add_rawfiles(example_file, prefix=None if prefix == 0 else str(prefix))
        calc_id = str(int(example_calc_id) + prefix)
        with upload.archive_file(calc_id, 'wt') as f:
            f.write('"archive"')
        with upload.archive_log_file(calc_id, 'wt') as f:
            f.write('archive')
        calc = dict(**example_calc)
        calc['calc_id'] = calc_id
        if prefix > 0:
            calc['mainfile'] = os.path.join(str(prefix), calc['mainfile'])
        if calc_spec == 'r':
            calc['restricted'] = True
        elif calc_spec == 'p':
            calc['restricted'] = False
        upload.metadata.insert(calc)
        prefix += 1

    if calc_specs.startswith('P'):
        public_only = True
        calc_specs = calc_specs[1:]
    else:
        public_only = False
    upload._is_authorized = lambda: not public_only

    assert len(upload.metadata) == len(calc_specs)
    return upload


class TestStagingUploadFiles(UploadFilesContract):

    @pytest.fixture(scope='function', params=['r', 'rr', 'pr', 'rp', 'p', 'pp'])
    def test_upload(self, request, test_upload_id: str) -> StagingUploadFiles:
        return create_staging_upload(test_upload_id, calc_specs=request.param)

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id) -> Generator[UploadFiles, None, None]:
        yield StagingUploadFiles(test_upload_id, create=True, is_authorized=lambda: True)

    @pytest.mark.parametrize('prefix', [None, 'prefix'])
    def test_add_rawfiles_zip(self, test_upload_id, prefix):
        test_upload = StagingUploadFiles(test_upload_id, create=True, is_authorized=lambda: True)
        test_upload.add_rawfiles(example_file, prefix=prefix)
        for filepath in example_file_contents:
            filepath = os.path.join(prefix, filepath) if prefix else filepath
            with test_upload.raw_file(filepath) as f:
                content = f.read()
                if filepath == example_file_mainfile:
                    assert len(content) > 0

    def test_write_archive(self, test_upload):
        assert json.load(test_upload.archive_file(example_calc_id, 'rt')) == 'archive'

    def test_calc_id(self, test_upload):
        assert test_upload.calc_id(example_file_mainfile) is not None

    def test_pack(self, test_upload):
        test_upload.pack()

    @pytest.mark.parametrize('with_mainfile', [True, False])
    def test_calc_files(self, test_upload: StagingUploadFiles, with_mainfile):
        for calc in test_upload.metadata:
            mainfile = calc['mainfile']
            calc_files = test_upload.calc_files(mainfile, with_mainfile=with_mainfile)
            assert len(list(calc_files)) == len(example_file_contents) - 0 if with_mainfile else 1
            if with_mainfile:
                for one, two in zip(calc_files, [mainfile] + sorted(example_file_contents[1:])):
                    assert one.endswith(two)
                    assert one.startswith(mainfile[:3])

    def test_delete(self, test_upload: StagingUploadFiles):
        test_upload.delete()
        assert not test_upload.exists()

    def test_update_metadata(self, test_upload):
        test_upload.metadata.update(example_calc_id, dict(data='updated'))
        test_upload.metadata.get(example_calc_id)['data'] == 'updated'


class TestArchiveBasedStagingUploadFiles(UploadFilesFixtures):
    def test_create(self, test_upload_id):
        test_upload = ArchiveBasedStagingUploadFiles(test_upload_id, create=True)
        shutil.copy(example_file, test_upload.upload_file_os_path)
        test_upload.extract()
        assert sorted(list(test_upload.raw_file_manifest())) == sorted(example_file_contents)
        assert os.path.exists(test_upload.upload_file_os_path)

    def test_local_path(self, test_upload_id):
        test_upload = ArchiveBasedStagingUploadFiles(test_upload_id, create=True, local_path=example_file)
        test_upload.extract()
        assert sorted(list(test_upload.raw_file_manifest())) == sorted(example_file_contents)
        assert os.path.exists(test_upload.upload_file_os_path)

    def test_invalid(self, test_upload_id):
        assert ArchiveBasedStagingUploadFiles(test_upload_id, create=True, local_path=example_file).is_valid
        assert not ArchiveBasedStagingUploadFiles(test_upload_id, create=True).is_valid


def create_public_upload(upload_id: str, calc_specs: str, **kwargs):
    staging_upload = create_staging_upload(upload_id, calc_specs)
    staging_upload.pack()
    staging_upload.delete()
    return PublicUploadFiles(upload_id, **kwargs)


class TestPublicUploadFiles(UploadFilesContract):

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id: str) -> Generator[UploadFiles, None, None]:
        yield create_public_upload(test_upload_id, calc_specs='', is_authorized=lambda: True)

    @pytest.fixture(scope='function', params=['r', 'rr', 'pr', 'rp', 'p', 'pp', 'Ppr', 'Prp'])
    def test_upload(self, request, test_upload_id: str) -> PublicUploadFiles:
        calc_specs = request.param
        if calc_specs.startswith('P'):
            public_only = True
            calc_specs = calc_specs[1:]
        else:
            public_only = False

        staging_upload = create_staging_upload(test_upload_id, calc_specs=calc_specs)
        staging_upload.pack()
        return PublicUploadFiles(test_upload_id, is_authorized=lambda: not public_only)
