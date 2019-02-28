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
import itertools
import zipfile
import re

from nomad import config
from nomad.files import DirectoryObject, PathObject
from nomad.files import Metadata, PublicMetadata, StagingMetadata
from nomad.files import StagingUploadFiles, PublicUploadFiles, UploadFiles, Restricted, \
    ArchiveBasedStagingUploadFiles


# example_file uses an artificial parser for faster test execution, can also be
# changed to examples_vasp.zip for using vasp parser
example_file = 'tests/data/proc/examples_template.zip'
example_directory = 'tests/data/proc/examples_template'
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

    def test_delete_prefix(self, test_bucket: str):
        dir_1 = DirectoryObject(test_bucket, 'test_directory1', create=True, prefix=True)
        dir_2 = DirectoryObject(test_bucket, 'test_directory2', create=True, prefix=True)
        dir_1.delete()
        dir_2.delete()

        prefix = os.path.dirname(dir_2.os_path)
        assert len(os.path.basename(prefix)) == 3
        assert not os.path.exists(prefix)


example_calc: Dict[str, Any] = {
    'calc_id': '0',
    'mainfile': 'examples_template/template.json',
    'data': 'value'
}
example_calc_id = example_calc['calc_id']


def generate_example_calc(calc_id: int, with_mainfile_prefix: bool, subdirectory: str = None, **kwargs):
    example_calc = dict(calc_id=str(calc_id), data='value')
    if with_mainfile_prefix:
        mainfile = '%d.template.json' % calc_id
    else:
        mainfile = 'template.json'

    if subdirectory is not None:
        mainfile = os.path.join(subdirectory, mainfile)

    example_calc['mainfile'] = mainfile
    example_calc.update(**kwargs)

    example_file = os.path.join(config.fs.tmp, 'example.zip')
    with zipfile.ZipFile(example_file, 'w', zipfile.ZIP_DEFLATED) as zf:
        for filepath in example_file_contents:
            filename = os.path.basename(filepath)
            arcname = filename
            if arcname == 'template.json' and with_mainfile_prefix:
                arcname = '%d.template.json' % calc_id

            if subdirectory is not None:
                arcname = os.path.join(subdirectory, arcname)
            zf.write(os.path.join(example_directory, filename), arcname)

    return example_calc, example_file


def assert_example_files(names, with_mainfile: bool = True):
    # TODO its complicated
    # To compare the files with the example_file_contents list we have to assume
    # - different subdirectories
    # - mainfile prefixes
    # - mainfiles among aux files
    is_multi = any(re.search(r'[0-9].t', name) for name in names)

    def normalized_file(name):
        name = re.sub(r'[0-9].t', 't', name)
        name = re.sub(r'^[0-9]\/', '', name)
        return name

    source = sorted(set(normalized_file(name) for name in names if not name.endswith('template.json') or with_mainfile or not is_multi))
    target = sorted(name for name in example_file_contents if not name.endswith('template.json') or with_mainfile)
    assert source == target


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
        assert len(test_upload.metadata) > 0
        for calc in test_upload.metadata:
            try:
                with test_upload.raw_file(calc['mainfile']) as f:
                    assert len(f.read()) > 0
                if not test_upload._is_authorized():
                    assert not test_upload.metadata.get(calc['calc_id']).get('with_embargo', False)
            except Restricted:
                assert not test_upload._is_authorized()
                assert test_upload.metadata.get(calc['calc_id']).get('with_embargo', False)

    @pytest.mark.parametrize('prefix', [None, 'examples'])
    def test_raw_file_manifest(self, test_upload: StagingUploadFiles, prefix: str):
        raw_files = list(test_upload.raw_file_manifest(path_prefix=prefix))
        assert_example_files(raw_files)

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
                assert not test_upload.metadata.get(example_calc_id).get('with_embargo', False)
        except Restricted:
            assert not test_upload._is_authorized()
            assert test_upload.metadata.get(example_calc_id).get('with_embargo', False)

    def test_metadata(self, test_upload):
        assert_example_calc(test_upload.metadata.get(example_calc_id))


def create_staging_upload(upload_id: str, calc_specs: str) -> StagingUploadFiles:
    """
    Create an upload according to given spec. Additional arguments are given to
    the StagingUploadFiles contstructor.

    Arguments:
        upload_id: The id that should be given to this test upload.
        calc_specs: A string that determines the properties of the given upload.
            With letters determining example calcs being public `p` or restricted `r`.
            The calcs will be copies of calcs in `example_file`.
            First calc is at top level, following calcs will be put under 1/, 2/, etc.
            All calcs with capital `P`/`R` will be put in the same directory under multi/.
    """
    upload = StagingUploadFiles(upload_id, create=True, is_authorized=lambda: True)

    prefix = 0
    for calc_spec in calc_specs:
        is_multi = calc_spec in ['R', 'P']
        calc_spec = calc_spec.lower()

        if is_multi or prefix == 0:
            directory = 'examples_template'
        else:
            directory = os.path.join(str(prefix), 'examples_template')

        calc, upload_file = generate_example_calc(
            prefix, with_mainfile_prefix=is_multi, subdirectory=directory, with_embargo=calc_spec == 'r')
        calc_id = calc['calc_id']

        upload.add_rawfiles(upload_file)

        with upload.archive_file(calc_id, 'wt') as f:
            f.write('"archive"')
        with upload.archive_log_file(calc_id, 'wt') as f:
            f.write('archive')

        upload.metadata.insert(calc)
        prefix += 1

    assert len(upload.metadata) == len(calc_specs)
    return upload


class TestStagingUploadFiles(UploadFilesContract):

    @pytest.fixture(scope='function', params=['r', 'rr', 'pr', 'rp', 'p', 'pp', 'RP', 'RR', 'PP'])
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
            assert_example_files(calc_files, with_mainfile=with_mainfile)

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

    @pytest.fixture(scope='function', params=itertools.product(
        ['r', 'rr', 'pr', 'rp', 'p', 'pp', 'RP', 'RR', 'PP'], [True, False]))
    def test_upload(self, request, test_upload_id: str) -> PublicUploadFiles:
        calc_specs, protected = request.param
        staging_upload = create_staging_upload(test_upload_id, calc_specs=calc_specs)
        staging_upload.pack()
        return PublicUploadFiles(test_upload_id, is_authorized=lambda: not protected)
