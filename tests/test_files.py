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

from typing import Generator, Any, Dict, Tuple
import os
import os.path
import shutil
import pytest
import json
import itertools
import zipfile
import re

from nomad import config
from nomad.datamodel import UploadWithMetadata, CalcWithMetadata
from nomad.files import DirectoryObject, PathObject
from nomad.files import StagingUploadFiles, PublicUploadFiles, UploadFiles, Restricted, \
    ArchiveBasedStagingUploadFiles

from tests.utils import assert_exception


CalcWithFiles = Tuple[CalcWithMetadata, str]
UploadWithFiles = Tuple[UploadWithMetadata, UploadFiles]
StagingUploadWithFiles = Tuple[UploadWithMetadata, StagingUploadFiles]
PublicUploadWithFiles = Tuple[UploadWithMetadata, PublicUploadFiles]

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


@pytest.fixture(scope='function', autouse=True)
def raw_files_on_all_tests(raw_files):
    """ Autouse fixture to apply raw_files to all tests. """
    pass


class TestObjects:

    @pytest.fixture(scope='function')
    def test_bucket(self):
        yield config.fs.staging

        bucket = os.path.join(config.fs.staging)
        if os.path.exists(bucket):
            shutil.rmtree(os.path.join(config.fs.staging))

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
        assert directory.os_path.endswith(os.path.join('te' if prefix else '', 'test'))

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
        assert len(os.path.basename(prefix)) == 2
        assert not os.path.exists(prefix)


example_calc: Dict[str, Any] = {
    'calc_id': '0',
    'mainfile': 'examples_template/template.json',
    'data': 'value'
}
example_calc_id = example_calc['calc_id']


def generate_example_calc(
        calc_id: int, with_mainfile_prefix: bool, subdirectory: str = None,
        **kwargs) -> CalcWithFiles:
    """ Generate an example calc with :class:`CalcWithMetadata` and rawfile. """

    example_calc = CalcWithMetadata(calc_id=str(calc_id))

    if with_mainfile_prefix:
        mainfile = '%d.template.json' % calc_id
    else:
        mainfile = 'template.json'

    if subdirectory is not None:
        mainfile = os.path.join(subdirectory, mainfile)

    example_calc.mainfile = mainfile
    example_calc.update(**kwargs)

    example_file = os.path.join(config.fs.tmp, 'example.zip')
    example_calc.files = []
    with zipfile.ZipFile(example_file, 'w', zipfile.ZIP_DEFLATED) as zf:
        for filepath in example_file_contents:
            filename = os.path.basename(filepath)
            arcname = filename
            if arcname == 'template.json' and with_mainfile_prefix:
                arcname = '%d.template.json' % calc_id

            if subdirectory is not None:
                arcname = os.path.join(subdirectory, arcname)
            example_calc.files.append(arcname)
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


class UploadFilesFixtures:

    @pytest.fixture(scope='function')
    def test_upload_id(self) -> Generator[str, None, None]:
        for bucket in [config.fs.staging, config.fs.public]:
            directory = DirectoryObject(bucket, 'test_upload', prefix=True)
            if directory.exists():
                directory.delete()
        yield 'test_upload'
        for bucket in [config.fs.staging, config.fs.public]:
            directory = DirectoryObject(bucket, 'test_upload', prefix=True)
            if directory.exists():
                directory.delete()


class UploadFilesContract(UploadFilesFixtures):

    @pytest.fixture(scope='function', params=['r'])
    def test_upload(self, request, test_upload_id) -> UploadWithFiles:
        raise NotImplementedError()

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id) -> UploadFiles:
        raise NotImplementedError()

    def test_create(self, empty_test_upload):
        assert UploadFiles.get(empty_test_upload.upload_id).__class__ == empty_test_upload.__class__

    def test_rawfile(self, test_upload: UploadWithFiles):
        upload, upload_files = test_upload
        for calc in upload.calcs:
            try:
                for file_path in calc.files:
                    with upload_files.raw_file(file_path) as f:
                        assert len(f.read()) > 0
                    if not upload_files._is_authorized():
                        assert not calc.with_embargo
            except Restricted:
                assert not upload_files._is_authorized()
                assert calc.with_embargo

    @pytest.mark.parametrize('prefix', [None, 'examples'])
    def test_raw_file_manifest(self, test_upload: UploadWithFiles, prefix: str):
        _, upload_files = test_upload
        raw_files = list(upload_files.raw_file_manifest(path_prefix=prefix))
        assert_example_files(raw_files)

    @pytest.mark.parametrize('prefix', [None, 'examples_template'])
    def test_raw_file_list(self, test_upload: UploadWithFiles, prefix: str):
        _, upload_files = test_upload
        raw_files = list(upload_files.raw_file_list(directory=prefix))
        if prefix is None:
            assert len(raw_files) == 0
        elif upload_files._is_authorized() or len(raw_files) > 0:
            assert '1.aux' in list(path for path, _ in raw_files)
            for file, size in raw_files:
                if file.endswith('.aux'):
                    assert size == 8
                else:
                    assert size > 0
            assert_example_files([os.path.join(prefix, path) for path, _ in raw_files])

    @pytest.mark.parametrize('test_logs', [True, False])
    def test_archive(self, test_upload: UploadWithFiles, test_logs: bool):
        upload, upload_files = test_upload
        calcs = upload.calcs_dict
        try:
            if test_logs:
                with upload_files.archive_log_file(example_calc_id, 'rt') as f:
                    assert f.read() == 'archive'
            else:
                f = upload_files.archive_file(example_calc_id, 'rt')
                assert json.load(f) == 'archive'

            if not upload_files._is_authorized():
                assert not calcs.get(example_calc_id).with_embargo
        except Restricted:
            assert not upload_files._is_authorized()
            assert calcs.get(example_calc_id).with_embargo


def create_staging_upload(upload_id: str, calc_specs: str) -> StagingUploadWithFiles:
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
    upload_files = StagingUploadFiles(upload_id, create=True, is_authorized=lambda: True)
    upload = UploadWithMetadata(upload_id=upload_id)
    calcs = []

    prefix = 0
    for calc_spec in calc_specs:
        is_multi = calc_spec in ['R', 'P']
        calc_spec = calc_spec.lower()

        if is_multi or prefix == 0:
            directory = 'examples_template'
        else:
            directory = os.path.join(str(prefix), 'examples_template')

        calc, calc_file = generate_example_calc(
            prefix, with_mainfile_prefix=is_multi, subdirectory=directory,
            with_embargo=calc_spec == 'r')

        upload_files.add_rawfiles(calc_file)

        with upload_files.archive_file(calc.calc_id, 'wt') as f:
            f.write('"archive"')
        with upload_files.archive_log_file(calc.calc_id, 'wt') as f:
            f.write('archive')

        calcs.append(calc)
        prefix += 1

    assert len(calcs) == len(calc_specs)
    upload.calcs = calcs
    return upload, upload_files


class TestStagingUploadFiles(UploadFilesContract):

    @pytest.fixture(scope='function', params=['r', 'rr', 'pr', 'rp', 'p', 'pp', 'RP', 'RR', 'PP'])
    def test_upload(self, request, test_upload_id: str) -> StagingUploadWithFiles:
        return create_staging_upload(test_upload_id, calc_specs=request.param)

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id) -> UploadFiles:
        return StagingUploadFiles(test_upload_id, create=True, is_authorized=lambda: True)

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

    def test_write_archive(self, test_upload: StagingUploadWithFiles):
        _, upload_files = test_upload
        assert json.load(upload_files.archive_file(example_calc_id, 'rt')) == 'archive'

    def test_calc_id(self, test_upload: StagingUploadWithFiles):
        _, upload_files = test_upload
        assert upload_files.calc_id(example_file_mainfile) is not None

    def test_pack(self, test_upload: StagingUploadWithFiles):
        upload, upload_files = test_upload
        upload_files.pack(upload)

    @pytest.mark.parametrize('with_mainfile', [True, False])
    def test_calc_files(self, test_upload: StagingUploadWithFiles, with_mainfile):
        upload, upload_files = test_upload
        for calc in upload.calcs:
            mainfile = calc.mainfile
            calc_files = upload_files.calc_files(mainfile, with_mainfile=with_mainfile)
            assert_example_files(calc_files, with_mainfile=with_mainfile)

    def test_delete(self, test_upload: StagingUploadWithFiles):
        _, upload_files = test_upload
        upload_files.delete()
        assert not upload_files.exists()

    def test_create_extracted_copy(self, test_upload: StagingUploadWithFiles):
        upload, upload_files = test_upload
        upload_files.create_extracted_copy()
        for calc in upload.calcs:
            assert os.path.exists(os.path.join(
                config.fs.coe_extracted, upload_files.upload_id, calc.mainfile))


class TestArchiveBasedStagingUploadFiles(UploadFilesFixtures):
    def test_create(self, test_upload_id):
        test_upload = ArchiveBasedStagingUploadFiles(
            test_upload_id, create=True, upload_path=example_file)
        test_upload.extract()
        assert sorted(list(test_upload.raw_file_manifest())) == sorted(example_file_contents)
        assert os.path.exists(test_upload.upload_path)

    def test_invalid(self, test_upload_id):
        assert ArchiveBasedStagingUploadFiles(
            test_upload_id, create=True, upload_path=example_file).is_valid
        assert not ArchiveBasedStagingUploadFiles(
            test_upload_id, create=True, upload_path='does not exist').is_valid


def create_public_upload(
        upload_id: str, calc_specs: str, **kwargs) -> PublicUploadWithFiles:

    upload, upload_files = create_staging_upload(upload_id, calc_specs)
    upload_files.pack(upload)
    upload_files.delete()
    return upload, PublicUploadFiles(upload_id, **kwargs)


class TestPublicUploadFiles(UploadFilesContract):

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id: str) -> UploadFiles:
        _, upload_files = create_public_upload(
            test_upload_id, calc_specs='', is_authorized=lambda: True)

        return upload_files

    @pytest.fixture(scope='function', params=itertools.product(
        ['r', 'rr', 'pr', 'rp', 'p', 'pp', 'RP', 'RR', 'PP'], [True, False]))
    def test_upload(self, request, test_upload_id: str) -> PublicUploadWithFiles:
        calc_specs, protected = request.param
        upload, upload_files = create_staging_upload(test_upload_id, calc_specs=calc_specs)
        upload_files.pack(upload)
        upload_files.delete()
        return upload, PublicUploadFiles(test_upload_id, is_authorized=lambda: not protected)

    def test_to_staging_upload_files(self, test_upload):
        upload, upload_files = test_upload
        assert upload_files.to_staging_upload_files() is None
        staging_upload_files = upload_files.to_staging_upload_files(create=True)
        assert staging_upload_files is not None
        assert str(staging_upload_files) == str(upload_files.to_staging_upload_files())

        upload_path = upload_files.os_path
        all_files = list(
            os.path.join(upload_path, f)
            for f in os.listdir(upload_path)
            if os.path.isfile(os.path.join(upload_path, f)))

        # We override the public files before packing to see what packing does to the files
        for f in all_files:
            with open(f, 'wt') as fh:
                fh.write('')

        staging_upload_files.pack(upload)
        staging_upload_files.delete()

        # We do a very simple check. We made all files empty, those that are rezipped
        # by pack, should not be empty anymore.
        new_sizes = list(os.path.getsize(f) for f in all_files)
        for f, new in zip(all_files, new_sizes):
            if 'archive' in f:
                assert new > 0
            else:
                assert new == 0

        assert upload_files.to_staging_upload_files() is None

    def test_repack(self, test_upload):
        upload, upload_files = test_upload
        for calc in upload.calcs:
            calc.with_embargo = False
        upload_files.re_pack(upload)
        assert_upload_files(upload, PublicUploadFiles, with_embargo=False)
        assert len(os.listdir(upload_files.os_path)) == 4
        with assert_exception(KeyError):
            StagingUploadFiles(upload_files.upload_id)


def assert_upload_files(
        upload: UploadWithMetadata, cls, no_archive: bool = False, **kwargs):
    """
    Asserts the files aspect of uploaded data after processing or publishing

    Arguments:
        upload_id: The id of the upload to assert
        cls: The :class:`UploadFiles` subclass that this upload should have
        n_calcs: The number of expected calcs in the upload
        **kwargs: Key, value pairs that each calc metadata should have
    """
    upload_files = UploadFiles.get(upload.upload_id, is_authorized=lambda: True)
    assert upload_files is not None
    assert isinstance(upload_files, cls)

    upload_files = UploadFiles.get(upload.upload_id)
    for calc in upload.calcs:
        try:
            with upload_files.raw_file(calc.mainfile) as f:
                f.read()

            try:
                with upload_files.archive_file(calc.calc_id) as f:
                    f.read()
                with upload_files.archive_log_file(calc.calc_id) as f:
                    f.read()
            except KeyError:
                assert no_archive

            assert not calc.with_embargo and isinstance(upload_files, PublicUploadFiles)
        except Restricted:
            assert calc.with_embargo or isinstance(upload_files, StagingUploadFiles)
