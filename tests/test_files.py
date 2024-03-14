#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

from typing import Generator, Any, Dict, Tuple, Iterable, List
from datetime import datetime
import os
import os.path
import shutil
import pytest
import itertools
import zipfile
import re
import pathlib

from nomad import config, datamodel, utils
from nomad.archive import to_json
from nomad.files import (
    DirectoryObject,
    PathObject,
    empty_zip_file_size,
    empty_archive_file_size,
)
from nomad.files import StagingUploadFiles, PublicUploadFiles, UploadFiles
from nomad.processing import Upload


EntryWithFiles = Tuple[datamodel.EntryMetadata, str]
UploadWithFiles = Tuple[str, List[datamodel.EntryMetadata], UploadFiles]
StagingUploadWithFiles = Tuple[str, List[datamodel.EntryMetadata], StagingUploadFiles]
PublicUploadWithFiles = Tuple[str, List[datamodel.EntryMetadata], PublicUploadFiles]

# example_file uses an artificial parser for faster test execution, can also be
# changed to examples_vasp.zip for using vasp parser
example_mainfile_raw_path = 'examples_template/template.json'

example_file = 'tests/data/proc/examples_template.zip'
example_directory = 'tests/data/proc/examples_template'
example_file_contents = [
    'examples_template/template.json',
    'examples_template/1.aux',
    'examples_template/2.aux',
    'examples_template/3.aux',
    'examples_template/4.aux',
]
example_file_aux = 'tests/data/proc/examples_template/1.aux'
example_file_mainfile = 'tests/data/proc/examples_template/template.json'
example_file_mainfile_different_atoms = (
    'tests/data/proc/templates/different_atoms/template.json'
)
example_file_unparsable = 'tests/data/proc/templates/unparsable/template.json'
example_file_vasp_with_binary = 'tests/data/proc/example_vasp_with_binary.zip'
example_file_corrupt_zip = 'tests/data/proc/examples_corrupt_zip.zip'
empty_file = 'tests/data/proc/empty.zip'
example_archive_contents = {
    'run': [],
    'metadata': {},
    'processing_logs': [{'entry': 'test'}],
}


@pytest.fixture(scope='function', autouse=True)
def raw_files_on_all_tests(raw_files_function):
    """Autouse fixture to apply raw_files to all tests."""
    pass


@pytest.fixture(scope='session')
def example_mainfile_contents():
    with zipfile.ZipFile(example_file, 'r') as zf:
        with zf.open(example_mainfile_raw_path) as f:
            return f.read().decode()


class TestObjects:
    @pytest.fixture(scope='function')
    def test_area(self):
        yield config.fs.staging

        if os.path.exists(config.fs.staging):
            shutil.rmtree(config.fs.staging)

    def test_file_dir_existing(self, test_area):
        file = PathObject(os.path.join(test_area, 'sub/test_id'))
        assert not os.path.exists(os.path.dirname(file.os_path))

    @pytest.mark.parametrize('dirpath', ['test', os.path.join('sub', 'test')])
    @pytest.mark.parametrize('create', [True, False])
    def test_directory(self, test_area: str, dirpath: str, create: bool) -> None:
        directory = DirectoryObject(os.path.join(test_area, dirpath), create=create)
        assert directory.exists() == create
        assert os.path.isdir(directory.os_path) == create

    @pytest.mark.parametrize('dirpath', ['test', os.path.join('sub', 'test')])
    @pytest.mark.parametrize('create', [True, False])
    @pytest.mark.parametrize('join_create', [True, False])
    def test_directory_join(
        self, test_area: str, dirpath: str, create: bool, join_create: bool
    ) -> None:
        directory = DirectoryObject(os.path.join(test_area, 'parent'), create=create)
        directory = directory.join_dir(dirpath, create=join_create)

        assert directory.exists() == join_create
        assert os.path.isdir(directory.os_path) == join_create

    @pytest.mark.parametrize('filepath', ['test', 'sub/test'])
    @pytest.mark.parametrize('create', [True, False])
    def test_directory_join_file_dir_create(
        self, test_area: str, filepath: str, create: bool
    ):
        directory = DirectoryObject(os.path.join(test_area, 'parent'), create=create)
        file = directory.join_file(filepath, create_dir=create)
        assert os.path.exists(directory.os_path) == create
        assert os.path.exists(os.path.dirname(file.os_path)) == create


example_entry: Dict[str, Any] = {
    'entry_id': '0',
    'mainfile': 'examples_template/template.json',
    'data': 'value',
}
example_entry_id = example_entry['entry_id']


def generate_example_entry(
    entry_id: int, with_mainfile_prefix: bool, subdirectory: str = None, **kwargs
) -> EntryWithFiles:
    """Generate an example entry with :class:`EntryMetadata` and rawfile."""

    example_entry = datamodel.EntryMetadata(domain='dft', entry_id=str(entry_id))

    if with_mainfile_prefix:
        mainfile = '%d.template.json' % entry_id
    else:
        mainfile = 'template.json'

    if subdirectory is not None:
        mainfile = os.path.join(subdirectory, mainfile)

    example_entry.mainfile = mainfile
    example_entry.m_update(**kwargs)

    example_file = os.path.join(config.fs.tmp, 'example.zip')
    example_entry.files = []
    with zipfile.ZipFile(example_file, 'w', zipfile.ZIP_DEFLATED) as zf:
        for filepath in example_file_contents:
            filename = os.path.basename(filepath)
            arcname = filename
            if arcname == 'template.json' and with_mainfile_prefix:
                arcname = '%d.template.json' % entry_id

            if subdirectory is not None:
                arcname = os.path.join(subdirectory, arcname)
            example_entry.files.append(arcname)
            zf.write(os.path.join(example_directory, filename), arcname)

    return example_entry, example_file


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

    source = sorted(
        set(
            normalized_file(name)
            for name in names
            if not name.endswith('template.json') or with_mainfile or not is_multi
        )
    )
    target = sorted(
        name
        for name in example_file_contents
        if not name.endswith('template.json') or with_mainfile
    )
    assert source == target


def assert_example_entry(entry):
    assert entry is not None
    assert entry['data'] == example_entry['data']


class UploadFilesFixtures:
    @pytest.fixture(scope='function')
    def test_upload_id(self) -> Generator[str, None, None]:
        for cls in [StagingUploadFiles, PublicUploadFiles]:
            directory = DirectoryObject(cls.base_folder_for('test_upload'))
            if directory.exists():
                directory.delete()
        yield 'test_upload'
        for cls in [StagingUploadFiles, PublicUploadFiles]:
            directory = DirectoryObject(cls.base_folder_for('test_upload'))
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
        assert (
            UploadFiles.get(empty_test_upload.upload_id).__class__
            == empty_test_upload.__class__
        )

    def test_os_path(self, test_upload: UploadWithFiles):
        upload_files = test_upload[2]
        assert upload_files.os_path is not None
        if upload_files.external_os_path:
            os_posix_path = pathlib.Path(upload_files.os_path).as_posix()
            assert upload_files.external_os_path.endswith(os_posix_path)

    def test_rawfile(self, test_upload: UploadWithFiles):
        _, entries, upload_files = test_upload
        for entry in entries:
            for file_path in entry.files:
                with upload_files.raw_file(file_path) as f:
                    assert len(f.read()) > 0

    def test_rawfile_size(self, test_upload: UploadWithFiles):
        _, entries, upload_files = test_upload
        for entry in entries:
            for file_path in entry.files:
                assert upload_files.raw_file_size(file_path) > 0

    @pytest.mark.parametrize('prefix', [None, 'examples'])
    def test_raw_directory_list_prefix(self, test_upload: UploadWithFiles, prefix: str):
        _, _, upload_files = test_upload
        path_infos = upload_files.raw_directory_list(
            recursive=True, files_only=True, path_prefix=prefix
        )
        raw_files = list(path_info.path for path_info in path_infos)
        assert_example_files(raw_files)

    @pytest.mark.parametrize('path', ['', 'examples_template'])
    def test_raw_directory_list(self, test_upload: UploadWithFiles, path: str):
        upload_id, _, upload_files = test_upload
        # Add file to root to test corner case
        append_raw_files(upload_id, 'tests/data/proc/examples_template/1.aux', '1.aux')
        # Test recursive call (but do not verify result)
        upload_files.raw_directory_list(path, files_only=False, recursive=True)
        # Test non-recursive call, verify result partially
        raw_files = list(upload_files.raw_directory_list(path, files_only=True))
        if not path:
            assert len(raw_files) == 1
            assert raw_files[0].size == 8
        else:
            assert '1.aux' in list(
                os.path.basename(path_info.path) for path_info in raw_files
            )
            for path_info in raw_files:
                if path_info.path.endswith('.aux'):
                    assert path_info.size == 8
                else:
                    assert path_info.size > 0
            assert_example_files([path_info.path for path_info in raw_files])

    @pytest.mark.parametrize('with_access', [False, True])
    def test_read_archive(self, test_upload: UploadWithFiles, with_access: str):
        _, _, upload_files = test_upload

        with upload_files.read_archive(example_entry_id) as archive:
            assert to_json(archive[example_entry_id]) == example_archive_contents


def create_staging_upload(
    upload_id: str, entry_specs: str, embargo_length: int = 0
) -> StagingUploadWithFiles:
    """
    Create an upload according to given spec. Additional arguments are given to
    the StagingUploadFiles contstructor.

    Arguments:
        upload_id: The id that should be given to this test upload.
        entry_specs: A string that determines the properties of the given upload.
            With letters determining example entries being public `p` or restricted `r`.
            The entries will be copies of entries in `example_file`.
            First entry is at top level, following entries will be put under 1/, 2/, etc.
            All entries with capital `P`/`R` will be put in the same directory under multi/.
    """
    upload_files = StagingUploadFiles(upload_id, create=True)
    entries = []

    prefix = 0
    for entry_spec in entry_specs:
        is_multi = entry_spec in ['R', 'P']
        entry_spec = entry_spec.lower()
        assert (entry_spec == 'r') == (embargo_length > 0)
        if is_multi or prefix == 0:
            directory = 'examples_template'
        else:
            directory = os.path.join(str(prefix), 'examples_template')

        entry, entry_file = generate_example_entry(
            prefix,
            with_mainfile_prefix=is_multi,
            subdirectory=directory,
            with_embargo=embargo_length > 0,
        )

        upload_files.add_rawfiles(entry_file)
        upload_files.write_archive(entry.entry_id, example_archive_contents)

        entries.append(entry)
        prefix += 1

    assert len(entries) == len(entry_specs)
    return upload_id, entries, upload_files


class TestStagingUploadFiles(UploadFilesContract):
    @pytest.fixture(scope='function', params=['r', 'rr', 'p', 'pp', 'RR', 'PP'])
    def test_upload(self, request, test_upload_id: str) -> StagingUploadWithFiles:
        embargo_length = 12 if 'r' in request.param.lower() else 0
        return create_staging_upload(
            test_upload_id, entry_specs=request.param, embargo_length=embargo_length
        )

    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id) -> UploadFiles:
        return StagingUploadFiles(test_upload_id, create=True)

    @pytest.mark.parametrize('target_dir', ['', 'subdir'])
    def test_add_rawfiles_zip(self, test_upload_id, target_dir):
        test_upload = StagingUploadFiles(test_upload_id, create=True)
        test_upload.add_rawfiles(example_file, target_dir=target_dir)
        for filepath in example_file_contents:
            filepath = os.path.join(target_dir, filepath) if target_dir else filepath
            with test_upload.raw_file(filepath) as f:
                content = f.read()
                if filepath == example_mainfile_raw_path:
                    assert len(content) > 0

    def test_pack(self, test_upload: StagingUploadWithFiles):
        _, entries, upload_files = test_upload
        upload_files.pack(entries, with_embargo=entries[0].with_embargo)

    @pytest.mark.parametrize('entry_specs', ['r', 'p'])
    def test_pack_potcar(self, entry_specs):
        embargo_length = 12 if 'r' in entry_specs.lower() else 0
        upload_id, entries, upload_files = create_staging_upload(
            'test_potcar', entry_specs=entry_specs, embargo_length=embargo_length
        )
        # Add potcar files: one stripped and one unstripped
        filenames = ('POTCAR', 'POTCAR.stripped')
        for filename in filenames:
            with open(
                os.path.join(
                    upload_files.os_path, 'raw', 'examples_template', filename
                ),
                'w',
            ) as f:
                f.write('some content')
        upload_files.pack(entries, with_embargo=embargo_length > 0)
        upload_files.delete()
        upload_files = PublicUploadFiles(upload_id)
        for filename in filenames:
            try:
                pf = upload_files.raw_file('examples_template/' + filename)
                pf.read()
                assert filename.endswith(
                    '.stripped'
                ), 'Non-stripped POTCAR file could be read'
            except KeyError:
                assert not filename.endswith(
                    '.stripped'
                ), 'Only non-stripped file should be removed'

    @pytest.mark.parametrize('with_mainfile', [True, False])
    def test_entry_files(self, test_upload: StagingUploadWithFiles, with_mainfile):
        _, entries, upload_files = test_upload
        for entry in entries:
            mainfile = entry.mainfile
            entry_files = upload_files.entry_files(
                mainfile, with_mainfile=with_mainfile
            )
            assert_example_files(entry_files, with_mainfile=with_mainfile)

    def test_delete(self, test_upload: StagingUploadWithFiles):
        _, _, upload_files = test_upload
        upload_files.delete()
        assert not upload_files.exists()

    def test_add_rawfiles(self, test_upload_id):
        test_upload = StagingUploadFiles(test_upload_id, create=True)
        assert test_upload.is_empty()
        test_upload.add_rawfiles(example_file)
        path_infos = test_upload.raw_directory_list(recursive=True, files_only=True)
        assert sorted(list(path_info.path for path_info in path_infos)) == sorted(
            example_file_contents
        )

    @pytest.mark.parametrize('prefix_size', [0, 2])
    def test_prefix_size(self, monkeypatch, prefix_size):
        monkeypatch.setattr('nomad.config.fs.prefix_size', prefix_size)
        upload_id = 'test_upload'
        upload_files = StagingUploadFiles(upload_id, create=True)
        if not prefix_size:
            assert upload_files.os_path == os.path.join(config.fs.staging, upload_id)
        else:
            prefix = upload_id[:prefix_size]
            assert upload_files.os_path == os.path.join(
                config.fs.staging, prefix, upload_id
            )
        upload_files.delete()

    def test_delete_prefix(self, monkeypatch):
        monkeypatch.setattr('nomad.config.fs.prefix_size', 2)
        upload_1 = StagingUploadFiles('test_upload_1', create=True)
        upload_2 = StagingUploadFiles('test_upload_2', create=True)
        upload_1.delete()
        upload_2.delete()

        prefix = os.path.dirname(upload_1.os_path)
        assert len(os.path.basename(prefix)) == 2
        assert not os.path.exists(prefix)


def create_public_upload(
    upload_id: str, entry_specs: str, embargo_length: int = 0, with_upload: bool = True
) -> PublicUploadWithFiles:
    _, entries, upload_files = create_staging_upload(
        upload_id, entry_specs, embargo_length
    )

    upload_files.pack(entries, with_embargo=embargo_length > 0)
    upload_files.delete()
    if with_upload:
        upload = Upload.get(upload_id)
        upload.publish_time = datetime.utcnow()
        assert upload.embargo_length == embargo_length, 'Wrong embargo_length provided'
        upload.save()
    return upload_id, entries, PublicUploadFiles(upload_id)


class TestPublicUploadFiles(UploadFilesContract):
    @pytest.fixture(scope='function')
    def empty_test_upload(self, test_upload_id: str) -> UploadFiles:
        _, _, upload_files = create_public_upload(
            test_upload_id, entry_specs='', with_upload=False
        )

        return upload_files

    @pytest.fixture(
        scope='function',
        params=itertools.product(['r', 'rr', 'p', 'pp', 'RR', 'PP'], [True, False]),
    )
    def test_upload(self, request, test_upload_id: str) -> PublicUploadWithFiles:
        entry_specs, both_accesses = request.param
        embargo_length = 12 if 'r' in entry_specs.lower() else 0
        _, entries, upload_files = create_staging_upload(
            test_upload_id, entry_specs=entry_specs, embargo_length=embargo_length
        )
        upload_files.pack(entries, with_embargo=embargo_length > 0)
        upload_files.delete()
        public_upload_files = PublicUploadFiles(test_upload_id)
        if both_accesses:
            # Artificially create an empty archive files and raw zip file with the opposite access
            # TODO: This should only be needed for an interim period
            other_access = 'public' if embargo_length else 'restricted'
            other_raw_zip_file_object = PublicUploadFiles._create_raw_zip_file_object(
                public_upload_files, other_access
            )
            other_msg_file_object = PublicUploadFiles._create_msg_file_object(
                public_upload_files, other_access
            )
            # Fill them with dummy content (we should never try to open them)
            with open(other_raw_zip_file_object.os_path, mode='wb') as f:
                f.write(b'-' * empty_zip_file_size)
            with open(other_msg_file_object.os_path, mode='wb') as f:
                f.write(b'-' * empty_archive_file_size)
        return test_upload_id, entries, PublicUploadFiles(test_upload_id)

    def test_to_staging_upload_files(self, test_upload):
        _, entries, upload_files = test_upload
        access = upload_files.access
        assert upload_files.to_staging_upload_files() is None
        staging_upload_files = upload_files.to_staging_upload_files(create=True)
        assert staging_upload_files is not None
        assert str(staging_upload_files) == str(upload_files.to_staging_upload_files())

        upload_path = upload_files.os_path
        all_files = list(
            os.path.join(upload_path, f)
            for f in os.listdir(upload_path)
            if os.path.isfile(os.path.join(upload_path, f))
        )

        # We override the public files before packing to see what packing does to the files
        for f in all_files:
            with open(f, 'wb') as fh:
                if access in os.path.basename(f):
                    fh.write(b'-' * 50)
                else:
                    fh.write(b'')

        staging_upload_files.pack(
            entries,
            with_embargo=entries[0].with_embargo,
            create=False,
            include_raw=False,
        )
        staging_upload_files.delete()

        # We do a very simple check. Files that are expected to be modified by pack
        # should have a different size, and files that have been removed should have the wrong access.
        for file_name in all_files:
            if not os.path.exists(file_name):
                assert access not in os.path.basename(file_name)
            elif access in os.path.basename(file_name):
                if 'archive' in file_name:
                    assert (
                        os.path.getsize(file_name) > 100
                    ), 'Archive files should have been packed'
                else:
                    assert (
                        os.path.getsize(file_name) == 50
                    ), 'Raw file should not have been changed'
            else:
                # other access
                assert (
                    os.path.getsize(file_name) <= empty_archive_file_size
                ), 'Files with other access should be empty'

        assert upload_files.to_staging_upload_files() is None

    def test_repack(self, test_upload):
        upload_id, entries, upload_files = test_upload
        for entry in entries:
            entry.with_embargo = False
        upload_files.re_pack(with_embargo=False)
        assert_upload_files(upload_id, entries, PublicUploadFiles, with_embargo=False)
        assert upload_files.access == 'public'
        with pytest.raises(KeyError):
            StagingUploadFiles(upload_files.upload_id)

    @pytest.mark.parametrize(
        'suffixes,suffix',
        [
            pytest.param(None, '', id='none'),
            pytest.param('v1', '-v1', id='single'),
            pytest.param(['v2', 'v1'], '-v2', id='fallback'),
        ],
    )
    def test_archive_version_suffix(
        self, monkeypatch, test_upload_id, suffixes, suffix
    ):
        monkeypatch.setattr('nomad.config.fs.archive_version_suffix', suffixes)
        _, entries, upload_files = create_staging_upload(
            test_upload_id, entry_specs='p'
        )
        upload_files.pack(entries, with_embargo=False)
        upload_files.delete()

        public_upload_files = PublicUploadFiles(test_upload_id)

        assert os.path.exists(
            public_upload_files.join_file('raw-public.plain.zip').os_path
        )
        assert os.path.exists(
            public_upload_files.join_file(f'archive-public{suffix}.msg.msg').os_path
        )

        assert_upload_files(test_upload_id, entries, PublicUploadFiles)

    def test_archive_version_suffix_fallback(self, monkeypatch, test_upload_id):
        monkeypatch.setattr('nomad.config.fs.archive_version_suffix', ['v1'])
        _, entries, upload_files = create_staging_upload(
            test_upload_id, entry_specs='p'
        )
        upload_files.pack(entries, with_embargo=False)

        monkeypatch.setattr('nomad.config.fs.archive_version_suffix', ['v2', 'v1'])
        v1_file = upload_files._archive_file_object(0, fallback=True)
        v2_file = upload_files._archive_file_object(0, fallback=False)
        assert os.path.basename(v1_file.os_path) == '0-v1.msg'
        assert os.path.basename(v2_file.os_path) == '0-v2.msg'

        upload_files.write_archive('0', {})
        assert v1_file.exists()
        assert v2_file.exists()
        assert (
            os.path.basename(
                upload_files._archive_file_object(0, fallback=True).os_path
            )
            == '0-v2.msg'
        )

        upload_files.delete()

        public_upload_files = PublicUploadFiles(test_upload_id)
        v2_file = PublicUploadFiles._create_msg_file_object(
            public_upload_files, public_upload_files.access, fallback=False
        )
        v1_file = PublicUploadFiles._create_msg_file_object(
            public_upload_files, public_upload_files.access, fallback=True
        )
        assert not v2_file.exists()
        assert os.path.basename(v1_file.os_path) == 'archive-public-v1.msg.msg'


def assert_upload_files(
    upload_id: str,
    entries: Iterable[datamodel.EntryMetadata],
    cls,
    no_archive: bool = False,
    **kwargs,
):
    """
    Asserts the files aspect of uploaded data after processing or publishing

    Arguments:
        upload_id: The id of the upload to assert
        cls: The :class:`UploadFiles` subclass that this upload should have
        no_archive:
        **kwargs: Key, value pairs that each entry metadata should have
    """
    upload_files = UploadFiles.get(upload_id)
    assert upload_files is not None
    assert isinstance(upload_files, cls)

    upload_files = UploadFiles.get(upload_id)
    for entry in entries:
        with upload_files.raw_file(entry.mainfile, 'rb') as f:
            f.read()

        try:
            with upload_files.read_archive(entry.entry_id) as archive:
                assert entry.entry_id in archive
        except KeyError:
            assert no_archive

    upload_files.close()


def create_test_upload_files(
    upload_id: str,
    archives: List[datamodel.EntryArchive] = None,
    published: bool = True,
    embargo_length: int = 0,
    raw_files: str = None,
    template_files: str = example_file,
    template_mainfile: str = example_mainfile_raw_path,
    additional_files_path: str = None,
) -> UploadFiles:
    """
    Creates an upload_files object and the underlying files for test/mock purposes.

    Arguments:
        upload_id: The upload id for the upload. Will generate a random UUID if None.
        archives: A list of class:`datamodel.EntryArchive` metainfo objects. This will
            be used to determine the mainfiles. Will create respective directories and
            copy the template entry to create raw files for each archive.
            Will also be used to fill the archives in the create upload.
        published: Creates a :class:`PublicUploadFiles` object with published files
            instead of a :class:`StagingUploadFiles` object with staging files. Default
            is published.
        embargo_length: The embargo length
        raw_files: A directory path. All files here will be copied into the raw files
            dir of the created upload files.
        template_files: A zip file with example files in it. One directory will be used
            as a template. It will be copied for each given archive.
        template_mainfile: Path of the template mainfile within the given template_files.
        additional_files_path: Path to additional files to add.
    """
    if upload_id is None:
        upload_id = utils.create_uuid()
    if archives is None:
        archives = []

    upload_files = StagingUploadFiles(upload_id, create=True)
    if raw_files:
        shutil.rmtree(upload_files._raw_dir.os_path)
        shutil.copytree(raw_files, upload_files._raw_dir.os_path)
    upload_files.add_rawfiles(template_files)
    if additional_files_path:
        upload_files.add_rawfiles(additional_files_path)

    upload_raw_files = upload_files.join_dir('raw')
    source = upload_raw_files.join_dir(os.path.dirname(template_mainfile)).os_path

    if archives is None:
        archives = []
    for archive in archives:
        # create a copy of the given template files for each archive
        mainfile = archive.metadata.mainfile
        mainfile_key = archive.metadata.mainfile_key
        assert (
            mainfile is not None
        ), 'Archives to create test upload must have a mainfile'
        target = upload_raw_files.join_file(os.path.dirname(mainfile)).os_path
        if not mainfile_key:
            if os.path.exists(target):
                for file_ in os.listdir(source):
                    shutil.copy(os.path.join(source, file_), target)
            else:
                shutil.copytree(source, target)
            os.rename(
                os.path.join(target, os.path.basename(template_mainfile)),
                os.path.join(target, os.path.basename(mainfile)),
            )

        # create an archive "file" for each archive
        entry_id = archive.metadata.entry_id
        assert (
            entry_id is not None
        ), 'Archives to create test upload must have an entry_id'
        upload_files.write_archive(entry_id, archive.m_to_dict())

    # remove the template
    shutil.rmtree(source)

    if published:
        upload_files.pack(
            [archive.metadata for archive in archives], with_embargo=embargo_length > 0
        )
        upload_files.delete()
        return UploadFiles.get(upload_id)

    return upload_files


def append_raw_files(upload_id: str, path_source: str, path_in_upload: str):
    """
    Used to append files to the raw files of an upload (published or not), for
    testing purposes.
    """
    upload_files = UploadFiles.get(upload_id)
    if isinstance(upload_files, PublicUploadFiles):
        zip_path = upload_files.raw_zip_file_object().os_path  # type: ignore
        with zipfile.ZipFile(zip_path, 'a') as zf:
            zf.write(path_source, path_in_upload)
    else:
        path = upload_files._raw_dir.os_path  # type: ignore
        shutil.copy(path_source, os.path.join(path, path_in_upload))


def test_test_upload_files(raw_files_infra):
    upload_id = utils.create_uuid()
    archives: datamodel.EntryArchive = []
    for index in range(0, 3):
        archive = datamodel.EntryArchive()
        metadata = archive.m_create(datamodel.EntryMetadata)
        metadata.entry_id = 'example_entry_id_%d' % index
        metadata.mainfile = 'test/test/entry_%d/mainfile_%d.json' % (index, index)
        archives.append(archive)

    upload_files = create_test_upload_files(upload_id, archives, embargo_length=0)

    try:
        assert_upload_files(
            upload_id, [archive.metadata for archive in archives], PublicUploadFiles
        )
    finally:
        if upload_files.exists():
            upload_files.delete()
