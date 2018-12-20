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

"""
Uploads contains classes and functions to create and maintain file structures
for uploads.

There are two different structures for uploads in two different states: *staging* and *public*.
Possible operations on uploads differ based on this state. Staging is used for
processing, heavily editing, creating hashes, etc. Public is supposed to be a
almost readonly (beside metadata) storage.

::
    fs/staging/<upload>/metadata/<calc>.json
                       /raw/**
                       /archive/<calc>.hdf5
                       /.frozen
                       /.public
                       /.restricted
    fs/public/<upload>/metadata.json
                      /metadata.json.lock
                      /raw-public.bagit.zip
                      /raw-restricted.bagit.zip
                      /archive-public.hdf5.zip
                      /archive-restricted.hdf5.zip
"""

from abc import ABCMeta
from typing import IO, Generator, Dict, Iterator, Iterable, Callable
from filelock import Timeout, FileLock
import ujson
import os.path
import os
import shutil
from zipfile import ZipFile, BadZipFile
from bagit import make_bag
import contextlib
import hashlib

from nomad import config, utils


class PathObject:
    """
    Object storage-like abstraction for paths in general.
    Arguments:
        bucket: The bucket to store this object in
        object_id: The object id (i.e. directory path)
        os_path: Override the "object storage" path with the given path.
        prefix: Add a 3-digit prefix directory, e.g. foo/test/ -> foo/tes/test
    """
    def __init__(self, bucket: str, object_id: str, os_path: str = None, prefix: bool = False) -> None:
        if os_path:
            self.os_path = os_path
        else:
            self.os_path = os.path.join(config.fs.objects, bucket, object_id)

        if prefix:
            segments = list(os.path.split(self.os_path))
            last = segments[-1]
            segments[-1] = last[:3]
            segments.append(last)
            self.os_path = os.path.join(*segments)

    def delete(self) -> None:
        shutil.rmtree(self.os_path)

    def exists(self) -> bool:
        return os.path.exists(self.os_path)

    def __repr__(self) -> str:
        return self.os_path


class DirectoryObject(PathObject):
    """
    Object storage-like abstraction for directories.
    Arguments:
        bucket: The bucket to store this object in
        object_id: The object id (i.e. directory path)
        create: True if the directory structure should be created. Default is False.
    """
    def __init__(self, bucket: str, object_id: str, create: bool = False, **kwargs) -> None:
        super().__init__(bucket, object_id, **kwargs)
        self._create = create
        if create and not os.path.isdir(self.os_path):
            os.makedirs(self.os_path)

    def join_dir(self, path, create: bool = None) -> 'DirectoryObject':
        if create is None:
            create = self._create
        return DirectoryObject(None, None, create=create, os_path=os.path.join(self.os_path, path))

    def join_file(self, path) -> PathObject:
        dirname = os.path.dirname(path)
        if dirname != '':
            return self.join_dir(dirname).join_file(os.path.basename(path))
        else:
            return PathObject(None, None, os_path=os.path.join(self.os_path, path))

    def exists(self) -> bool:
        return os.path.isdir(self.os_path)


class MetadataTimeout(Exception):
    pass


class Metadata(metaclass=ABCMeta):
    """
    An ABC for a contextmanager that encapsulates access to a set of calc metadata.
    Allows to add, update, read metadata. Subclasses might deal with concurrent access.
    """
    def __enter__(self) -> 'Metadata':
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        return None

    def open(self):
        pass

    def close(self):
        pass

    def insert(self, calc: dict) -> None:
        """ Insert a calc, using hash as key. """
        raise NotImplementedError()

    def update(self, calc_hash: str, updates: dict) -> dict:
        """ Updating a calc, using hash as key and running dict update with the given data. """
        raise NotImplementedError()

    def get(self, calc_id: str) -> dict:
        """ Retrive the calc metadata for a given calc. """
        raise NotImplementedError()

    def __iter__(self) -> Iterator[dict]:
        raise NotImplementedError()

    def __len__(self) -> int:
        raise NotImplementedError()


class StagingMetadata(Metadata):
    """
    A Metadata implementation based on individual .json files per calc stored in a given
    directory.
    Arguments:
        directory: The DirectoryObject for the directory to store the metadata in.
    """
    def __init__(self, directory: DirectoryObject) -> None:
        self._dir = directory

    def __enter__(self) -> 'Metadata':
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        return None

    def open(self):
        pass

    def close(self):
        pass

    def insert(self, calc: dict) -> None:
        id = calc['hash']
        path = self._dir.join_file('%s.json' % id)
        assert not path.exists()
        with open(path.os_path, 'wt') as f:
            ujson.dump(calc, f)

    def update(self, calc_hash: str, updates: dict) -> dict:
        metadata = self.get(calc_hash)
        metadata.update(updates)
        path = self._dir.join_file('%s.json' % calc_hash)
        with open(path.os_path, 'wt') as f:
            ujson.dump(metadata, f)
        return metadata

    def get(self, calc_id: str) -> dict:
        try:
            with open(self._dir.join_file('%s.json' % calc_id).os_path, 'rt') as f:
                return ujson.load(f)
        except FileNotFoundError:
            raise KeyError()

    def __iter__(self) -> Iterator[dict]:
        for root, _, files in os.walk(self._dir.os_path):
            for file in files:
                with open(os.path.join(root, file), 'rt') as f:
                    yield ujson.load(f)

    def __len__(self) -> int:
        return len(os.listdir(self._dir.os_path))


class PublicMetadata(Metadata):
    """
    A Metadata implementation based on a single .json file. It loads and write
    the metadata to the given path and uses a lock to deal with concurrent access.

    Arguments:
        path: The parent directory for the metadata and lock file.
        lock_timeout: Max timeout before __enter__ raises MetadataTimeout while waiting
            for an available lock on the metadata file. Default is 1s.
    """
    def __init__(self, path: str, lock_timeout=1) -> None:
        self._db_file = os.path.join(path, 'metadata.json')
        self._lock_file = os.path.join(path, 'metadata.json.lock')
        self._lock: FileLock = FileLock(self._lock_file, timeout=lock_timeout)
        self._modified = False
        self.data: Dict[str, dict] = None

    def __enter__(self) -> 'Metadata':
        self.open()
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.close()
        return None

    def open(self):
        assert self.data is None, "Metadata is already open."

        try:
            self._lock.acquire()
        except Timeout:
            raise MetadataTimeout()

        if os.path.exists(self._db_file):
            with open(self._db_file, 'rt') as f:
                self.data = ujson.load(f)
        else:
            self.data = {}
            self._modified = True

    def close(self):
        assert self.data is not None, "Metadata is not open."
        if self._modified:
            with open(self._db_file, 'wt') as f:
                ujson.dump(self.data, f, ensure_ascii=False)
        self.data = None
        self._lock.release()

    def insert(self, calc: dict) -> None:
        assert self.data is not None, "Metadata is not open."

        id = calc['hash']
        assert id not in self.data
        self.data[id] = calc
        self._modified = True

    def update(self, calc_hash: str, updates: dict) -> dict:
        assert self.data is not None, "Metadata is not open."
        if calc_hash not in self.data:
            raise KeyError()

        self.data[calc_hash].update(updates)
        self._modified = True
        return self.data[calc_hash]

    def get(self, calc_hash: str) -> dict:
        assert self.data is not None, "Metadata is not open."

        return self.data[calc_hash]

    def __iter__(self) -> Iterator[dict]:
        assert self.data is not None, "Metadata is not open."
        return self.data.values().__iter__()

    def __len__(self) -> int:
        assert self.data is not None, "Metadata is not open."
        return len(self.data)


class Restricted(Exception):
    pass


class UploadFiles(metaclass=ABCMeta):
    def __init__(self, upload_id: str, archive_ext: str = 'json') -> None:
        self.logger = utils.get_logger(__name__, upload_id=upload_id)
        self.upload_id = upload_id
        self._archive_ext = archive_ext

    @property
    def metadata(self) -> Metadata:
        """ The calc metadata for this upload. """
        raise NotImplementedError

    @contextlib.contextmanager
    def raw_file(self, file_path: str, read: bool = True, access: Callable[[dict], bool] = None) -> Generator[IO, None, None]:
        """
        Opens a raw file and returns a file-like objects.
        Arguments:
            file_path: The path to the file relative to the upload.
            read: Open for read or write. Default is True=read.
            access: Function that evaluates calc metadata to bool and determines
                restricted data might be accessed.
        Raises:
            KeyError: If the file does not exist.
            Restricted: If the file is restricted and access function is not given or evaluated to False.
        """
        raise NotImplementedError()

    @contextlib.contextmanager
    def archive_file(self, calc_hash: str, read: bool = True) -> Generator[IO, None, None]:
        """
        Opens a archive file and returns a file-like objects.
        Arguments:
            calc_hash: The hash identifying the calculation.
            read: Open for read or write. Default is True=read.
        Raises:
            KeyError: If the calc does not exist.
        """
        raise NotImplementedError()


class StagingUploadFiles(UploadFiles):
    def __init__(self, upload_id: str, create: bool = False, archive_ext: str = 'json') -> None:
        super().__init__(upload_id=upload_id, archive_ext=archive_ext)

        self._upload_dir = DirectoryObject(
            config.files.staging_bucket, upload_id, create=create, prefix=True)
        if not create and not self._upload_dir.exists():
            raise KeyError()
        self._raw_dir = self._upload_dir.join_dir('raw')
        self._archive_dir = self._upload_dir.join_dir('archive')
        self._frozen_file = self._upload_dir.join_file('.frozen')
        self._restricted_dir = self._upload_dir.join_dir('.restricted', create=False)
        self._public_dir = self._upload_dir.join_dir('.public', create=False)

        metadata_dir = self._upload_dir.join_dir('metadata')
        self._metadata = StagingMetadata(metadata_dir)

    @property
    def metadata(self) -> Metadata:
        return self._metadata

    @contextlib.contextmanager
    def _file(self, path, read: bool) -> Generator[IO, None, None]:
        try:
            with open(path, 'rb' if read else 'wb') as f:
                yield f
        except FileNotFoundError:
            raise KeyError()

    @contextlib.contextmanager
    def raw_file(self, file_path: str, read: bool = True) -> Generator[IO, None, None]:
        path = os.path.join(self._raw_dir.os_path, file_path)
        with self._file(path, read) as f:
            yield f

    @contextlib.contextmanager
    def archive_file(self, calc_hash: str, read: bool = True) -> Generator[IO, None, None]:
        path = os.path.join(self._archive_dir.os_path, '%s.%s' % (calc_hash, self._archive_ext))
        with self._file(path, read) as f:
            yield f

    def add_rawfiles(self, path: str, move: bool = False, prefix: str = None) -> None:
        """
        Add rawfiles to the upload. The given file will be copied, moved, or extracted.
        Arguments:
            path: Path to a directory, file, or zip file. Zip files will be extracted.
            move: Whether the file should be moved instead of copied. Zips will be extracted and then deleted.
            prefix: Optional path prefix for the added files.
        """
        assert not self.is_frozen
        assert os.path.exists(path)
        target_dir = self._raw_dir
        if prefix is not None:
            target_dir = target_dir.join_dir(prefix, create=True)
        ext = os.path.splitext(path)[1]
        if ext == '.zip':
            try:
                with ZipFile(path) as zf:
                    zf.extractall(target_dir.os_path)
                if move:
                    os.remove(path)
                return
            except BadZipFile:
                pass

        if move:
            shutil.move(path, target_dir.os_path)
        else:
            shutil.copy(path, target_dir.os_path)

    @property
    def is_frozen(self) -> bool:
        """ Returns True if this upload is already *bagged*. """
        return self._frozen_file.exists()

    def pack(self, bagit_metadata: dict = None) -> None:
        """
        Replaces the staging upload data with a public upload record by packing all
        data into files. It is only available if upload *is_bag*.
        This is potentially a long running operation.
        Arguments:
            bagit_metadata: Additional data added to the bagit metadata.
        """
        # freeze the upload
        assert not self.is_frozen, "Cannot pack an upload that is packed, or packing."
        with open(self._frozen_file.os_path, 'wt') as f:
            f.write('frozen')
        os.makedirs(self._public_dir.os_path)

        # copy raw -> .restricted
        shutil.copytree(self._raw_dir.os_path, self._restricted_dir.os_path)

        # move public data .restricted -> .public
        for calc in self.metadata:
            if not calc.get('restricted', False):
                mainfile: str = calc['mainfile']
                assert mainfile is not None
                for filepath in self.calc_files(mainfile):
                    os.rename(
                        self._restricted_dir.join_file(filepath).os_path,
                        self._public_dir.join_file(filepath).os_path)

        # create bags
        make_bag(self._restricted_dir.os_path, bag_info=bagit_metadata, checksums=['sha512'])
        make_bag(self._public_dir.os_path, bag_info=bagit_metadata, checksums=['sha512'])

        # zip bags
        def zip_dir(zip_filepath, path):
            root_len = len(path)
            with ZipFile(zip_filepath, 'w') as zf:
                for root, _, files in os.walk(path):
                    for file in files:
                        filepath = os.path.join(root, file)
                        zf.write(filepath, filepath[root_len:])

        packed_dir = self._upload_dir.join_dir('.packed', create=True)

        zip_dir(packed_dir.join_file('raw-restricted.bagit.zip').os_path, self._restricted_dir.os_path)
        zip_dir(packed_dir.join_file('raw-public.bagit.zip').os_path, self._public_dir.os_path)

        # zip archives
        def create_zipfile(prefix: str) -> ZipFile:
            file = packed_dir.join_file('archive-%s.%s.zip' % (prefix, self._archive_ext))
            return ZipFile(file.os_path, mode='w')

        archive_public_zip = create_zipfile('public')
        archive_restricted_zip = create_zipfile('restricted')
        for calc in self.metadata:
            archive_filename = '%s.%s' % (calc['hash'], self._archive_ext)
            archive_zip = archive_restricted_zip if calc.get('restricted', False) else archive_public_zip
            archive_zip.write(self._archive_dir.join_file(archive_filename).os_path, archive_filename)

        # pack metadata
        with PublicMetadata(packed_dir.os_path) as packed_metadata:
            for calc in self.metadata:
                packed_metadata.insert(calc)

        # move to public bucket
        target_dir = DirectoryObject(config.files.public_bucket, self.upload_id, create=False, prefix=True)
        assert not target_dir.exists()
        shutil.move(packed_dir.os_path, target_dir.os_path)

    @property
    def all_rawfiles(self) -> Generator[str, None, None]:
        """ Returns: A generator of all file paths of all raw files. """
        for root, _, files in os.walk(self._raw_dir.os_path):
            for file in files:
                yield os.path.join(root, file)

    def calc_files(self, mainfile: str) -> Iterable[str]:
        """
        Returns all the auxfiles and mainfile for a given mainfile. This implements
        nomad's logic about what is part of a calculation and what not.
        """
        dir = os.path.dirname(self._raw_dir.join_file(mainfile).os_path)
        return sorted(path for path in os.listdir(dir) if os.path.isfile(path))

    def calc_hash(self, mainfile: str) -> str:
        """
        Calculates a hash for the given calc.
        Arguments:
            mainfile: The mainfile path relative to the upload that identifies the calc in the folder structure.
        Returns:
            The calc hash
        Raises:
            KeyError: If the mainfile does not exist.
        """
        hash = hashlib.sha512()
        for filepath in self.calc_files(mainfile):
            with open(filepath, 'rb') as f:
                for data in iter(lambda: f.read(65536), b''):
                    hash.update(data)

        return utils.websave_hash(hash.digest(), utils.default_hash_len)

    def upload_hash(self) -> str:
        """ Returns: A hash for the whole upload. It is only available if upload *is_bag*. """
        pass


class PublicUploadFiles(UploadFiles):
    def __init__(self, upload_id: str, *args, **kwargs) -> None:
        super().__init__(upload_id, *args, **kwargs)

        self._upload_dir = DirectoryObject(
            config.files.public_bucket, upload_id, create=False, prefix=True)
        self._metadata = PublicMetadata(self._upload_dir.os_path)

    @property
    def metadata(self) -> Metadata:
        return self._metadata

    @contextlib.contextmanager
    def _file(self, prefix: str, ext: str, path: str) -> Generator[IO, None, None]:
        for access in ['public', 'restricted']:
            try:
                zip_file = self._upload_dir.join_file('%s-%s.%s.zip' % (prefix, access, ext))
                with ZipFile(zip_file.os_path) as zf:
                    with zf.open(path, 'r') as f:
                        yield f
                        return
            except KeyError:
                pass

        raise KeyError()

    @contextlib.contextmanager
    def raw_file(self, file_path: str, read: bool = True) -> Generator[IO, None, None]:
        assert read
        with self._file('raw', 'bagit', 'data/' + file_path) as f:
            yield f

    @contextlib.contextmanager
    def archive_file(self, calc_hash: str, read: bool = True) -> Generator[IO, None, None]:
        assert read
        with self._file('archive', self._archive_ext, '%s.%s' % (calc_hash, self._archive_ext)) as f:
            yield f


    def repack(self) -> None:
        """
        Replaces the existing public/restricted data file pairs with new ones, based
        on current restricted information in the metadata. Should be used after updating
        the restrictions on calculations. This is potentially a long running operation.
        """
        pass
