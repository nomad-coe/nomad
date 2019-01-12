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
from zipfile import ZipFile, BadZipFile, is_zipfile
from bagit import make_bag
import hashlib
import io

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

    @property
    def size(self) -> int:
        """ Returns the os determined file size. """
        return os.stat(self.os_path).st_size

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


class UploadFiles(DirectoryObject, metaclass=ABCMeta):

    _archive_ext = 'json'

    def __init__(
            self, bucket: str, upload_id: str,
            is_authorized: Callable[[], bool] = lambda: False,
            create: bool = False) -> None:
        self.logger = utils.get_logger(__name__, upload_id=upload_id)

        super().__init__(bucket, upload_id, create=create, prefix=True)

        if not create and not self.exists():
            raise KeyError()

        self.upload_id = upload_id
        self._is_authorized = is_authorized

    @staticmethod
    def get(upload_id: str, *args, **kwargs) -> 'UploadFiles':
        if DirectoryObject(config.files.staging_bucket, upload_id, prefix=True).exists():
            return StagingUploadFiles(upload_id, *args, **kwargs)
        elif DirectoryObject(config.files.public_bucket, upload_id, prefix=True).exists():
            return PublicUploadFiles(upload_id, *args, **kwargs)
        else:
            return None

    @property
    def metadata(self) -> Metadata:
        """ The calc metadata for this upload. """
        raise NotImplementedError

    def raw_file(self, file_path: str, *args, **kwargs) -> IO:
        """
        Opens a raw file and returns a file-like objects. Additional args, kwargs are
        delegated to the respective `open` call.
        Arguments:
            file_path: The path to the file relative to the upload.
        Raises:
            KeyError: If the file does not exist.
            Restricted: If the file is restricted and upload access evaluated to False.
        """
        raise NotImplementedError()

    def raw_file_manifest(self, path_prefix: str = None) -> Generator[str, None, None]:
        """
        Returns the path for all raw files in the archive (with a given prefix).
        Arguments:
            path_prefix: An optional prefix; only returns those files that have the prefix.
        Returns:
            An iterable over all (matching) raw files.
        """
        raise NotImplementedError()

    def archive_file(self, calc_hash: str, *args, **kwargs) -> IO:
        """
        Opens a archive file and returns a file-like objects. Additional args, kwargs are
        delegated to the respective `open` call.
        Arguments:
            calc_hash: The hash identifying the calculation.
        Raises:
            KeyError: If the calc does not exist.
            Restricted: If the file is restricted and upload access evaluated to False.
        """
        raise NotImplementedError()

    def archive_log_file(self, calc_hash: str, *args, **kwargs) -> IO:
        """
        Opens a archive log file and returns a file-like objects. Additional args, kwargs are
        delegated to the respective `open` call.
        Arguments:
            calc_hash: The hash identifying the calculation.
        Raises:
            KeyError: If the calc does not exist.
            Restricted: If the file is restricted and upload access evaluated to False.
        """
        raise NotImplementedError()


class StagingUploadFiles(UploadFiles):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(config.files.staging_bucket, *args, **kwargs)

        self._raw_dir = self.join_dir('raw')
        self._archive_dir = self.join_dir('archive')
        self._frozen_file = self.join_file('.frozen')

        metadata_dir = self.join_dir('metadata')
        self._metadata = StagingMetadata(metadata_dir)

        self._size = 0

    @property
    def size(self) -> int:
        return self._size

    @property
    def metadata(self) -> Metadata:
        return self._metadata

    def _file(self, path_object: PathObject, *args, **kwargs) -> IO:
        try:
            return open(path_object.os_path, *args, **kwargs)
        except FileNotFoundError:
            raise KeyError()

    def raw_file(self, file_path: str, *args, **kwargs) -> IO:
        if not self._is_authorized():
            raise Restricted
        return self._file(self.raw_file_object(file_path), *args, **kwargs)

    def raw_file_object(self, file_path: str) -> PathObject:
        return self._raw_dir.join_file(file_path)

    def archive_file(self, calc_hash: str, *args, **kwargs) -> IO:
        if not self._is_authorized():
            raise Restricted
        return self._file(self.archive_file_object(calc_hash), *args, **kwargs)

    def archive_log_file(self, calc_hash: str, *args, **kwargs) -> IO:
        if not self._is_authorized():
            raise Restricted
        return self._file(self.archive_log_file_object(calc_hash), *args, **kwargs)

    def archive_file_object(self, calc_hash: str) -> PathObject:
        return self._archive_dir.join_file('%s.%s' % (calc_hash, self._archive_ext))

    def archive_log_file_object(self, calc_hash: str) -> PathObject:
        return self._archive_dir.join_file('%s.log' % calc_hash)

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
        self._size += os.stat(path).st_size
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

        # create tmp dirs for restricted and public raw data
        restricted_dir = self.join_dir('.restricted', create=False)
        public_dir = self.join_dir('.public', create=True)

        # copy raw -> .restricted
        shutil.copytree(self._raw_dir.os_path, restricted_dir.os_path)

        # move public data .restricted -> .public
        for calc in self.metadata:
            if not calc.get('restricted', True):
                mainfile: str = calc['mainfile']
                assert mainfile is not None
                for filepath in self.calc_files(mainfile):
                    os.rename(
                        restricted_dir.join_file(filepath).os_path,
                        public_dir.join_file(filepath).os_path)

        # create bags
        make_bag(restricted_dir.os_path, bag_info=bagit_metadata, checksums=['sha512'])
        make_bag(public_dir.os_path, bag_info=bagit_metadata, checksums=['sha512'])

        # zip bags
        def zip_dir(zip_filepath, path):
            root_len = len(path)
            with ZipFile(zip_filepath, 'w') as zf:
                for root, _, files in os.walk(path):
                    for file in files:
                        filepath = os.path.join(root, file)
                        zf.write(filepath, filepath[root_len:])

        packed_dir = self.join_dir('.packed', create=True)

        zip_dir(packed_dir.join_file('raw-restricted.bagit.zip').os_path, restricted_dir.os_path)
        zip_dir(packed_dir.join_file('raw-public.bagit.zip').os_path, public_dir.os_path)

        # zip archives
        def create_zipfile(prefix: str) -> ZipFile:
            file = packed_dir.join_file('archive-%s.%s.zip' % (prefix, self._archive_ext))
            return ZipFile(file.os_path, mode='w')

        archive_public_zip = create_zipfile('public')
        archive_restricted_zip = create_zipfile('restricted')

        for calc in self.metadata:
            archive_zip = archive_restricted_zip if calc.get('restricted', False) else archive_public_zip

            archive_filename = '%s.%s' % (calc['hash'], self._archive_ext)
            archive_zip.write(self._archive_dir.join_file(archive_filename).os_path, archive_filename)

            archive_log_filename = '%s.%s' % (calc['hash'], 'log')
            log_file = self._archive_dir.join_file(archive_log_filename)
            if log_file.exists():
                archive_zip.write(log_file.os_path, archive_log_filename)

        archive_restricted_zip.close()
        archive_public_zip.close()

        # pack metadata
        with PublicMetadata(packed_dir.os_path) as packed_metadata:
            for calc in self.metadata:
                packed_metadata.insert(calc)

        # move to public bucket
        target_dir = DirectoryObject(config.files.public_bucket, self.upload_id, create=False, prefix=True)
        assert not target_dir.exists()
        shutil.move(packed_dir.os_path, target_dir.os_path)

    def raw_file_manifest(self, path_prefix: str = None) -> Generator[str, None, None]:
        upload_prefix_len = len(self._raw_dir.os_path) + 1
        for root, _, files in os.walk(self._raw_dir.os_path):
            for file in files:
                path = os.path.join(root, file)[upload_prefix_len:]
                if path_prefix is None or path.startswith(path_prefix):
                    yield path

    def calc_files(self, mainfile: str, with_mainfile: bool = True) -> Iterable[str]:
        """
        Returns all the auxfiles and mainfile for a given mainfile. This implements
        nomad's logic about what is part of a calculation and what not.
        Arguments:
            mainfile: The mainfile relative to upload
            with_mainfile: Do include the mainfile, default is True
        """
        mainfile_object = self._raw_dir.join_file(mainfile)
        if not mainfile_object.exists():
            raise KeyError()

        mainfile = os.path.basename(mainfile)
        calc_dir = os.path.dirname(mainfile_object.os_path)
        calc_relative_dir = calc_dir[len(self._raw_dir.os_path) + 1:]
        return sorted(
            os.path.join(calc_relative_dir, path) for path in os.listdir(calc_dir)
            if os.path.isfile(os.path.join(calc_dir, path)) and (with_mainfile or path != mainfile))

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
            with open(self._raw_dir.join_file(filepath).os_path, 'rb') as f:
                for data in iter(lambda: f.read(65536), b''):
                    hash.update(data)

        return utils.websave_hash(hash.digest(), utils.default_hash_len)


class ArchiveBasedStagingUploadFiles(StagingUploadFiles):
    """
    :class:`StagingUploadFiles` based on a single uploaded archive file (.zip)

    Arguments:
        local_path: Optional override for the path used to store/access the uploaded file.
    """

    formats = ['zip']
    """ A human readable list of supported file formats. """

    def __init__(self, upload_id: str, local_path: str = None, *args, **kwargs) -> None:
        super().__init__(upload_id, *args, **kwargs)
        self._local_path = local_path
        self._upload_file = self.join_file('.upload.zip')

    @property
    def upload_file_os_path(self):
        if self._local_path is not None:
            return self._local_path
        else:
            return self._upload_file.os_path

    @property
    def is_valid(self) -> bool:
        if not os.path.exists(self.upload_file_os_path):
            return False
        elif not os.path.isfile(self.upload_file_os_path):
            return False
        else:
            return is_zipfile(self.upload_file_os_path)

    def extract(self) -> None:
        assert next(self.raw_file_manifest(), None) is None, 'can only extract once'
        super().add_rawfiles(self.upload_file_os_path)

    def add_rawfiles(self, path: str, move: bool = False, prefix: str = None) -> None:
        assert False, 'do not add_rawfiles to a %s' % self.__class__.__name__


class PublicUploadFiles(UploadFiles):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(config.files.public_bucket, *args, **kwargs)

        self._metadata = PublicMetadata(self.os_path)

    @property
    def metadata(self) -> Metadata:
        return self._metadata

    def _file(self, prefix: str, ext: str, path: str, *args, **kwargs) -> IO:
        mode = kwargs.get('mode') if len(args) == 0 else args[0]
        if 'mode' in kwargs:
            del(kwargs['mode'])
        mode = mode if mode else 'rb'

        for access in ['public', 'restricted']:
            try:
                zip_file = self.join_file('%s-%s.%s.zip' % (prefix, access, ext))
                with ZipFile(zip_file.os_path) as zf:
                    f = zf.open(path, 'r', **kwargs)
                    if access == 'restricted' and not self._is_authorized():
                        raise Restricted
                    if 't' in mode:
                        return io.TextIOWrapper(f)
                    else:
                        return f
            except FileNotFoundError:
                pass
            except KeyError:
                pass

        raise KeyError()

    def raw_file(self, file_path: str, *args, **kwargs) -> IO:
        return self._file('raw', 'bagit', 'data/' + file_path, *args, *kwargs)

    def raw_file_manifest(self, path_prefix: str = None) -> Generator[str, None, None]:
        for access in ['public', 'restricted']:
            try:
                zip_file = self.join_file('raw-%s.bagit.zip' % access)
                with ZipFile(zip_file.os_path) as zf:
                    for full_path in zf.namelist():
                        path = full_path[5:]  # remove data/
                        if path_prefix is None or path.startswith(path_prefix):
                            yield path
            except FileNotFoundError:
                pass

    def archive_file(self, calc_hash: str, *args, **kwargs) -> IO:
        return self._file('archive', self._archive_ext, '%s.%s' % (calc_hash, self._archive_ext), *args, **kwargs)

    def archive_log_file(self, calc_hash: str, *args, **kwargs) -> IO:
        return self._file('archive', self._archive_ext, '%s.log' % calc_hash, *args, **kwargs)

    def repack(self) -> None:
        """
        Replaces the existing public/restricted data file pairs with new ones, based
        on current restricted information in the metadata. Should be used after updating
        the restrictions on calculations. This is potentially a long running operation.
        """
        pass
