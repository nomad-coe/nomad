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

from typing import IO, Generator, Dict, Any, Iterable
from filelock import Timeout, FileLock
import ujson
import os.path
import os
import shutil
from zipfile import ZipFile, BadZipFile
from bagit import make_bag

from nomad import config, utils


class PathObject:
    def __init__(self, bucket: str, object_id: str, os_path: str = None) -> None:
        if os_path:
            self.os_path = os_path
        else:
            self.os_path = os.path.join(config.fs.objects, bucket, object_id)

    def delete(self) -> None:
        shutil.rmtree(self.os_path)

    def exists(self) -> bool:
        return os.path.exists(self.os_path)

    def __repr__(self) -> str:
        return self.os_path


class FileObject(PathObject):
    def __init__(self, bucket: str, object_id: str, create: bool = False, **kwargs) -> None:
        super().__init__(bucket, object_id, **kwargs)
        dirname = os.path.dirname(self.os_path)
        if create and not os.path.isdir(dirname):
            os.makedirs(dirname)


class DirectoryObject(PathObject):
    def __init__(self, bucket: str, object_id: str, create: bool = False, **kwargs) -> None:
        super().__init__(bucket, object_id, **kwargs)
        self._create = create
        if create and not os.path.isdir(self.os_path):
            os.makedirs(self.os_path)

    def join_dir(self, path, create: bool = None) -> 'DirectoryObject':
        if create is None:
            create = self._create
        return DirectoryObject(None, None, create=create, os_path=os.path.join(self.os_path, path))

    def join_file(self, path) -> 'FileObject':
        return FileObject(None, None, os_path=os.path.join(self.os_path, path))


class MetadataTimeout(Exception):
    pass


class Metadata():
    """
    A contextmanager that wraps around a metadata dictionary. It loads and write
    the metadata to the given path and uses a lock to deal with concurrent access.

    Arguments:
        path: The parent directory for the metadata and lock file.
        lock_timeout: Max timeout before __enter__ raises MetadataTimeout while waiting
            for an available lock on the metadata file. Default is 1s.
        calc_id_key: The key used for ensuring uniqueness on calc metadata. Default is 'hash'.
    """
    def __init__(self, path: str, lock_timeout=1, calc_id_key='hash') -> None:
        self._db_file = os.path.join(path, 'metadata.json')
        self._lock_file = os.path.join(path, 'metadata.json.lock')
        self._lock: FileLock = FileLock(self._lock_file, timeout=lock_timeout)
        self._modified = False
        self._calc_id_key = calc_id_key
        self.data: Dict[Any, dict] = None

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
        """ Insert a calc, using hash as key. """
        assert self.data is not None, "Metadata is not open."

        id = calc[self._calc_id_key]
        assert id not in self.data
        self.data[id] = calc
        self._modified = True

    def update(self, calc: dict) -> None:
        """ Updating a calc, using hash as key and running dict update with the given data. """
        assert self.data is not None, "Metadata is not open."

        id = calc[self._calc_id_key]
        if id not in self.data:
            raise KeyError()
        self.data[id].update(calc)
        self._modified = True

    def get(self, calc_id: Any) -> dict:
        """ Retrive the calc metadata for a given calc. """
        assert self.data is not None, "Metadata is not open."

        return self.data[calc_id]


class UploadFiles():
    def __init__(self, upload_id: str) -> None:
        self.logger = utils.get_logger(__name__, upload_id=upload_id)

    def get_metadata(self, calc_hash: str) -> dict:
        """
        Returns: the metadata for the given calc.
        Arguments:
            calc_hash: The hash identifying the calculation.
        Raises:
            KeyError: If the calc does not exist.
        """
        raise NotImplementedError()

    def raw_file(self, file_path: str, *args, **kwargs) -> Generator[IO, None, None]:
        """
        Opens a raw file and returns a file-like objects. Additional arguments are
        based to the respective low-level open function.
        Arguments:
            file_path: The path to the file relative to the upload.
        Raises:
            KeyError: If the file does not exist.
        """
        raise NotImplementedError()

    def archive_file(self, calc_hash: str, extension: str, *args, **kwargs) -> Generator[IO, None, None]:
        """
        Opens a archive file and returns a file-like objects. Additional arguments are
        based to the respective low-level open function.
        Arguments:
            calc_hash: The hash identifying the calculation.
        Raises:
            KeyError: If the calc does not exist.
        """
        raise NotImplementedError()

    def all_metadata(self) -> Iterable[Dict[Any, dict]]:
        """ Returns: An iterable with the metadata of all calcs """
        raise NotImplementedError()

    def update_metadata(self, calc_hash: str, updates: dict) -> None:
        """
        Allows to update calculation metadata.
        Raises:
            KeyError: If the calc does not exist.
        """
        raise NotImplementedError()


class StagingUploadFiles(UploadFiles):
    def __init__(self, upload_id: str, create: bool = False) -> None:
        super().__init__(upload_id=upload_id)
        self._upload_dir = DirectoryObject(config.files.staging_bucket, upload_id, create=create)
        if not create and not self._upload_dir.exists():
            raise KeyError()
        self._raw_dir = self._upload_dir.join_dir('raw')
        self._archive_dir = self._upload_dir.join_dir('archive')
        self._metadata_dir = self._upload_dir.join_dir('metadata')
        self._frozen_file = self._upload_dir.join_file('.frozen')
        self._restricted_dir = self._upload_dir.join_dir('.restricted', create=False)
        self._public_dir = self._upload_dir.join_dir('.public', create=False)

    def raw_file(self, file_path: str, *args, **kwargs) -> Generator[IO, None, None]:
        path = os.path.join(self._raw_dir.os_path, file_path)
        try:
            with open(path, *args, **kwargs) as f:
                yield f
        except FileNotFoundError:
            raise KeyError()

    def archive_file(self, calc_hash: str, extension: str, *args, **kwargs) -> Generator[IO, None, None]:
        path = os.path.join(self._archive_dir.os_path, '%s.%s' % (calc_hash, extension))
        if not os.path.exists(path):
            raise KeyError()
        try:
            with open(path, *args, **kwargs) as f:
                yield f
        except FileNotFoundError:
            raise KeyError()

    def all_metadata(self) -> Iterable[dict]:
        pass

    def get_metadata(self, calc_hash: str) -> dict:
        pass

    def update_metadata(self, calc_hash: str, updates: dict) -> None:
        pass

    def add_rawfiles(self, path: str) -> None:
        """
        Add rawfiles to the upload. The given file will be moved, or extracted.
        Arguments:
            path: Path to a directory, file, or zip file. Zip files will be extracted.
        """
        assert not self.is_frozen
        assert os.path.exists(path)
        ext = os.path.splitext(path)[1]
        if ext == 'zip':
            try:
                with ZipFile(path) as zf:
                    zf.extractall(self._raw_dir.os_path)
                return
            except BadZipFile:
                pass

        shutil.move(path, self._raw_dir.os_path)

    def insert_metadata(self, calc_hash: str, calc: dict) -> None:
        """ Allows to add calculation metadata. """
        metadata_path = self._metadata_dir.join_file(calc_hash)
        with open(metadata_path.os_path, 'wt') as f:
            ujson.dump(calc, f, ensure_ascii=False)

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
        for calc in self.all_metadata():
            if not calc.get('restricted', False):
                mainfile: str = calc['mainfile']
                dirname = os.path.dirname(mainfile)
                target_dir = self._public_dir.join_dir(dirname, create=False)
                if target_dir.exists():
                    # TODO this is an indicator that one calculation was uploaded in a
                    # subdirectory to another calculation. The packing gets more complex
                    # if we want to support this.
                    self.logger.error('nested calculation raw data', calc_hash=calc['hash'])
                    continue
                shutil.move(os.path.join(self._restricted_dir.os_path, os.path.dirname(mainfile)), target_dir.os_path)

        # create bags
        make_bag(self._restricted_dir, bag_info=bagit_metadata, checksums=['sha512'])
        make_bag(self._public_dir, bag_info=bagit_metadata, checksums=['sha512'])

        # zip bags
        def zip_dir(zip_filepath, path):
            with ZipFile(zip_filepath, 'w') as zf:
                for root, _, files in os.walk(path):
                    for file in files:
                        zf.write(os.path.join(root, file))

        packed_dir = self._upload_dir.join_dir('.packed', create=True)

        zip_dir(packed_dir.join_file('raw-restricted.bagit.zip').os_path, self._restricted_dir)
        zip_dir(packed_dir.join_file('raw-public.bagit.zip').os_path, self._restricted_dir)

        # zip archives
        # archive_public_zip = ZipFile(packed_dir.join_file('archive-public.json.zip'))
        # archive_restricted_zip = ZipFile(packed_dir.join_file('archive-restricted.json.zip'))

        # try:
        #     for calc in self.all_metadata():
        #         if not calc.get('restricted', False):
        #             # public
        #             pass
        #         else:
        #             # restricted
        #             pass
        # finally:
        #     archive_public_zip.close()
        #     archive_restricted_zip.close()

        # move metadata

    def all_files(self) -> Generator[str, None, None]:
        """ Returns: A generator of all file paths of all raw files. """
        pass

    def calc_hash(self, mainfile: str) -> str:
        """
        Calculates a hash for the given calc. It is only available if upload *is_bag*.
        Arguments:
            mainfile: The mainfile path relative to the upload that identifies the calc in the folder structure.
        Returns:
            The calc hash
        Raises:
            KeyError: If the mainfile does not exist.
        """

    def upload_hash(self) -> str:
        """ Returns: A hash for the whole upload. It is only available if upload *is_bag*. """
        pass


class PublicUploadFiles(UploadFiles):
    def __init__(self, upload_hash: str) -> None:
        pass

    def repack(self) -> None:
        """
        Replaces the existing public/restricted data file pairs with new ones, based
        on current restricted information in the metadata. Should be used after updating
        the restrictions on calculations. This is potentially a long running operation.
        """
        pass
