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
This file storage abstraction uses an *object storage*-like metaphor to store
objects on the file system. Objects are organized in *buckets* and object ids
are basically paths. All major file system operations for dealing with
uploaded files, archive, files, raw files, etc. should be part of this module to
allow later introduction of real object storage systems.

.. note:: This module still uses ``os.path``. As long as the whole nomad runs on a
    POSIX (or Windows) os everything should be fine. This means respective paths in the
    dbs, and indices. In the future, this should be replaced with abstract path representations
    ala ``PathLib``.

.. autoclass:: File
    :members:
.. autoclass:: ZippedFile
    :members:
.. autoclass:: ObjectFile
    :members:
.. autoclass:: UploadFile
    :members:
.. autoclass:: ArchiveFile
    :members:
.. autoclass:: DataContainer
    :members:
.. autoclass:: BaggedDataContainer
    :members:
.. autoclass:: ZippedDataContainer
    :members:
"""
from abc import ABC
from typing import List, Generator, IO, TextIO, cast, Dict, Any
import os
import os.path
from zipfile import ZipFile, BadZipFile, is_zipfile
import shutil
from contextlib import contextmanager
import gzip
import io
import bagit
import json

from nomad import config, utils


class File:
    """
    Base class for handling a file. Allows to open (read, write) and delete files.

    Arguments:
        os_path: The path to the file in the os filesystem.

    Attributes:
        logger: A structured logger with bucket and object information.
        path: The abstract path of the file.
    """
    def __init__(self, os_path: str = None) -> None:
        self.os_path = os_path

        self.logger = self.bind_logger(utils.get_logger(__name__))

    def bind_logger(self, logger):
        """ Adds context information to the given logger and returns it. """
        return logger.bind(path=self.os_path)

    @contextmanager
    def open(self, *args, **kwargs) -> Generator[IO, None, None]:
        """ Opens the object with he given mode, etc. """
        self.logger.debug('open file')
        try:
            with open(self.os_path, *args, **kwargs) as f:
                yield f
        except FileNotFoundError:
            raise KeyError()

    def delete(self) -> None:
        """ Deletes the file. """
        try:
            os.remove(self.os_path)
            self.logger.debug('file deleted')
        except FileNotFoundError:
            raise KeyError()

    def exists(self) -> bool:
        """ Returns true if object exists. """
        return os.path.exists(self.os_path)

    @property
    def size(self) -> int:
        """ Returns the os determined file size. """
        return os.stat(self.os_path).st_size

    @property
    def path(self) -> str:
        return self.os_path


class ZippedFile(File):
    """ A file contained in a .zip archive. """
    def __init__(self, zip_os_path: str, filename: str) -> None:
        self.filename = filename
        super().__init__(zip_os_path)

    def bind_logger(self, logger):
        return super().bind_logger(logger).bind(filename=self.filename)

    @contextmanager
    def open(self, *args, **kwargs) -> Generator[IO, None, None]:
        self.logger.debug('open file')
        try:
            with ZipFile(self.os_path) as zip_file:
                yield zip_file.open(self.filename, *args, **kwargs)
        except FileNotFoundError:
            raise KeyError()
        except KeyError as e:
            raise e
        except Exception as e:
            msg = 'Could not read upload.'
            self.logger.error(msg, exc_info=e)
            raise FileError(msg, e)

    def delete(self) -> None:
        assert False, "A file in a zip archive cannot be deleted."

    @property
    def size(self) -> int:
        with ZipFile(self.os_path) as zip_file:
            return zip_file.getinfo(self.filename).file_size

    @property
    def path(self) -> str:
        return os.path.join(
            os.path.dirname(self.os_path),
            os.path.basename(self.os_path),
            self.filename)


class Objects:
    @classmethod
    def _os_path(cls, bucket: str, name: str, ext: str = None) -> str:
        if ext is not None and ext != '':
            file_name = ".".join([name, ext])
        elif name is None or name == '':
            file_name = ''
        else:
            file_name = name

        path_segments = file_name.split('/')
        path = os.path.join(*([config.fs.objects, bucket] + path_segments))
        directory = os.path.dirname(path)

        if not os.path.isdir(directory):
            os.makedirs(directory)

        return os.path.abspath(path)

    @classmethod
    def delete_all(cls, bucket: str, prefix: str = ''):
        """ Delete all files with given prefix, prefix must denote a directory. """
        try:
            shutil.rmtree(cls._os_path(bucket, prefix, ext=None))
        except FileNotFoundError:
            pass


class ObjectFile(File):
    """
    Base class for file objects. Allows to open (read, write) and delete objects.
    File objects filesystem location is govern by its bucket, object_id, and ext.
    This object store location can be overridden with a local_path.

    Arguments:
        bucket (str): The 'bucket' for this object.
        object_id (str): The object_id for this object. Might contain `/` to structure
            the bucket further. Will be mapped to directories in the filesystem.
        ext (str): Optional extension for the object file in the filesystem.

    Attributes:
        logger: A structured logger with bucket and object information.
        has_local_path: True, if this object is stored somewhere else in the fs.
    """
    def __init__(self, bucket: str, object_id: str, ext: str = None, local_path: str = None) -> None:
        self.bucket = bucket
        self.object_id = object_id
        self.ext = ext

        self.has_local_path = local_path is not None
        path = Objects._os_path(self.bucket, self.object_id, self.ext)
        path = local_path if self.has_local_path else path

        super().__init__(path)

    def bind_logger(self, logger):
        """ Adds context information to the given logger and returns it. """
        return super().bind_logger(logger).bind(bucket=self.bucket, object=self.object_id)

    def delete(self) -> None:
        """ Deletes the file, if it has not a localpath. Localpath files are never deleted.  """
        # Do not delete local files, no matter what
        if not self.has_local_path:
            super().delete()


class FileError(Exception):
    def __init__(self, msg, cause):
        super().__init__(msg, cause)


class UploadFile(ObjectFile):
    """
    Instances of ``UploadFile`` represent an uploaded file in the *'object storage'*.

    Currently only user ``.zip`` files are supported.

    Uploads can be extracted to tmp storage (open/close), the list of files in
    the upload is provided, and files can be opened for read. Extracting uploads
    is optional, all functions in this module are also available without extracting.
    Extracts are automatically bagged with *bagit*.

    This class is a context manager, that extracts the file when using a ``with``
    statement with instances of this class.

    UploadFiles are stored in their own *bucket*. But, storage can be overridden
    by providing a ``local_path``. This is useful when the file is already stored
    in nomad's distributed file system, e.g. for bulk processing of already uploaded
    files.

    Uploads can be persistet as :class:`ZippedDataContainers` for permanent repository
    raw data storage.

    Arguments:
        upload_id: The upload of this uploaded file.
        local_path: Optional override for the path used to store/access the uploaded file.

    Attributes:
        is_extracted: True if the upload is extracted.
        upload_extract_dir: The path of the tmp directory with the extracted contents.
        filelist: A list of filenames relative to the .zipped upload root.
    """

    formats = ['zip']
    """ A human readable list of supported file formats. """

    def __init__(self, upload_id: str, local_path: str = None) -> None:
        super().__init__(
            bucket=config.files.uploads_bucket,
            object_id=upload_id,
            ext='zip',
            local_path=local_path)

        self._extract_dir: str = os.path.join(config.fs.tmp, 'uploads_extracted', upload_id)
        self._bagged_container: DataContainer = None
        if os.path.isdir(self._extract_dir):
            self._bagged_container = BaggedDataContainer(self._extract_dir)

    def bind_logger(self, logger):
        return super().bind_logger(logger).bind(upload_id=self.object_id)

    # There is not good way to capsule decorators in a class:
    # https://medium.com/@vadimpushtaev/decorator-inside-python-class-1e74d23107f6
    class Decorators:
        @classmethod
        def handle_errors(cls, decorated):
            def wrapper(self, *args, **kwargs):
                try:
                    return decorated(self, *args, **kwargs)
                except Exception as e:
                    msg = 'Could not %s upload.' % decorated.__name__
                    self.logger.error(msg, upload_id=self.object_id, exc_info=e)
                    raise FileError(msg, e)
            return wrapper

    @contextmanager
    def _zip(self):
        assert self.exists(), "Can only access uploaded file if it exists."
        zip_file = None
        try:
            zip_file = ZipFile(self.os_path)
            yield zip_file
        except BadZipFile as e:
            raise FileError('Upload is not a zip file', e)
        finally:
            if zip_file is not None:
                zip_file.close()

    @property
    def filelist(self) -> List[str]:
        if self.is_extracted:
            return self._bagged_container.manifest
        else:
            with self._zip() as zip_file:
                return [
                    zip_info.filename for zip_info in zip_file.filelist
                    if not zip_info.filename.endswith('/')]

    @property
    def is_extracted(self) -> bool:
        return self._bagged_container is not None

    @Decorators.handle_errors
    def upload_hash(self) -> str:
        assert self.is_extracted
        return self._bagged_container.hash

    @Decorators.handle_errors
    def extract(self) -> None:
        """
        'Opens' the upload. This means the upload files get extracted and bagged to tmp.

        Raises:
            UploadFileError: If some IO went wrong.
            KeyError: If the upload does not exist.
        """
        os.makedirs(os.path.join(config.fs.tmp, 'uploads_extracted'), exist_ok=True)

        with self._zip() as zip_file:
            zip_file.extractall(self._extract_dir)

        self.logger.debug('extracted uploaded file')

        self._bagged_container = BaggedDataContainer.create(self._extract_dir)
        self.logger.debug('bagged uploaded file')

    def persist(self, object_id: str = None):
        """
        Persists the extracted and bagged upload to the repository raw data bucket.
        """
        assert self.is_extracted
        if object_id is None:
            object_id = self.upload_hash()

        return ZippedDataContainer.create(
            self._extract_dir, Objects._os_path(config.files.repository_bucket, object_id))

    @Decorators.handle_errors
    def remove_extract(self) -> None:
        """
        Closes the upload. This means the tmp. files are deleted.

        Raises:
            UploadFileError: If some IO went wrong.
            KeyError: If the upload does not exist.
        """
        try:
            shutil.rmtree(self._extract_dir)
        except FileNotFoundError:
            raise KeyError()

        self.logger.debug('removed uploaded file extract')

    def __enter__(self):
        self.extract()
        return self

    def __exit__(self, exc_type, exc, exc_tb):
        self.remove_extract()

    def get_file(self, filename: str) -> File:
        """
        Returns a :class:`File` instance as a handle to the file with the given name.
        Only works on extracted uploads. The given filename must be one of the
        name in ``self.filelist``.
        """
        assert self.is_extracted
        return self._bagged_container.get_file(filename)

    @property
    def is_valid(self):
        return is_zipfile(self.os_path)

    def get_siblings(self, filename: str) -> Generator[str, None, None]:
        """
        Returns the names of all files that share the same prefix (object id),
        respectively are part of the same directory (incl. files in sub directories).
        In nomad terms, the aux files the this file. Returned siblings are relative
        to this files directory.
        """
        dirname = os.path.dirname(filename)
        dirname_len = len(dirname) + 1
        for other in self.filelist:
            if other.startswith(dirname) and other != filename:
                yield other[dirname_len:]

    def get_sibling_file(self, filename: str, sibling: str) -> File:
        sibling_name = os.path.join(os.path.dirname(filename), sibling)
        return self.get_file(sibling_name)


class RepositoryFile(ObjectFile):
    """
    Represents a repository file. A repository file is a persistet bagged upload, incl.
    the upload metadata. It is used to serve raw data.
    """
    def __init__(self, upload_hash: str) -> None:
        super().__init__(
            bucket=config.files.repository_bucket,
            object_id=upload_hash,
            ext='zip')

        self._zipped_container = ZippedDataContainer(self.os_path)

    def get_file(self, path: str) -> ZippedFile:
        return self._zipped_container.get_file(path)


class ArchiveFile(ObjectFile):
    """
    Represents the archive file for an individual calculation. Allows to write the
    archive, read the archive, delete the archive.

    Archive files are stored in their own *bucket*.
    """
    def __init__(self, archive_id: str) -> None:
        super().__init__(
            bucket=config.files.archive_bucket,
            object_id=archive_id,
            ext='json.gz' if config.files.compress_archive else 'json')

    def bind_logger(self, logger):
        upload_hash, calc_hash = self.object_id.split('/')
        return super().bind_logger(logger).bind(
            archive_id=self.object_id, upload_hash=upload_hash, calc_hash=calc_hash)

    @contextmanager
    def write_archive_json(self) -> Generator[TextIO, None, None]:
        """ Context manager that yields a file-like to write the archive json. """
        with self.open('wb') as binary_out:
            if config.files.compress_archive:
                gzip_wrapper = cast(TextIO, gzip.open(binary_out, 'wt'))
                out = gzip_wrapper
            else:
                text_wrapper = io.TextIOWrapper(binary_out, encoding='utf-8')
                out = text_wrapper

            try:
                yield out
            finally:
                out.flush()
                out.close()

        self.logger.debug('archive file written')

    @contextmanager
    def read_archive_json(self) -> Generator[TextIO, None, None]:
        """ Context manager that yields a file-like to read the archive json. """
        with self.open(mode='rb') as binary_in:
            try:
                if config.files.compress_archive:
                    gzip_wrapper = cast(TextIO, gzip.open(binary_in, 'rt'))
                    in_file = gzip_wrapper
                else:
                    text_wrapper = io.TextIOWrapper(binary_in, encoding='utf-8')
                    in_file = text_wrapper
            except FileNotFoundError:
                raise KeyError()

            try:
                yield in_file
            finally:
                in_file.close()

        self.logger.debug('archive file read')

    @staticmethod
    def delete_archives(upload_hash: str):
        """ Delete all archives of one upload with the given hash. """
        bucket = config.files.archive_bucket
        Objects.delete_all(bucket, upload_hash)

        utils.get_logger(__name__, bucket=bucket, upload_hash=upload_hash) \
            .debug('archive files deleted')


class ArchiveLogFile(ObjectFile):
    """
    Represents a log file that was created for processing a single calculation to create
    an archive.
    Logfiles are stored within the *archive_bucket* alongside the archive files.
    """
    def __init__(self, archive_id: str) -> None:
        super().__init__(
            bucket=config.files.archive_bucket,
            object_id=archive_id,
            ext='log')


class DataContainer(ABC):
    """
    An abstract baseclass for a *data container*. A data container is a persistent
    bundle of related files, like the calculation raw data of a user upload.

    A container has a *manifest* and arbitrary *metadata*.
    """
    @property
    def manifest(self) -> List[str]:
        """
        A readonly list of paths to files within the container relative to the containers
        payload directory.
        """
        pass

    @property
    def metadata(self) -> Dict[str, Any]:
        """
        The modifiable metadata of this manifest. On the top-level its a string keyed
        dictionary. The values can be arbitrary, but have to be JSON-serializable.
        Modifications have to be saved (:func:`save_metadata`).
        """
        pass

    def save_metadata(self) -> None:
        """ Persists metadata changes. """
        pass

    def get_file(self, manifest_path: str) -> File:
        """
        Returns a file-like for the given manifest path.
        """
        pass

    @property
    def hash(self) -> str:
        return self.metadata['Nomad-Hash']


class BaggedDataContainer(DataContainer):
    """
    A *data container* based on *bagit*. Once created no more files can be added.
    """
    def __init__(self, path: str) -> None:
        self.path = path
        self.bag = bagit.Bag(path)
        self._metadata = None
        self.payload_directory = os.path.join(path, 'data')

    @staticmethod
    def create(path: str) -> 'BaggedDataContainer':
        """
        Makes a bag from the given directory and returns the respective BaggedDataContainer
        instance.
        """
        bag = bagit.make_bag(path, checksums=['sha512'])
        hashes = [
            value['sha512'] for key, value in bag.entries.items()
            if key.startswith('data/')
        ]
        bag.info['Nomad-Hash'] = utils.hash(''.join(hashes))
        bag.save()
        return BaggedDataContainer(path)

    @property
    def metadata(self):
        if self._metadata is None:
            self._metadata = BaggedDataContainer._load_bagit_metadata(self.bag.info)
        return self._metadata

    @staticmethod
    def _load_bagit_metadata(info):
        metadata = info
        for key, value in metadata.items():
            if key not in bagit.STANDARD_BAG_INFO_HEADERS:
                try:
                    metadata[key] = json.loads(value)
                except Exception:
                    pass
        return metadata

    def save_metadata(self):
        metadata = self.bag.info
        for key, value in metadata.items():
            if key not in bagit.STANDARD_BAG_INFO_HEADERS and not isinstance(value, str):
                metadata[key] = json.dumps(value)
        self.bag.save()

    @property
    def manifest(self):
        return [path[5:] for path in self.bag.entries.keys() if path.startswith('data/')]

    def get_file(self, path):
        return File(os.path.join(self.payload_directory, path))


class ZippedDataContainer(File, DataContainer):
    """
    A *bagit*-based data container that has been zipped. Its metadata cannot be changed
    anymore.
    """
    def __init__(self, os_path: str) -> None:
        super(ZippedDataContainer, self).__init__(os_path)
        self._metadata = None

    @staticmethod
    def create(path: str, target: str = None) -> 'ZippedDataContainer':
        if not target:
            target = path

        assert os.path.isdir(path)
        archive_file = shutil.make_archive(target, 'zip', path)
        return ZippedDataContainer(archive_file)

    @contextmanager
    def _zip(self):
        assert self.exists(), "Can only access uploaded file if it exists."
        zip_file = None
        try:
            zip_file = ZipFile(self.os_path)
            yield zip_file
        except BadZipFile as e:
            raise FileError('Upload is not a zip file', e)
        finally:
            if zip_file is not None:
                zip_file.close()

    @property
    def manifest(self):
        with self._zip() as zip_file:
            return [
                zip_info.filename[5:] for zip_info in zip_file.filelist
                if not zip_info.filename.endswith('/') and zip_info.filename.startswith('data/')]

    @property
    def metadata(self):
        if self._metadata is None:
            self._metadata = self._load_metadata()
        return self._metadata

    def _load_metadata(self):
        with ZippedFile(self.os_path, 'bag-info.txt').open('r') as metadata_file:
            metadata_contents = metadata_file.read()

        metadata_file = io.StringIO(metadata_contents.decode("utf-8"))
        tags = {}
        for name, value in bagit._parse_tags(metadata_file):
            if name not in tags:
                tags[name] = value
                continue

            if not isinstance(tags[name], list):
                tags[name] = [tags[name], value]
            else:
                tags[name].append(value)

        print(tags)
        return BaggedDataContainer._load_bagit_metadata(tags)

    def get_file(self, path):
        return ZippedFile(self.path, 'data/' + path)
