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

Uploads
-------
.. autoclass:: File
    :members:
.. autoclass:: UploadFile
    :members:
.. autoclass:: ArchiveFile
    :members:

"""
from typing import List, Any, Generator, IO, TextIO, cast
import os
import os.path
from zipfile import ZipFile, BadZipFile, is_zipfile
import shutil
from contextlib import contextmanager
import gzip
import io
import shutil

from nomad import config, utils


class Objects:
    @classmethod
    def _os_path(cls, bucket: str, name: str, ext: str) -> str:
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


class File:
    """
    Base class for file objects. Allows to open (read, write) and delete objects.

    Arguments:
        bucket (str): The 'bucket' for this object.
        object_id (str): The object_id for this object. Might contain `/` to structure
            the bucket further. Will be mapped to directories in the filesystem.
        ext (str): Optional extension for the object file in the filesystem.

    Attributes:
        logger: A structured logger with bucket and object information.
    """
    def __init__(self, bucket: str, object_id: str, ext: str = None) -> None:
        self.bucket = bucket
        self.object_id = object_id
        self.ext = ext

        self.logger = self.bind_logger(utils.get_logger(__name__))

    def bind_logger(self, logger):
        """ Adds context information to the given logger and returns it. """
        return logger.bind(bucket=self.bucket, object=self.object_id)

    def open(self, *args, **kwargs) -> IO:
        """ Opens the object with he given mode, etc. """
        self.logger.debug('open file')
        try:
            return open(self.os_path, *args, **kwargs)
        except FileNotFoundError:
            raise KeyError()

    def delete(self) -> None:
        """ Deletes the file with the given object id. """
        try:
            os.remove(self.os_path)
            self.logger.debug('file deleted')
        except FileNotFoundError:
            raise KeyError()

    def exists(self) -> bool:
        """ Returns true if object exists. """
        return os.path.exists(self.os_path)

    @property
    def os_path(self) -> str:
        """ The path of the object in the os filesystem. """
        return Objects._os_path(self.bucket, self.object_id, self.ext)


class FileError(Exception):
    def __init__(self, msg, cause):
        super().__init__(msg, cause)


class UploadFile(File):
    """
    Instances represent an uploaded file in the *object storage*. Class is a conext
    manager and supports the `with` statements.
    In conext the upload will be extracted and contained files can be opened.
    Some functions are only available for extracted uploads.

    Uploads are stored in their own *bucket*.

    Arguments:
        upload_id: The upload of this uploaded file.
        local_path: Optional override for the path used to store/access the upload
            on the server. This can be usefull to create uploads for files that
            were not uploaded but put to the server in another way, e.g. offline
            processing, syncing with other data, etc.

    Attributes:
        upload_extract_dir: The path of the tmp directory with the extracted contents.
        filelist: A list of filenames relative to the .zipped upload root.
    """

    formats = ['zip']
    """ A human readable list of supported file formats. """

    def __init__(self, upload_id: str, local_path: str = None) -> None:
        super().__init__(
            bucket=config.files.uploads_bucket,
            object_id=upload_id,
            ext='zip')

        self.upload_extract_dir: str = os.path.join(config.fs.tmp, 'uploads_extracted', upload_id)
        self.filelist: List[str] = None
        self._local_path = local_path

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

    @property
    def os_path(self):
        return self._local_path if self._local_path is not None else super().os_path

    @Decorators.handle_errors
    def hash(self) -> str:
        """ Calculates the first 28 bytes of a websafe base64 encoded SHA512 of the upload. """
        with self.open('rb') as f:
            return utils.hash(f)

    @Decorators.handle_errors
    def extract(self) -> None:
        """
        'Opens' the upload. This means the upload files get extracted to tmp.

        Raises:
            UploadFileError: If some IO went wrong.
            KeyError: If the upload does not exist.
        """
        os.makedirs(os.path.join(config.fs.tmp, 'uploads_extracted'), exist_ok=True)

        zipFile = None
        try:
            zipFile = ZipFile(self.os_path)
            zipFile.extractall(self.upload_extract_dir)
            self.filelist = [
                zipInfo.filename for zipInfo in zipFile.filelist
                if not zipInfo.filename.endswith('/')]
        except BadZipFile as e:
            raise FileError('Upload is not a zip file', e)
        finally:
            if zipFile is not None:
                zipFile.close()

        self.logger.debug('extracted uploaded file')

    @Decorators.handle_errors
    def remove_extract(self) -> None:
        """
        Closes the upload. This means the tmp. files are deleted.

        Raises:
            UploadFileError: If some IO went wrong.
            KeyError: If the upload does not exist.
        """
        try:
            shutil.rmtree(self.upload_extract_dir)
        except FileNotFoundError:
            raise KeyError()

        self.logger.debug('removed uploaded file extract')

    def __enter__(self):
        self.extract()
        return self

    def __exit__(self, exc_type, exc, exc_tb):
        self.remove_extract()

    @Decorators.handle_errors
    def open_file(self, filename: str, *args, **kwargs) -> IO[Any]:
        """ Opens a file within an open upload and returns a file like. """
        return open(self.get_path(filename), *args, **kwargs)

    def get_path(self, filename: str) -> str:
        """ Returns the tmp directory relative version of a filename. """
        return os.path.join(self.upload_extract_dir, filename)

    def delete(self) -> None:
        """ Deletes the file with the given object id. """
        # Do not delete local files, no matter what
        if self._local_path is None:
            try:
                os.remove(self.os_path)
                self.logger.debug('file deleted')
            except FileNotFoundError:
                raise KeyError()

    @property
    def is_valid(self):
        return is_zipfile(self.os_path)


class ArchiveFile(File):
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
        if config.files.compress_archive:
            binary_out = self.open('wb')
            gzip_wrapper = cast(TextIO, gzip.open(binary_out, 'wt'))
            out = gzip_wrapper
        else:
            binary_out = self.open('wb')
            text_wrapper = io.TextIOWrapper(binary_out, encoding='utf-8')
            out = text_wrapper

        try:
            yield out
        finally:
            out.flush()
            out.close()
            binary_out.close()

        self.logger.debug('archive file written')

    @contextmanager
    def read_archive_json(self) -> Generator[TextIO, None, None]:
        """ Context manager that yields a file-like to read the archive json. """
        try:
            if config.files.compress_archive:
                binary_in = self.open(mode='rb')
                gzip_wrapper = cast(TextIO, gzip.open(binary_in, 'rt'))
                in_file = gzip_wrapper
            else:
                binary_in = self.open(mode='rb')
                text_wrapper = io.TextIOWrapper(binary_in, encoding='utf-8')
                in_file = text_wrapper
        except FileNotFoundError:
            raise KeyError()

        try:
            yield in_file
        finally:
            in_file.close()
            binary_in.close()

        self.logger.debug('archive file read')

    @staticmethod
    def delete_archives(upload_hash: str):
        """ Delete all archives of one upload with the given hash. """
        bucket = config.files.archive_bucket
        Objects.delete_all(bucket, upload_hash)

        utils.get_logger(__name__, bucket=bucket, upload_hash=upload_hash) \
            .debug('archive files deleted')


class ArchiveLogFile(File):
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
