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
This file storage abstraction currently uses the object storage API
http://minio.io to manage and organize files. Object storage
organizes files in *buckets/ids*, with small amounts of *buckets* and virtually
unlimited numbers of *ids*. *Ids* can contain delimiters like `/` to mimic
filesystem structure. There is a 1024 utf-8 character limit on *id* length.

The file storage is organized in multiple buckets:

* *uploads*: used for uploaded user code input/output archives. Currently only .zip files
  are suported

Presigned URLs
--------------
Users (or GUI clients) can upload files directly to the object storage system. To avoid
authentication hassly, presigned URLs can be created that can be used directly to safely
*PUT* files.

.. autofunction:: nomad.files.get_presigned_upload_url
.. autofunction:: nomad.files.create_curl_upload_cmd

Uploads
-------
.. autoclass:: Upload
    :members:

"""
from typing import Callable, List, Any, Generator, IO, TextIO, cast, AnyStr
import os
import os.path
from zipfile import ZipFile, BadZipFile
import shutil
from contextlib import contextmanager
import gzip
import io
import shutil

from nomad import config, utils

# _client = None

# if _client is None and 'sphinx' not in sys.modules:
#     _client = Minio('%s:%s' % (config.minio.host, config.minio.port),
#                     access_key=config.minio.accesskey,
#                     secret_key=config.minio.secret,
#                     secure=False)

#     # ensure all neccessary buckets exist
#     def ensure_bucket(name):
#         try:
#             _client.make_bucket(bucket_name=name)
#             logger.info('Created uploads bucket', bucket=name)
#         except minio.error.BucketAlreadyOwnedByYou:
#             pass

#     ensure_bucket(config.files.uploads_bucket)
#     ensure_bucket(config.files.archive_bucket)


# def get_presigned_upload_url(upload_id: str) -> str:
#     """
#     Generates a presigned upload URL. Presigned URL allows users (and their client programs)
#     to safely *PUT* a single file without further authorization or API to the *uploads* bucket
#     using the given ``upload_id``. Example usages for presigned URLs include
#     browser based uploads or simple *curl* commands (see also :func:`create_curl_upload_cmd`).

#     Arguments:
#         upload_id: The upload id for the uploaded file.

#     Returns:
#         The presigned URL string.
#     """
#     return _client.presigned_put_object(config.files.uploads_bucket, upload_id)


# def create_curl_upload_cmd(presigned_url: str, file_dummy: str='<ZIPFILE>') -> str:
#     """Creates a readymade curl command for uploading.

#     Arguments:
#         presigned_url: The presigned URL to base the command on.

#     Kwargs:
#         file_dummy: A placeholder for the file that the user/client has to replace.

#     Returns:
#         The curl shell command with correct method, url, headers, etc.
#     """
#     return 'curl "%s" --upload-file %s' % (presigned_url, file_dummy)


# def upload_put_handler(func: Callable[[str], None]) -> Callable[[], None]:
#     def upload_notifications(events: List[Any]) -> Generator[str, None, None]:
#         for event in events:
#             for event_record in event['Records']:
#                 try:
#                     event_name = event_record['eventName']
#                     if event_name == 's3:ObjectCreated:Put':
#                         upload_id = event_record['s3']['object']['key']
#                         logger.debug('Received bucket upload event', upload_id=upload_id)
#                         yield upload_id
#                         break  # only one per record, pls
#                     else:
#                         logger.debug('Unhanled bucket event', bucket_event_name=event_name)
#                 except KeyError:
#                     logger.warning(
#                         'Unhandled bucket event due to unexprected event format',
#                         bucket_event_record=event_record)

#     def wrapper(*args, **kwargs) -> None:
#         logger.info('Start listening to uploads notifications.')

#         _client.remove_all_bucket_notification(config.files.uploads_bucket)
#         events = _client.listen_bucket_notification(
#             config.files.uploads_bucket,
#             events=['s3:ObjectCreated:*'])

#         upload_ids = upload_notifications(events)
#         for upload_id in upload_ids:
#             try:
#                 func(upload_id)
#             except StopIteration:
#                 # Using StopIteration to allow clients to stop handling of events.
#                 logger.debug('Handling of upload notifications was stopped via StopIteration.')
#                 return
#             except Exception:
#                 pass

#     return wrapper

class Objects:
    """
    Object store like abstraction based on a regular file system.
    """
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

        return path

    @classmethod
    def open(cls, bucket: str, name: str, ext: str=None, *args, **kwargs) -> IO:
        """ Open an object like you would a file, e.g. with 'rb', etc. """
        try:
            return open(cls._os_path(bucket, name, ext), *args, **kwargs)
        except FileNotFoundError:
            raise KeyError()

    @classmethod
    def delete(cls, bucket: str, name: str, ext: str=None) -> None:
        """ Delete a single object. """
        try:
            os.remove(cls._os_path(bucket, name, ext))
        except FileNotFoundError:
            raise KeyError()

    @classmethod
    def delete_all(cls, bucket: str, prefix: str=''):
        """ Delete all files with given prefix, prefix must denote a directory. """
        try:
            shutil.rmtree(cls._os_path(bucket, prefix, ext=None))
        except FileNotFoundError:
            pass

    @classmethod
    def exists(cls, bucket: str, name: str, ext: str=None) -> bool:
        """ Returns True if object exists. """
        return os.path.exists(cls._os_path(bucket, name, ext))


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
    def __init__(self, bucket: str, object_id: str, ext: str=None) -> None:
        self.bucket = bucket
        self.object_id = object_id
        self.ext = ext

        self.logger = utils.get_logger(__name__, bucket=bucket, object=object_id)

    def open(self, *args, **kwargs) -> IO:
        """ Opens the object with he given mode, etc. """
        self.logger.debug('open file')
        return Objects.open(self.bucket, self.object_id, self.ext, *args, **kwargs)

    def delete(self) -> None:
        """ Deletes the file with the given object id. """
        try:
            Objects.delete(self.bucket, self.object_id, self.ext)
            self.logger.debug('file deleted')
        except FileNotFoundError:
            raise KeyError()

    def exists(self) -> bool:
        """ Returns true if object exists. """
        return Objects.exists(self.bucket, self.object_id, self.ext)

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

    Attributes:
        upload_extract_dir: The path of the tmp directory with the extracted contents.
        filelist: A list of filenames relative to the .zipped upload root.
    """
    def __init__(self, upload_id: str) -> None:
        super().__init__(
            bucket=config.files.uploads_bucket,
            object_id=upload_id,
            ext='zip')

        self.upload_extract_dir: str = os.path.join(config.fs.tmp, 'uploads_extracted', upload_id)
        self.filelist: List[str] = None

    # There is not good way to capsule decorators in a class:
    # https://medium.com/@vadimpushtaev/decorator-inside-python-class-1e74d23107f6
    class Decorators:
        @classmethod
        def handle_errors(cls, decorated):
            def wrapper(self, *args, **kwargs):
                try:
                    return decorated(self, *args, **kwargs)
                except Exception as e:
                    msg = 'Could not %s upload %s.' % (decorated.__name__, self.upload_id)
                    self.logger.error(msg, exc_info=e)
                    raise FileError(msg, e)
            return wrapper

    @Decorators.handle_errors
    def hash(self) -> str:
        """ Calculates the first 28 bytes of a websafe base64 encoded SHA512 of the upload. """
        with self.open('rb') as f:
            return utils.hash(f)

    @Decorators.handle_errors
    def extract(self) -> None:
        """
        'Opens' the upload. This means the uploaed files get extracted to tmp.

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


class ArchiveFile(File):
    """
    Represents the archive file for an individual calculation. Allows to write the
    archive, read the archive, delete the archive. Archive files are stored in
    their own *bucket*.
    """
    def __init__(self, archive_id: str) -> None:
        super().__init__(
            bucket=config.files.archive_bucket,
            object_id=archive_id,
            ext='json.gz' if config.files.compress_archive else 'json')

    @contextmanager
    def write_archive_json(self) -> Generator[TextIO, None, None]:
        """ Context manager that yiels a file-like to write the archive json. """
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

    @contextmanager
    def read_archive_json(self) -> Generator[TextIO, None, None]:
        """ Context manager that yiels a file-like to read the archive json. """
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

    @staticmethod
    def delete_archives(upload_hash: str):
        """ Delete all archives of one upload with the given hash. """
        bucket = config.files.archive_bucket
        Objects.delete_all(bucket, upload_hash)
