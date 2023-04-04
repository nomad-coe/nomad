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

'''
Contains classes and functions to create and maintain file structures
for uploads, and some generic file utilities.

There are two different structures for uploads in two different states: *staging* and *public*.
Possible operations on uploads differ based on this state. Staging is used for
processing, heavily editing, creating hashes, etc. Public is supposed to be a
almost readonly (beside metadata) storage.

.. code-block:: sh

    fs/staging/<upload>/raw/**
                       /archive/<entry_id>.msg
    fs/public/<upload>/raw-{access}.plain.zip
                      /archive-{access}.msg.msg

Where `access` is either "public" (non-embargoed) or "restricted" (embargoed).

There is an implicit relationship between files, based on them being in the same
directory. Each directory with at least one *mainfile* is an *entry directory*
and all the files are *aux* files to that mainfile. This is independent of whether the
respective files actually contributes data or not. An entry directory might
contain multiple mainfiles. E.g., user simulated multiple states of the same system, have
one entry based on the other, etc. In this case the other mainfile is an *aux file* to the
original mainfile, and vice versa.
'''

from abc import ABCMeta
from typing import IO, Set, Dict, Iterable, Iterator, List, Tuple, Any, NamedTuple
from functools import lru_cache
from pydantic import BaseModel
from datetime import datetime
import os.path
import os
import shutil
import tarfile
import zipstream
import hashlib
import io
import json
import yaml
import magic
import zipfile

from nomad import config, utils, datamodel
from nomad.config.models import BundleImportSettings, BundleExportSettings
from nomad.archive import write_archive, read_archive, ArchiveReader


decompress_file_extensions = ('.zip', '.tgz', '.gz', '.tar.gz', '.tar.bz2', '.tar', '.eln')
bundle_info_filename = 'bundle_info.json'

# Used to check if zip-files/archive files are empty
# TODO: These should not be needed when we move on to only keep files with the right access
empty_zip_file_size = 22
empty_archive_file_size = 32


def auto_decompress(path: str):
    '''
    Returns the decompression format ('zip', 'tar' or 'error') if `path` specifies a file
    which should be automatically decompressed before adding it to an upload. If `path`
    does *not* specify a file which we think should be decompressed, or if it specifies a
    directory, we return None.

    The value 'error' means that we think this file should be decompressed, but that we cannot
    decompress it. This indicates that the file is corrupted, has a bad file format, or has
    the wrong file extension.

    Note, some files, like for example excel files, are actually zip files, and we don't want
    to extract such files. Therefore, we only auto decompress if the file has an extension
    we recognize as decompressable, like ".zip", ".tar" etc.
    '''
    if os.path.isdir(path):
        return None
    basename_lower = os.path.basename(path).lower()
    for extension in decompress_file_extensions:
        if basename_lower.endswith(extension):
            if tarfile.is_tarfile(path):
                return 'tar'
            elif zipfile.is_zipfile(path):
                return 'zip'
            return 'error'
    return None


def copytree(src, dst):
    '''
    A close on ``shutils.copytree`` that does not try to copy the stats on all files.
    This is unecessary for our usecase and also causes permission denies for unknown
    reasons.
    '''
    os.makedirs(dst, exist_ok=False)

    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            copytree(s, d)
        else:
            shutil.copyfile(s, d)


def create_tmp_dir(prefix: str) -> str:
    '''
    Creates a temporary directory in the directory specified by `config.fs.tmp`. The name
    of the directory will first be set to `prefix`, but if that name is already taken, a
    suffix will be added to ensure a completely clean, new directory is created. If prefix
    contains a '/', it will be replaced with '_', to ensure the validity of the path.
    The full path to the created directory is returned.
    '''
    assert prefix
    prefix = prefix.replace(os.path.sep, '_')
    assert is_safe_basename(prefix)
    for index in range(1, 100):
        dir_name = prefix if index == 1 else f'{prefix}_{index}'
        path = os.path.join(config.fs.tmp, dir_name)
        try:
            os.makedirs(path)
            return path
        except FileExistsError:
            pass  # Try again with different suffix
    raise RuntimeError('Could not create temporary directory - too many directories with same prefix?')


def is_safe_basename(basename: str) -> bool:
    '''
    Checks if `basename` is a *safe* base name (file/folder name). We consider it safe if
    it is not empty, does not contain any '/', and is not equal to '.' or '..'
    '''
    if not basename or '/' in basename or basename == '.' or basename == '..':
        return False
    return True


def is_safe_relative_path(path: str) -> bool:
    '''
    Checks if path is a *safe* relative path. We consider it safe if it does not start with
    '/' or use '.' or '..' elements (which could be open for security leaks if allowed).
    It may end with a single '/', indicating that a folder is referred. For referring to
    the base folder, the empty string should be used (not '.' etc).
    '''
    if type(path) != str:
        return False
    if path == '':
        return True
    if path.startswith('/') or '//' in path or '\n' in path:
        return False
    for element in path.split('/'):
        if element == '.' or element == '..':
            return False
    return True


class PathObject:
    '''
    Object storage-like abstraction for paths in general.
    Attributes:
        os_path: The full os path of the object.
    '''
    def __init__(self, os_path: str):
        self.os_path = os_path

    def delete(self) -> None:
        if os.path.isfile(self.os_path):
            os.remove(self.os_path)
        else:
            shutil.rmtree(self.os_path)

    def exists(self) -> bool:
        return os.path.exists(self.os_path)

    @property
    def size(self) -> int:
        ''' The os determined file size. '''
        return os.stat(self.os_path).st_size

    def __repr__(self) -> str:
        return self.os_path


class DirectoryObject(PathObject):
    '''
    Object storage-like abstraction for directories.
    '''
    def __init__(self, os_path: str, create: bool = False):
        self.os_path = os_path
        if create and not os.path.isdir(self.os_path):
            os.makedirs(self.os_path)

    def join_dir(self, path, create: bool = False) -> 'DirectoryObject':
        return DirectoryObject(os.path.join(self.os_path, path), create)

    def join_file(self, path, create_dir: bool = False) -> PathObject:
        if create_dir:
            dirname = os.path.dirname(path)
            if dirname:
                dir_os_path = os.path.join(self.os_path, dirname)
                if not os.path.exists(dir_os_path):
                    os.makedirs(dir_os_path)
        return PathObject(os.path.join(self.os_path, path))

    def exists(self) -> bool:
        return os.path.isdir(self.os_path)


class RawPathInfo(NamedTuple):
    '''
    Stores basic info about a file or folder located at a specific raw path.
    '''
    path: str
    is_file: bool
    size: int
    access: str


class StreamedFile(BaseModel):
    '''
    Convenience class for representing a streamed file, together with information about
    file size and an associated path.
    '''
    f: Any
    path: str
    size: int


class FileSource(metaclass=ABCMeta):
    '''
    An abstract class which represents a generic "file source", from which some number of files
    can be retrieved. There are several different ways to create a file source, see subclasses.
    The files in the source are associated with paths and have known sizes.
    '''
    def to_streamed_files(self) -> Iterable[StreamedFile]:
        '''
        Retrieves the files in the source as :class:`StreamedFile` objects.
        The caller should close the streams when consumed.
        '''
        raise NotImplementedError()

    def to_zipstream(self) -> Iterator[bytes]:
        ''' Returns a zip stream with the files from this FileSource. '''
        return create_zipstream(self.to_streamed_files())

    def to_zipfile(self, path, overwrite: bool = False):
        '''
        Generates a zip file from the files in this FileSource and stores it to disk. The
        zipfile content is created by calling :func:`to_zipstream`.
        '''
        assert not os.path.isdir(path), 'Exporting to zip file requires a file path, not directory.'
        assert overwrite or not os.path.exists(path), '`path` already exists. Use `overwrite` to overwrite.'
        with open(path, 'wb') as f:
            for chunk in self.to_zipstream():
                f.write(chunk)

    def to_disk(self, destination_dir: str, move_files: bool = False, overwrite: bool = False):
        '''
        Writes the files from this FileSource to disk, uncompressed. The default implementation
        makes use of :func:`to_streamed_files`. The `destination_dir` should be a directory
        (it will be created if it does not exist). The `move_files` argument instructs
        the method to move the source files if possible.
        '''
        if not os.path.exists(destination_dir):
            os.makedirs(destination_dir)
        assert os.path.isdir(destination_dir), '`destination_dir` is not a directory'
        for streamed_file in self.to_streamed_files():
            assert is_safe_relative_path(streamed_file.path), 'Unsafe relative path encountered'
            os_path = os.path.join(destination_dir, streamed_file.path)
            dir_path = os.path.dirname(os_path)
            if os.path.exists(os_path):
                assert overwrite, 'Target already exists and `overwrite` is False'
                PathObject(os_path).delete()
            os.makedirs(dir_path, exist_ok=True)
            with open(os_path, 'wb') as output_file:
                with streamed_file.f:
                    for chunk in streamed_file.f:
                        output_file.write(chunk)

    def close(self):
        ''' Perform "closing" of the source, if applicable. '''
        pass


class BrowsableFileSource(FileSource, metaclass=ABCMeta):
    '''
    A :class:`FileSource` which can be "browsed", like a folder on disk or a zip archive.
    '''
    def open(self, path, mode='rb') -> IO:
        ''' Opens a file by the specified path. '''
        raise NotImplementedError()

    def directory_list(self, path: str) -> List[str]:
        '''
        Returns a list of directory contents, located in the directory denoted by `path`
        in this file source.
        '''
        raise NotImplementedError()

    def sub_source(self, path: str) -> 'BrowsableFileSource':
        '''
        Creates a new instance of :class:`BrowsableFileSource` which just contains the
        files located under the specified path.
        '''
        raise NotImplementedError()


class StreamedFileSource(FileSource):
    '''
    A :class:`FileSource` created from a single :class:`StreamedFile`.
    '''
    def __init__(self, streamed_file: StreamedFile):
        self.streamed_file = streamed_file

    def to_streamed_files(self) -> Iterable[StreamedFile]:
        yield self.streamed_file


class DiskFileSource(BrowsableFileSource):
    '''
    A :class:`FileSource` corresponding to a single file or a folder on disk. The object
    is identified by a `base_path` and a `relative path`. The `base_path` should be a folder,
    the `relative_path` is optional, and used for selecting only a specific file or folder
    located under `base_folder`. The paths of the files retrieved from this source are given
    relative to the `base_path`.
    '''
    def __init__(self, base_path: str, relative_path: str = None):
        assert os.path.isdir(base_path)
        if relative_path:
            assert is_safe_relative_path(relative_path), 'Unsafe relative_path received'
            self.full_path = os.path.join(base_path, relative_path)
            assert os.path.exists(self.full_path)
        else:
            self.full_path = base_path
        self.base_path = base_path
        self.relative_path = relative_path

    def to_streamed_files(self) -> Iterable[StreamedFile]:
        if os.path.isfile(self.full_path):
            # Single file
            yield StreamedFile(
                path=self.relative_path,
                f=open(self.full_path, 'rb'),
                size=os.stat(self.full_path).st_size)
        else:
            # Directory - crawl it and its subfolders for files
            for dirpath, __, filenames in os.walk(self.full_path):
                for filename in filenames:
                    sub_full_path = os.path.join(dirpath, filename)
                    sub_relative_path = os.path.relpath(sub_full_path, self.base_path)
                    yield StreamedFile(
                        path=sub_relative_path,
                        f=open(sub_full_path, 'rb'),
                        size=os.stat(sub_full_path).st_size)

    def to_disk(self, destination_dir: str, move_files: bool = False, overwrite: bool = False):
        if self.relative_path:
            destination_path = os.path.join(destination_dir, self.relative_path)
        else:
            destination_path = destination_dir
        destination_parent = os.path.dirname(destination_path)
        os.makedirs(destination_parent, exist_ok=True)
        if os.path.exists(destination_path):
            assert overwrite, f'Target {destination_path} already exists and `overwrite` is False'
            PathObject(destination_path).delete()
        # All looks good. Copy or move the source to the destination
        if move_files:
            shutil.move(self.full_path, destination_path)
        else:
            if os.path.isfile(self.full_path):
                shutil.copyfile(self.full_path, destination_path)
            else:
                copytree(self.full_path, destination_path)

    def open(self, path, mode='rb') -> IO:
        assert is_safe_relative_path(path)
        return open(os.path.join(self.base_path, path), mode)

    def directory_list(self, path: str) -> List[str]:
        assert is_safe_relative_path(path)
        sub_path = os.path.join(self.base_path, path)
        return os.listdir(sub_path)

    def sub_source(self, path: str) -> 'DiskFileSource':
        assert is_safe_relative_path(path)
        return DiskFileSource(self.base_path, path)


class ZipFileSource(BrowsableFileSource):
    '''
    Allows us to "wrap" a :class:`zipfile.ZipFile` object and use it as a :class:`BrowsableFileSource`,
    i.e. it denotes a resource (single file or folder) stored in a ZipFile.
    '''
    def __init__(self, zip_file: zipfile.ZipFile, sub_path: str = ''):
        assert is_safe_relative_path(sub_path)
        self.zip_file = zip_file
        self.sub_path = sub_path
        self._namelist: List[str] = zip_file.namelist()

    def to_streamed_files(self) -> Iterable[StreamedFile]:
        path_prefix = '' if not self.sub_path else self.sub_path + os.path.sep
        for path in self._namelist:
            if path == self.sub_path or (path.startswith(path_prefix) and not path.endswith(os.path.sep)):
                yield StreamedFile(
                    path=path,
                    f=self.zip_file.open(path, 'r'),
                    size=self.zip_file.getinfo(path).file_size)

    def open(self, path, mode='rb') -> IO:
        assert 'r' in mode, 'Mode must be a read mode'
        for c in mode:
            assert c in ('r', 'b', 't'), f'Invalid mode for open command: {mode}'
        f = self.zip_file.open(path, 'r')
        if 't' in mode:
            return io.TextIOWrapper(f)
        return f

    def directory_list(self, path: str) -> List[str]:
        path_prefix = '' if not path else path + os.path.sep
        found = set()
        for path2 in self._namelist:
            if path2.startswith(path_prefix):
                found.add(path2.split(os.path.sep)[0])
        return sorted(found)

    def sub_source(self, path: str) -> 'ZipFileSource':
        assert is_safe_relative_path(path), 'Unsafe path provided'
        if self.sub_path:
            assert path.startswith(self.sub_path + os.path.sep), 'Provided `path` is not a sub path.'
        return ZipFileSource(self.zip_file, path)

    def close(self):
        self.zip_file.close()


class CombinedFileSource(FileSource):
    '''
    Class for defining a :class:`FileSource` by combining multiple "subsources" into one.
    '''
    def __init__(self, file_sources=Iterable[FileSource]):
        ''' file_sources: an Iterable for getting FileSources. '''
        self.file_sources = file_sources

    def to_streamed_files(self) -> Iterable[StreamedFile]:
        for file_source in self.file_sources:
            for streamed_file in file_source.to_streamed_files():
                yield streamed_file

    def to_disk(self, destination_dir: str, move_files: bool = False, overwrite: bool = False):
        for file_source in self.file_sources:
            file_source.to_disk(destination_dir, move_files, overwrite)


class StandardJSONEncoder(json.JSONEncoder):
    """ Our standard JSONEncoder with support for marshalling of datetime objects """
    def default(self, obj):  # pylint: disable=E0202
        if isinstance(obj, datetime):
            return {'$datetime': obj.isoformat()}
        return json.JSONEncoder.default(self, obj)


class StandardJSONDecoder(json.JSONDecoder):
    """ Our standard JSONDecoder, with support for marshalling of datetime objects """
    def __init__(self, *args, **kargs):
        json.JSONDecoder.__init__(self, object_hook=self.dict_to_object, *args, **kargs)

    def dict_to_object(self, d):
        v = d.get('$datetime')
        if v is not None and len(d) == 1:
            return datetime.fromisoformat(v)
        return d


def json_to_streamed_file(json_dict: Dict[str, Any], path: str) -> StreamedFile:
    ''' Converts a json dictionary structure to a :class:`StreamedFile`. '''
    json_bytes = json.dumps(json_dict, indent=2, cls=StandardJSONEncoder).encode()
    return StreamedFile(
        path=path,
        f=io.BytesIO(json_bytes),
        size=len(json_bytes))


def create_zipstream_content(streamed_files: Iterable[StreamedFile]) -> Iterable[Dict]:
    '''
    Generator which "casts" a sequence of StreamedFiles to a sequence of dictionaries, of
    the form which is required by the `zipstream` library, i.e. dictionaries with keys
    `arcname`, `iterable` and `buffer_size`. Useful for generating zipstreams.
    '''
    for streamed_file in streamed_files:

        def content_generator():
            with streamed_file.f as f:
                while True:
                    data = f.read(1024 * 64)
                    if not data:
                        break
                    yield data

        yield dict(
            arcname=streamed_file.path,
            iterable=content_generator(),
            buffer_size=streamed_file.size)


def create_zipstream(
        streamed_files: Iterable[StreamedFile],
        compress: bool = False) -> Iterator[bytes]:
    '''
    Creates a zip stream, i.e. a streamed zip file.
    '''
    compression = zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED
    zip_stream = zipstream.ZipFile(mode='w', compression=compression, allowZip64=True)
    zip_stream.paths_to_write = create_zipstream_content(streamed_files)

    return iter(zip_stream)


class UploadFiles(DirectoryObject, metaclass=ABCMeta):
    ''' Abstract base class for upload files. '''
    def __init__(self, upload_id: str, create: bool = False):
        self.logger = utils.get_logger(__name__, upload_id=upload_id)

        super().__init__(os_path=self.base_folder_for(upload_id), create=create)

        if not create and not self.exists():
            raise KeyError(upload_id)

        self.upload_id = upload_id

    @classmethod
    def file_area(cls):
        '''
        Full path to where the upload files of this class are stored (i.e. either
        staging or public file area).
        '''
        raise NotImplementedError()

    @property
    def external_os_path(self):
        '''
        Full path to where the upload files of this class are stored on the server.
        This is equal to `self.os_path` if no external path substitutes for staging
        and public area are configured. This is helpful, when nomad is run in a container
        and the mounted path used by nomad are different from the actual paths on the
        host server.
        '''
        raise NotImplementedError()

    @classmethod
    def base_folder_for(cls, upload_id: str) -> str:
        '''
        Full path to the base folder for the upload files (of this class) for the
        specified upload_id.
        '''
        os_path = cls.file_area()
        if config.fs.prefix_size:
            os_path = os.path.join(os_path, upload_id[:config.fs.prefix_size])
        os_path = os.path.join(os_path, upload_id)
        return os_path

    @classmethod
    def exists_for(cls, upload_id: str) -> bool:
        '''
        If an UploadFiles object (of this class) has been created for this upload_id.
        '''
        return os.path.exists(cls.base_folder_for(upload_id))

    def to_staging_upload_files(self, create: bool = False, include_archive: bool = False) -> 'StagingUploadFiles':
        ''' Casts to or creates corresponding staging upload files or returns None. '''
        raise NotImplementedError()

    @staticmethod
    def get(upload_id: str, create: bool = False) -> 'UploadFiles':
        if StagingUploadFiles.exists_for(upload_id):
            return StagingUploadFiles(upload_id, create)
        elif PublicUploadFiles.exists_for(upload_id):
            return PublicUploadFiles(upload_id, create)
        else:
            return None

    def is_empty(self) -> bool:
        ''' If this upload has no content yet. '''
        raise NotImplementedError()

    def raw_path_exists(self, path: str) -> bool:
        '''
        Returns True if the specified path is a valid raw path (either file or directory)
        '''
        raise NotImplementedError()

    def raw_path_is_file(self, path: str) -> bool:
        '''
        Returns True if the specified path points to a file (rather than a directory).
        '''
        raise NotImplementedError()

    def raw_directory_list(
            self, path: str = '', recursive=False, files_only=False, path_prefix=None) -> Iterable[RawPathInfo]:
        '''
        Returns an iterable of RawPathInfo, one for each element (file or folder) in
        the directory specified by `path`. If `recursive` is set to True, subdirectories are
        also crawled. If `files_only` is set, only the file objects found are returned.
        If path is not a valid directory, the result will be empty. Selecting empty string
        as path (which is the default value) gives the content of the whole raw directory.
        The `path_prefix` argument can be used to filter out elements where the path starts
        with a specific prefix.
        '''
        raise NotImplementedError()

    def raw_file(self, file_path: str, *args, **kwargs) -> IO:
        '''
        Opens a raw file and returns a file-like object. Additional args, kwargs are
        delegated to the respective `open` call.
        Arguments:
            file_path: The path to the file relative to the upload.
        Raises:
            KeyError: If the file does not exist.
        '''
        raise NotImplementedError()

    def raw_file_size(self, file_path: str) -> int:
        '''
        Returns:
            The size of the given raw file.
        '''
        raise NotImplementedError()

    def raw_file_mime_type(self, file_path: str) -> str:
        assert self.raw_path_is_file(file_path), 'Provided path does not specify a file, or is invalid.'
        raw_file = self.raw_file(file_path, 'br')
        buffer = raw_file.read(2048)
        mime_type = magic.from_buffer(buffer, mime=True)
        raw_file.close()
        if not mime_type:
            mime_type = 'application/octet-stream'
        return mime_type

    def read_archive(self, entry_id: str, use_blocked_toc: bool = True) -> ArchiveReader:
        '''
        Returns an :class:`nomad.archive.ArchiveReader` that contains the
        given entry_id.
        '''
        raise NotImplementedError()

    def close(self):
        ''' Release possibly held system resources (e.g. file handles). '''
        pass

    def delete(self) -> None:
        shutil.rmtree(self.os_path, ignore_errors=True)
        if config.fs.prefix_size > 0:
            # If using prefix, also remove the parent directory if empty
            parent_directory = os.path.dirname(self.os_path)
            if not os.listdir(parent_directory):
                try:
                    os.rmdir(parent_directory)
                except Exception as e:
                    utils.get_logger(__name__).error(
                        'could not remove empty prefix dir', directory=parent_directory, exc_info=e)

    def files_to_bundle(self, export_settings: BundleExportSettings) -> Iterable[FileSource]:
        '''
        A generator of :class:`FileSource` objects, defining the files/folders to be included in an
        upload bundle when *exporting*. The arguments allows for further filtering of what to include.

        Note, this only yields files to copy from the regular upload directory, not "special" files,
        like the bundle_info.json file, which is created by the :class:`BundleExporter`.
        '''
        raise NotImplementedError()

    @classmethod
    def files_from_bundle(
            cls, bundle_file_source: BrowsableFileSource, import_settings: BundleImportSettings) -> Iterable[FileSource]:
        '''
        Returns an Iterable of :class:`FileSource`, defining the files/folders to be included in an
        upload bundle when *importing*. Only the files specified by the import_settings are included.
        '''
        raise NotImplementedError()


class StagingUploadFiles(UploadFiles):
    def __init__(self, upload_id: str, create: bool = False):
        super().__init__(upload_id, create)

        self._raw_dir = self.join_dir('raw', create)
        self._archive_dir = self.join_dir('archive', create)
        self._frozen_file = self.join_file('.frozen')

        self._size = 0

    @classmethod
    def file_area(cls):
        return config.fs.staging

    @property
    def external_os_path(self):
        if not config.fs.staging_external:
            return self.os_path

        return self.os_path.replace(config.fs.staging, config.fs.staging_external)

    def to_staging_upload_files(self, create: bool = False, include_archive: bool = False) -> 'StagingUploadFiles':
        return self

    @property
    def size(self) -> int:
        return self._size

    def _file(self, path_object: PathObject, *args, **kwargs) -> IO:
        try:
            return open(path_object.os_path, *args, **kwargs)
        except FileNotFoundError:
            raise KeyError(path_object.os_path)
        except IsADirectoryError:
            raise KeyError(path_object.os_path)

    def is_empty(self) -> bool:
        return not os.path.exists(self._raw_dir.os_path) or not os.listdir(self._raw_dir.os_path)

    def raw_path_exists(self, path: str) -> bool:
        if not is_safe_relative_path(path):
            return False
        return os.path.exists(os.path.join(self._raw_dir.os_path, path))

    def raw_path_is_file(self, path: str) -> bool:
        if not is_safe_relative_path(path):
            return False
        return os.path.isfile(os.path.join(self._raw_dir.os_path, path))

    def raw_create_directory(self, path: str):
        assert path and is_safe_relative_path(path), 'Bad path provided'
        os.makedirs(os.path.join(self._raw_dir.os_path, path))

    def raw_directory_list(
            self, path: str = '', recursive=False, files_only=False, path_prefix=None) -> Iterable[RawPathInfo]:
        if not is_safe_relative_path(path):
            return
        os_path = os.path.join(self._raw_dir.os_path, path)
        if not os.path.isdir(os_path):
            is_file = os.path.isfile(os_path)
            if is_file:
                yield RawPathInfo(path=path, is_file=True, size=os.stat(os_path).st_size, access='unpublished')
            return
        for element_name in sorted(os.listdir(os_path)):
            element_raw_path = os.path.join(path, element_name)
            element_os_path = os.path.join(os_path, element_name)
            is_file = os.path.isfile(element_os_path)
            if not is_file:
                # Crawl sub directory.
                dir_size = 0
                for sub_path_info in self.raw_directory_list(element_raw_path, True, files_only):
                    if sub_path_info.is_file:
                        dir_size += sub_path_info.size
                    if recursive:
                        if not path_prefix or sub_path_info.path.startswith(path_prefix):
                            yield sub_path_info

            if not files_only or is_file:
                size = os.stat(element_os_path).st_size if is_file else dir_size
                if not path_prefix or element_raw_path.startswith(path_prefix):
                    yield RawPathInfo(
                        path=element_raw_path,
                        is_file=is_file,
                        size=size,
                        access='unpublished')

    def raw_file(self, file_path: str, *args, **kwargs) -> IO:
        assert is_safe_relative_path(file_path)
        return self._file(self.raw_file_object(file_path), *args, **kwargs)

    def raw_file_size(self, file_path: str) -> int:
        assert is_safe_relative_path(file_path)
        return self.raw_file_object(file_path).size

    def raw_file_object(self, file_path: str) -> PathObject:
        assert is_safe_relative_path(file_path)
        return self._raw_dir.join_file(file_path)

    def write_archive(self, entry_id: str, data: Any) -> int:
        ''' Writes the data as archive file and returns the archive file size. '''
        archive_file_object = self.archive_file_object(entry_id)
        try:
            write_archive(archive_file_object.os_path, 1, data=[(entry_id, data)])
        except Exception as e:
            # in case of failure, remove the possible corrupted archive file
            if archive_file_object.exists():
                archive_file_object.delete()

            raise e

        return self.archive_file_object(entry_id).size

    def read_archive(self, entry_id: str, use_blocked_toc: bool = True) -> ArchiveReader:
        try:
            return read_archive(self.archive_file_object(entry_id).os_path, use_blocked_toc=use_blocked_toc)

        except FileNotFoundError:
            raise KeyError(entry_id)

    def archive_file_object(self, entry_id: str) -> PathObject:
        version_suffix = '-' + config.fs.archive_version_suffix if config.fs.archive_version_suffix else ''
        return self._archive_dir.join_file(f'{entry_id}{version_suffix}.msg')

    def add_rawfiles(
            self, path: str, target_dir: str = '', cleanup_source_file_and_dir: bool = False,
            updated_files: Set[str] = None) -> None:
        '''
        Adds the file or folder specified by `path` to this upload, in the raw directory
        specified by `target_dir`. If `path` denotes a zip or tar archive file, it will
        first be extracted to a temporary directory. The file(s) are *merged* with the
        existing upload files, i.e. new files are added, replacing old files if there
        already exists file(s) by the same names, the rest of the old files are left
        untouched.

        Cleanup
        The method is responsible for trying to clean up temporarily extracted files.
        If `cleanup_source_file_and_dir` is True, the source file (defined by `path`), and
        its parent directory (which we also assume is temporary) are also cleaned up.
        Note: the cleanup steps are always carried out, also if the operation fails.

        Arguments:
            path: OS path to a file or folder to add.
            target_dir: A raw path (i.e. path relative to the raw directory) defining
                where the resource defined by `path` should be put. If `target_dir` is not
                specified, it defaults to the empty string, i.e. the upload's raw dir.
            cleanup_source_file_and_dir: If true, the source file/folder (defined by `path`) is
                deleted when done (regardless of success or failure). Additionally, the parent
                folder is deleted if it's empty or if the operation failed. Use when the file/folder
                to add is stored in a temporary directory.
            updated_files: An optional set of paths. If provided with the call, the raw
                path of all files added or updated by the operation will be added to this set.
        '''
        tmp_dir = None
        try:
            assert not self.is_frozen
            assert os.path.exists(path), f'{path} does not exist'
            assert is_safe_relative_path(target_dir)
            self._size += os.stat(path).st_size

            is_dir = os.path.isdir(path)
            decompress = auto_decompress(path)
            if decompress:
                tmp_dir = create_tmp_dir(self.upload_id + '_unzip')
                if decompress == 'zip':
                    with zipfile.ZipFile(path) as zf:
                        zf.extractall(tmp_dir)
                elif decompress == 'tar':
                    with tarfile.open(path) as tf:
                        tf.extractall(tmp_dir)
                elif decompress == 'error':
                    # Unknown / bad file format
                    assert False, 'Cannot extract file. Bad file format or file extension?'

            # Determine what to merge
            elements_to_merge: Iterable[Tuple[str, List[str], List[str]]] = []
            if is_dir:
                # Directory
                source_dir = path
                elements_to_merge = os.walk(source_dir)
            elif decompress:
                # Zipped archive
                source_dir = tmp_dir
                elements_to_merge = os.walk(source_dir)
            else:
                # Single, non-compressed file
                source_dir = os.path.dirname(path)
                elements_to_merge = [(source_dir, [], [os.path.basename(path)])]

            # Do the merge
            os_target_dir = os.path.join(self._raw_dir.os_path, target_dir)
            if not os.path.isdir(os_target_dir):
                os.makedirs(os_target_dir)
            for source_root, dirs, files in elements_to_merge:
                elements = dirs + files
                for element in elements:
                    element_source_path = os.path.join(source_root, element)
                    element_relative_path = os.path.relpath(element_source_path, source_dir)
                    element_target_path = os.path.join(os_target_dir, element_relative_path)
                    if os.path.islink(element_source_path):
                        continue  # Skip links, could pose security risk
                    if os.path.exists(element_target_path):
                        if os.path.isfile(element_target_path) != os.path.isfile(element_source_path):
                            assert False, f'Cannot merge a file with a directory or vice versa: {element_relative_path}'

                    # Copy or move the element
                    if os.path.isdir(element_source_path):
                        # Directory - just create corresponding directory in the target if needed.
                        if not os.path.exists(element_target_path):
                            os.makedirs(element_target_path)
                    else:
                        # File - copy or move it
                        if cleanup_source_file_and_dir or decompress:
                            # Move the file
                            shutil.move(element_source_path, element_target_path)
                        else:
                            # Copy the file
                            shutil.copyfile(element_source_path, element_target_path)
                        if updated_files is not None:
                            updated_files.add(os.path.join(target_dir, element_relative_path))
        except Exception:
            if cleanup_source_file_and_dir:
                parent_dir = os.path.dirname(path)
                if os.path.exists(parent_dir):
                    shutil.rmtree(parent_dir)
            raise
        finally:
            # Cleanup
            if tmp_dir and os.path.exists(tmp_dir):
                shutil.rmtree(tmp_dir)
            if cleanup_source_file_and_dir:
                if os.path.exists(path):
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
                parent_dir = os.path.dirname(path)
                if os.path.exists(parent_dir) and not os.listdir(parent_dir):
                    shutil.rmtree(parent_dir)

    def delete_rawfiles(self, path, updated_files: Set[str] = None):
        assert is_safe_relative_path(path)
        raw_os_path = os.path.join(self.os_path, 'raw')
        os_path = os.path.join(raw_os_path, path)
        if not os.path.exists(os_path):
            return
        if os.path.isfile(os_path):
            # Deleting a file
            if updated_files is not None:
                updated_files.add(path)
            os.remove(os_path)
        else:
            # Deleting a directory
            if updated_files is not None:
                for dirname, _, filenames in os.walk(os_path):
                    for filename in filenames:
                        file_os_path = os.path.join(dirname, filename)
                        file_raw_path = os.path.relpath(file_os_path, raw_os_path)
                        updated_files.add(file_raw_path)
            shutil.rmtree(os_path)
        if path == '':
            # Special case - deleting everything, i.e. the entire raw folder. Need to recreate.
            os.makedirs(os_path)

    def copy_or_move_rawfile(self, path_to_existing_file, path_to_target_file, copy_or_move, updated_files: Set[str] = None):
        assert is_safe_relative_path(path_to_existing_file)
        assert is_safe_relative_path(path_to_target_file)
        os_path_exisitng = os.path.join(self._raw_dir.os_path, path_to_existing_file)
        os_path_target = os.path.join(self._raw_dir.os_path, path_to_target_file)
        if not os.path.exists(os_path_exisitng):
            return
        if os.path.isfile(os_path_exisitng):
            # copying or moving a file
            if os.path.exists(os_path_target):
                raise ValueError('A file with the same name already exists.')
            if copy_or_move.lower() == 'copy':
                shutil.copyfile(os_path_exisitng, os_path_target)
            elif copy_or_move.lower() == 'move':
                shutil.move(os_path_exisitng, os_path_target)

            if updated_files is not None:
                updated_files.add(path_to_target_file)
                # if both the new and old name are the same then no new entry will be
                # added to the set. but if different, we add the old one so that later on
                # when self.matchall is called in data.py, the old filename is removed
                # from mongo database
                updated_files.add(path_to_existing_file)

        elif os.path.isdir(os_path_target):
            raise ValueError('Copying a directory is not possible.')

    @lru_cache()
    def metadata_file_cached(self, path_dir: str = ''):
        '''
        Gets the content of the metadata file located in the directory defined by `path_dir`.
        The `path_dir` should be relative to the `raw` folder.
        '''
        path_incl_filename = os.path.join(path_dir, config.process.metadata_file_name)
        for ext in config.process.metadata_file_extensions:
            path_incl_filename_ext = path_incl_filename + '.' + ext
            full_path = os.path.join(self._raw_dir.os_path, path_incl_filename_ext)
            if os.path.isfile(full_path):
                try:
                    with open(full_path) as f:
                        if full_path.endswith('.json'):
                            return json.load(f)
                        elif full_path.endswith('.yaml') or full_path.endswith('.yml'):
                            return yaml.load(f, Loader=getattr(yaml, 'FullLoader'))
                        else:
                            return {}
                except Exception as e:
                    # ignore the file contents if the file is not parsable, just warn.
                    self.logger.warn(
                        'could not parse nomad.yaml/json', path=path_incl_filename_ext, exc_info=e)
        return {}

    @property
    def is_frozen(self) -> bool:
        ''' Returns True if this upload is already *bagged*. '''
        return self._frozen_file.exists()

    def pack(
            self, entries: List[datamodel.EntryMetadata], with_embargo: bool, create: bool = True,
            include_raw: bool = True, include_archive: bool = True) -> None:
        '''
        Packs raw and/or archive files, to create the contents in the public file area.
        This method should be called when an upload is published, or when a
        published upload has been reprocessed.

        If the public upload files directory does not exist, it will be created.
        If the target archive file or raw file zip exists, they will be overwritten.
        If an archive file or raw file zip with the wrong access exists, they will be deleted.
        This is potentially a long running operation.

        Arguments:
            entries: A list of EntryMetadata to pack in the archive files
            with_embargo: If the upload is embargoed (determines which "access" is used in
                the file names)
            create: if the public upload files directory should be created. True by default.
            include_raw: determines if the raw data should be packed. True by default.
            include_archive: determines of the archive data should be packed. True by default.
        '''
        self.logger.info('started to pack upload')

        # freeze the upload
        assert not self.is_frozen, "Cannot pack an upload that is packed, or packing."
        with open(self._frozen_file.os_path, 'wt') as f:
            f.write('frozen')

        # Check embargo flag consistency
        for entry in entries:
            assert entry.with_embargo == with_embargo

        access = 'restricted' if with_embargo else 'public'
        other_access = 'public' if with_embargo else 'restricted'  # The "inverted" access

        # Get or create a target dir in the public area
        target_dir = DirectoryObject(PublicUploadFiles.base_folder_for(self.upload_id), create=create)
        if os.listdir(target_dir.os_path):
            # Target dir contains files. Check that the target access is identical
            assert PublicUploadFiles(self.upload_id).access == access, 'Inconsistent access'

        # zip archives
        if include_archive:
            with utils.timer(self.logger, 'packed msgpack archive') as log_data:
                number_of_entries = self._pack_archive_files(target_dir, entries, access, other_access)
                log_data.update(number_of_entries=number_of_entries)

        # zip raw files
        if include_raw:
            with utils.timer(self.logger, 'packed raw files'):
                self._pack_raw_files(target_dir, access, other_access)

    def _pack_archive_files(
            self, target_dir: DirectoryObject, entries: List[datamodel.EntryMetadata], access: str, other_access: str):
        number_of_entries = len(entries)

        def create_iterator():
            for entry in entries:
                archive_file = self.archive_file_object(entry.entry_id)
                if archive_file.exists():
                    data = read_archive(archive_file.os_path)[entry.entry_id].to_dict()
                    yield (entry.entry_id, data)
                else:
                    yield (entry.entry_id, {})

        try:
            file_object = PublicUploadFiles._create_msg_file_object(target_dir, access)
            write_archive(file_object.os_path, number_of_entries, create_iterator())
            # Remove the file with the opposite access, if it exists
            other_file_object = PublicUploadFiles._create_msg_file_object(target_dir, other_access)
            if other_file_object.exists():
                other_file_object.delete()  # This file should be empty, if it exists
        except Exception as e:
            self.logger.error('exception during packing archives', exc_info=e)
            raise

        return number_of_entries

    def _pack_raw_files(self, target_dir: DirectoryObject, access: str, other_access: str):
        try:
            raw_zip_file_object = PublicUploadFiles._create_raw_zip_file_object(target_dir, access)
            with zipfile.ZipFile(raw_zip_file_object.os_path, mode='w') as raw_zip:
                for path_info in self.raw_directory_list(recursive=True):
                    basename = os.path.basename(path_info.path)
                    if basename.startswith('POTCAR'):
                        if not basename.endswith('.stripped'):
                            continue  # Skip the unstripped POTCAR files when publishing
                        if basename.endswith('.stripped.stripped'):
                            continue  # Skip redundantly stripped POTCAR files (created due to bug #979) when publishing
                    raw_zip.write(self._raw_dir.join_file(path_info.path).os_path, path_info.path)
            # Remove the zip file with the opposite access, if it exists
            other_raw_zip_file_object = PublicUploadFiles._create_raw_zip_file_object(target_dir, other_access)
            if other_raw_zip_file_object.exists():
                other_raw_zip_file_object.delete()  # This file should be empty, if it exists
        except Exception as e:
            self.logger.error('exception during packing raw files', exc_info=e)
            raise

    def entry_files(self, mainfile: str, with_mainfile: bool = True, with_cutoff: bool = True) -> Iterable[str]:
        '''
        Returns all the auxfiles and mainfile for a given mainfile. This implements
        nomad's logic about what is part of an entry and what not. The mainfile
        is the first element, the rest is sorted.
        Arguments:
            mainfile: The mainfile path relative to upload
            with_mainfile: Do include the mainfile, default is True
        '''
        mainfile_object = self._raw_dir.join_file(mainfile)
        if not mainfile_object.exists():
            raise KeyError(mainfile)

        mainfile_basename = os.path.basename(mainfile)
        entry_dir = os.path.dirname(mainfile_object.os_path)
        entry_relative_dir = entry_dir[len(self._raw_dir.os_path) + 1:]

        file_count = 0
        aux_files: List[str] = []
        dir_elements = os.listdir(entry_dir)
        dir_elements.sort()
        for dir_element in dir_elements:
            if dir_element != mainfile_basename and os.path.isfile(os.path.join(entry_dir, dir_element)):
                aux_files.append(os.path.join(entry_relative_dir, dir_element))
                file_count += 1

            if with_cutoff and file_count > config.process.auxfile_cutoff:
                # If there are too many of them, its probably just a directory with lots of
                # mainfiles/entries. In this case it does not make any sense to provide thousands of
                # aux files.
                break

        aux_files = sorted(aux_files)

        if with_mainfile:
            return [mainfile] + aux_files
        else:
            return aux_files

    def entry_hash(self, mainfile: str, mainfile_key: str) -> str:
        '''
        Calculates a hash for the given entry based on file contents and aux file contents.
        Arguments:
            mainfile: The mainfile path relative to the upload that identifies the entry in
                the folder structure.
            mainfile_key: The mainfile_key of the entry (if any)
        Returns:
            The calculated hash
        Raises:
            KeyError: If the mainfile does not exist.
        '''
        hash = hashlib.sha512()
        for filepath in self.entry_files(mainfile):
            with open(self._raw_dir.join_file(filepath).os_path, 'rb') as f:
                for data in iter(lambda: f.read(65536), b''):
                    hash.update(data)
        if mainfile_key:
            hash.update(mainfile_key.encode('utf8'))
        return utils.make_websave(hash)

    def files_to_bundle(self, export_settings: BundleExportSettings) -> Iterable[FileSource]:
        # Defines files for upload bundles of staging uploads.
        if export_settings.include_raw_files:
            yield DiskFileSource(self.os_path, 'raw')
        if export_settings.include_archive_files:
            yield DiskFileSource(self.os_path, 'archive')

    @classmethod
    def files_from_bundle(
            cls, bundle_file_source: BrowsableFileSource, import_settings: BundleImportSettings) -> Iterable[FileSource]:
        # Files to import for a staging upload
        if import_settings.include_raw_files:
            yield bundle_file_source.sub_source('raw')
        if import_settings.include_archive_files:
            yield bundle_file_source.sub_source('archive')
        if import_settings.include_bundle_info:
            yield bundle_file_source.sub_source(bundle_info_filename)


class PublicUploadFiles(UploadFiles):

    def __init__(self, upload_id: str, create: bool = False):
        super().__init__(upload_id, create)
        self._directories: Dict[str, Dict[str, RawPathInfo]] = None
        self._raw_zip_file_object: PathObject = None
        self._raw_zip_file: zipfile.ZipFile = None
        self._archive_msg_file_object: PathObject = None
        self._archive_msg_file: ArchiveReader = None
        self._access: str = None
        self._missing_raw_files: bool = None

    @classmethod
    def file_area(cls):
        return config.fs.public

    @property
    def external_os_path(self):
        if not config.fs.public_external:
            return self.os_path

        self.os_path.replace(config.fs.public, config.fs.public_external)

    def close(self):
        if self._raw_zip_file is not None:
            self._raw_zip_file.close()

        if self._archive_msg_file is not None:
            self._archive_msg_file.close()

    @property
    def access(self):
        '''
        Which "access" is used, either 'public' (uploads without embargo) or 'restricted'
        (uploads with embargo). This is reflected in the names of the files holding the
        raw data and the archive data. The reason for this is so that it should be easy to
        see, by just looking at the files, if a published upload is embargoed or not.

        The access is determined by inspecting which files exist/contain data. If both
        public and restricted files exist/contain data, or if neither exists/contain data,
        a KeyError will be thrown (this should not happen if the upload is correctly packed).
        The inspection of the files is only done on the first call, and the cached result
        is used in subsequent calls. The only way to change the access is to call :func:`re_pack`.
        '''
        if self._access:
            return self._access
        # Determine access by inspecting the files
        files_found = False
        for access in ('public', 'restricted'):
            raw_zip_file_object = PublicUploadFiles._create_raw_zip_file_object(self, access)
            archive_msg_file_object = PublicUploadFiles._create_msg_file_object(self, access)
            found = (raw_zip_file_object.exists() and raw_zip_file_object.size > empty_zip_file_size) or (
                archive_msg_file_object.exists() and archive_msg_file_object.size > empty_archive_file_size)
            if found:
                if files_found:
                    self._access = self._raw_zip_file_object = self._archive_msg_file_object = None
                    raise KeyError('Inconsistency: both public and restricted files found')
                files_found = True
                self._raw_zip_file_object = raw_zip_file_object
                self._archive_msg_file_object = archive_msg_file_object
                self._access = access
        if not files_found:
            raise KeyError('Neither public nor restricted files found')
        return self._access

    @staticmethod
    def _create_raw_zip_file_object(target_dir: DirectoryObject, access: str) -> PathObject:
        return target_dir.join_file(f'raw-{access}.plain.zip')

    def raw_zip_file_object(self) -> PathObject:
        '''
        Gets the raw zip file, either public or restricted, depending on which one is used.
        If both public and restricted files exist, or if none of them exist, a KeyError will
        be thrown.
        '''
        self.access  # Invoke to initialize
        return self._raw_zip_file_object

    def _open_raw_zip_file(self) -> zipfile.ZipFile:
        if self._raw_zip_file:
            return self._raw_zip_file

        zip_path = self.raw_zip_file_object().os_path
        self._raw_zip_file = zipfile.ZipFile(zip_path)

        return self._raw_zip_file

    @property
    def missing_raw_files(self):
        if self._missing_raw_files is None:
            self._missing_raw_files = not os.path.exists(self.raw_zip_file_object().os_path)
        return self._missing_raw_files

    @staticmethod
    def _create_msg_file_object(target_dir: DirectoryObject, access: str) -> PathObject:
        version_suffix = '-' + config.fs.archive_version_suffix if config.fs.archive_version_suffix else ''
        return target_dir.join_file(f'archive-{access}{version_suffix}.msg.msg')

    def msg_file_object(self) -> PathObject:
        '''
        Gets the msg file, either public or restricted, depending on which one is used.
        If both public and restricted files exist, or if none of them exist, a KeyError will
        be thrown.
        '''
        self.access  # Invoke to initialize
        return self._archive_msg_file_object

    def _open_msg_file(self, use_blocked_toc: bool = True) -> ArchiveReader:
        if self._archive_msg_file is not None:
            if not self._archive_msg_file.is_closed():
                return self._archive_msg_file

        msg_file_object = self.msg_file_object()

        if not msg_file_object.exists():
            raise FileNotFoundError()

        archive = read_archive(msg_file_object.os_path, use_blocked_toc=use_blocked_toc)
        assert archive is not None
        self._archive_msg_file = archive

        return archive

    def to_staging_upload_files(self, create: bool = False, include_archive: bool = False) -> 'StagingUploadFiles':
        exists = StagingUploadFiles.exists_for(self.upload_id)
        if exists:
            if create:
                raise FileExistsError('Staging upload does already exist')
            return StagingUploadFiles(self.upload_id)
        elif not create:
            return None

        staging_upload_files = StagingUploadFiles(self.upload_id, create=True)
        # Extract files
        raw_zip_file = self.raw_zip_file_object()
        if raw_zip_file.exists():
            staging_upload_files.add_rawfiles(raw_zip_file.os_path)

        if include_archive:
            with self._open_msg_file() as archive:
                for entry_id, data in archive.items():
                    entry_id = entry_id.strip()
                    staging_upload_files.write_archive(entry_id, data.to_dict())

        return staging_upload_files

    def _parse_content(self):
        '''
        Parses the content of files and folders and caches it in self._directories for
        faster future access.
        '''
        if self._directories is None:
            self._directories = dict()
            self._directories[''] = {}  # Root folder
            directory_sizes: Dict[str, int] = {}
            # Add file RawPathInfo objects and calculate directory sizes
            try:
                zf = self._open_raw_zip_file()
                for path in zf.namelist():
                    file_name = os.path.basename(path)
                    directory_path = os.path.dirname(path)
                    size = zf.getinfo(path).file_size if file_name else 0

                    if directory_path:
                        # Ensure that all parent directories are added
                        sub_path = ''
                        for directory in directory_path.split(os.path.sep):
                            sub_path_next = os.path.join(sub_path, directory)
                            if sub_path_next not in self._directories:
                                self._directories[sub_path_next] = {}
                                directory_sizes.setdefault(sub_path_next, 0)
                            directory_sizes[sub_path_next] += size
                            sub_path = sub_path_next

                    if file_name:
                        directory_content = self._directories[directory_path]
                        directory_content[file_name] = RawPathInfo(
                            path=path,
                            is_file=True,
                            size=size,
                            access=self.access)
                        self._directories[path] = directory_content[file_name]
            except FileNotFoundError:
                pass
            # Add directories with the calculated sizes.
            for path, size in directory_sizes.items():
                basename = os.path.basename(path)
                directory_path = os.path.dirname(path)
                self._directories[directory_path][basename] = RawPathInfo(
                    path=path, is_file=False, size=size, access=self.access)

    def is_empty(self) -> bool:
        self._parse_content()
        return not self._directories.get('')

    def raw_path_exists(self, path: str) -> bool:
        if not is_safe_relative_path(path):
            return False
        if self.missing_raw_files:
            return not path  # We consider the empty path (i.e. root) to always "exists".
        self._parse_content()
        explicit_directory_path = path.endswith(os.path.sep)
        path = path.rstrip(os.path.sep)
        base_name = os.path.basename(path)
        directory_path = os.path.dirname(path)
        directory_content = self._directories.get(directory_path)
        if directory_content is not None:
            if not base_name:
                return True
            if base_name in directory_content:
                path_info = directory_content[base_name]
                if explicit_directory_path and path_info.is_file:
                    return False
                return True
        return False

    def raw_path_is_file(self, path: str) -> bool:
        if not is_safe_relative_path(path) or self.missing_raw_files:
            return False
        self._parse_content()
        base_name = os.path.basename(path)
        directory_path = os.path.dirname(path)
        if not base_name:
            return False  # Requested path is an explicit directory path
        directory_content = self._directories.get(directory_path)
        if directory_content and base_name in directory_content:
            path_info = directory_content[base_name]
            return path_info.is_file
        return False

    def raw_directory_list(
            self, path: str = '', recursive=False, files_only=False, path_prefix=None) -> Iterable[RawPathInfo]:
        if not is_safe_relative_path(path):
            return
        if not path and self.missing_raw_files:
            return
        self._parse_content()
        path = path.rstrip(os.path.sep)
        directory_content = self._directories.get(path)
        if directory_content is not None:
            for __, path_info in sorted(directory_content.items()):
                if not files_only or path_info.is_file:
                    if not path_prefix or path_info.path.startswith(path_prefix):
                        yield path_info
                if recursive and not path_info.is_file:
                    for sub_path_info in self.raw_directory_list(path_info.path, recursive, files_only):
                        if not path_prefix or sub_path_info.path.startswith(path_prefix):
                            yield sub_path_info

    def raw_file(self, file_path: str, *args, **kwargs) -> IO:
        assert is_safe_relative_path(file_path)
        mode = kwargs.get('mode') if len(args) == 0 else args[0]
        if 'mode' in kwargs:
            del(kwargs['mode'])
        mode = mode if mode else 'rb'

        try:
            zf = self._open_raw_zip_file()
            f = zf.open(file_path, 'r', **kwargs)
            if 't' in mode:
                return io.TextIOWrapper(f)
            else:
                return f
        except FileNotFoundError:
            pass
        except IsADirectoryError:
            pass
        except KeyError:
            pass

        raise KeyError(file_path)

    def raw_file_size(self, file_path: str) -> int:
        assert is_safe_relative_path(file_path)
        try:
            zf = self._open_raw_zip_file()
            info = zf.getinfo(file_path)
            return info.file_size
        except FileNotFoundError:
            pass
        except KeyError:
            pass

        raise KeyError(file_path)

    def read_archive(self, entry_id: str, use_blocked_toc: bool = True) -> Any:
        try:
            archive = self._open_msg_file(use_blocked_toc)
            if entry_id in archive:
                return archive
        except FileNotFoundError:
            pass

        raise KeyError(entry_id)

    def re_pack(self, with_embargo: bool) -> None:
        '''
        Repacks the files when changing the embargo flag on the upload. That is: when lifting the
        embargo the file names of the raw zip file and the archive file change from containing
        the keyword "restricted" to "public". Adding embargo to a non-embargoed published
        upload is also supported, but only admins should be allowed to do this. The existing
        files are just renamed, so this should be a rather quick operation. The upload must
        be correctly packed (i.e. there cannot be non-empty public and restricted files
        at the same time).
        '''
        if (self.access == 'restricted' and with_embargo) or (self.access == 'public' and not with_embargo):
            return  # Nothing to do
        self.close()
        new_access = 'restricted' if with_embargo else 'public'
        msg_file_object = self.msg_file_object()
        msg_file_object_new = PublicUploadFiles._create_msg_file_object(self, new_access)
        if msg_file_object.exists():
            if msg_file_object_new.exists():
                msg_file_object_new.delete()  # We have checked that the file is empty anyway
            os.rename(msg_file_object.os_path, msg_file_object_new.os_path)
        raw_zip_file_object = self.raw_zip_file_object()
        raw_zip_file_object_new = PublicUploadFiles._create_raw_zip_file_object(self, new_access)
        if raw_zip_file_object.exists():
            if raw_zip_file_object_new.exists():
                raw_zip_file_object_new.delete()  # We have checked that the file is empty anyway
            os.rename(raw_zip_file_object.os_path, raw_zip_file_object_new.os_path)

        # Clear the cached values
        self._access = None
        self._raw_zip_file = self._raw_zip_file_object = None
        self._archive_msg_file = self._archive_msg_file_object = None

    def files_to_bundle(self, export_settings: BundleExportSettings) -> Iterable[FileSource]:
        # Defines files for upload bundles of published uploads.
        for filename in sorted(os.listdir(self.os_path)):
            if filename.startswith('raw-') and export_settings.include_raw_files:
                yield DiskFileSource(self.os_path, filename)
            if filename.startswith('archive-') and export_settings.include_archive_files:
                yield DiskFileSource(self.os_path, filename)

    @classmethod
    def files_from_bundle(
            cls, bundle_file_source: BrowsableFileSource, import_settings: BundleImportSettings) -> Iterable[FileSource]:
        for filename in bundle_file_source.directory_list(''):
            if filename.startswith('raw-') and import_settings.include_raw_files:
                yield bundle_file_source.sub_source(filename)
            if filename.startswith('archive-') and import_settings.include_archive_files:
                yield bundle_file_source.sub_source(filename)
            if filename == bundle_info_filename and import_settings.include_bundle_info:
                yield bundle_file_source.sub_source(filename)
