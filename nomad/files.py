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

"""
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
"""

from __future__ import annotations

import hashlib
import io
import json
import os
import shutil
import tempfile
import zipfile
from abc import ABC, abstractmethod
from collections.abc import Callable, Iterable, Iterator
from contextlib import contextmanager, suppress
from datetime import datetime
from functools import cached_property
from pathlib import Path
from typing import IO, Any, Literal, NamedTuple

import magic
import yaml
import zipstream
from deprecation import deprecated
from fsspec import AbstractFileSystem
from fsspec.implementations.local import LocalFileSystem
from fsspec.implementations.tar import TarFileSystem
from fsspec.implementations.zip import ZipFileSystem
from pathvalidate import sanitize_filename, sanitize_filepath
from pydantic import BaseModel

from nomad import datamodel, utils
from nomad.archive import (
    ArchiveReader,
    combine_archive,
    read_archive,
    to_json,
    write_archive,
)
from nomad.common import get_compression_format, is_safe_relative_path
from nomad.config import config
from nomad.config.models.config import BundleExportSettings, BundleImportSettings

bundle_info_filename = 'bundle_info.json'

empty_zip_file_size = 22
empty_archive_file_size = 32
empty_hdf5_file_size = 96

UPath = Path


def mkdtemp(prefix: str):
    return tempfile.mkdtemp(None, sanitize_filename(prefix), config.fs.tmp)


class PathObject:
    """
    Object storage-like abstraction for paths in general.
    Attributes:
        os_path: The full os path of the object.
    """

    def __init__(self, os_path: str, *, fs: AbstractFileSystem | None = None):
        self.os_path = os_path
        self._fs = fs or LocalFileSystem()

    def delete(self):
        if self.exists():
            self._fs.rm(self.os_path, recursive=True)

    def exists(self):
        return self._fs.exists(self.os_path)

    @property
    def size(self):
        return self._fs.size(self.os_path)

    def __repr__(self) -> str:
        return self.os_path


class DirectoryObject(PathObject):
    """
    Object storage-like abstraction for directories.
    """

    def __init__(
        self,
        os_path: str,
        create: bool = False,
        *,
        fs: AbstractFileSystem | None = None,
    ):
        super().__init__(os_path, fs=fs)
        if create:
            self._fs.mkdirs(self.os_path, exist_ok=True)

    def join_dir(self, path, create: bool = False) -> DirectoryObject:
        return DirectoryObject(os.path.join(self.os_path, path), create, fs=self._fs)

    def join_file(self, path, create_dir: bool = False) -> PathObject:
        target_path = os.path.join(self.os_path, path)

        if create_dir and (target_folder := os.path.dirname(target_path)):
            self._fs.mkdirs(target_folder, exist_ok=True)

        return PathObject(target_path, fs=self._fs)

    def exists(self) -> bool:
        return self._fs.isdir(self.os_path)

    def zip_fp(self, access: str):
        """
        Utilities
        """
        return self.join_file(f'raw-{access}.plain.zip')

    def msg_fp(self, access: str, fallback: bool = False):
        def versioned_file_name(version_suffix):
            return f'archive-{access}{version_suffix}.msg.msg'

        return _versioned_archive_file_object(self, versioned_file_name, fallback)

    def h5_fp(self, access: str):
        return self.join_file(f'archive-{access}.h5')


class RawPathInfo(NamedTuple):
    """
    Stores basic info about a file or folder located at a specific raw path.
    """

    path: str
    is_file: bool
    size: int
    access: str


class RawDirPage(NamedTuple):
    """
    A paginated slice of raw directory metadata.
    """

    content: list[RawPathInfo]
    total: int
    total_size: int | None = None


class StreamedFile(BaseModel):
    """
    Convenience class for representing a streamed file, together with information about
    file size and an associated path.
    """

    src: Any = None
    path: str
    size: int


class FileSource(ABC):
    """
    An abstract class which represents a generic "file source", from which some number of files
    can be retrieved. There are several different ways to create a file source, see subclasses.
    The files in the source are associated with paths and have known sizes.
    """

    def __init__(self, fs: AbstractFileSystem | None = None):
        self._fs = fs or LocalFileSystem()

    @abstractmethod
    def to_streamed_files(self) -> Iterable[StreamedFile]:
        """
        Retrieves the files in the source as :class:`StreamedFile` objects.
        The caller should close the streams when consumed.
        """
        ...

    def to_zipfile(self, path, overwrite: bool = False):
        """
        Generates a zip file from the files in this FileSource and stores it to disk. The
        zipfile content is created by calling :func:`to_zipstream`.
        """
        assert not self._fs.isdir(path), (
            'Exporting to zip file requires a file path, not directory.'
        )
        assert overwrite or not self._fs.exists(path), (
            '`path` already exists. Use `overwrite` to overwrite.'
        )
        with self._fs.open(path, 'wb') as f:
            for chunk in create_zipstream(self.to_streamed_files()):
                f.write(chunk)

    def to_disk(
        self, destination_dir: str, move_files: bool = False, overwrite: bool = False
    ):
        """
        Writes the files from this FileSource to disk, uncompressed. The default implementation
        makes use of :func:`to_streamed_files`. The `destination_dir` should be a directory
        (it will be created if it does not exist). The `move_files` argument instructs
        the method to move the source files if possible.
        """
        dest_path = UPath(destination_dir)
        self._fs.mkdirs(dest_path, exist_ok=True)
        for streamed_file in self.to_streamed_files():
            full_path = dest_path / streamed_file.path
            if full_path.exists():
                assert overwrite, 'Target already exists and `overwrite` is False'
            self._fs.mkdirs(full_path.parent, exist_ok=True)
            with (
                self._fs.open(full_path.as_posix(), 'wb') as output_file,
                streamed_file.src,
            ):
                shutil.copyfileobj(streamed_file.src, output_file)

    def close(self):
        """Perform "closing" of the source, if applicable."""
        pass


class BrowsableFileSource(FileSource, ABC):
    """
    A :class:`FileSource` which can be "browsed", like a folder on disk or a zip archive.
    """

    @abstractmethod
    def open(self, path, mode='rb') -> IO:
        """Opens a file by the specified path."""
        ...

    @abstractmethod
    def find(self, path: str) -> list[str]:
        """
        Returns a list of directory contents, located in the directory denoted by `path`
        in this file source.
        """
        ...

    @abstractmethod
    def child(self, path: str) -> BrowsableFileSource:
        """
        Creates a new instance of :class:`BrowsableFileSource` which just contains the
        files located under the specified path.
        """
        ...


class StreamedFileSource(FileSource):
    """
    A :class:`FileSource` created from a single :class:`StreamedFile`.
    """

    def __init__(
        self, streamed_file: StreamedFile, fs: AbstractFileSystem | None = None
    ):
        super().__init__(fs)
        self._file = streamed_file

    def to_streamed_files(self) -> Iterable[StreamedFile]:
        yield self._file


class DiskFileSource(BrowsableFileSource):
    """
    A :class:`FileSource` corresponding to a single file or a folder on disk. The object
    is identified by a `base_path` and a `relative path`. The `base_path` should be a folder,
    the `relative_path` is optional, and used for selecting only a specific file or folder
    located under `base_folder`. The paths of the files retrieved from this source are given
    relative to the `base_path`.
    """

    def __init__(
        self,
        base_path: str,
        relative_path: str | None = None,
        fs: AbstractFileSystem | None = None,
    ):
        super().__init__(fs)
        assert self._fs.isdir(base_path)
        if relative_path:
            relative_path = sanitize_filepath(relative_path)
            assert is_safe_relative_path(relative_path), 'Unsafe relative_path received'
            self.full_path = os.path.join(base_path, relative_path)
            assert self._fs.exists(self.full_path)
        else:
            self.full_path = base_path
        self.base_path = base_path
        self.relative_path = relative_path

    def to_streamed_files(self) -> Iterable[StreamedFile]:
        for target_path in self._fs.find(self.full_path):
            yield StreamedFile(
                path=os.path.relpath(target_path, self.base_path),
                src=self._fs.open(target_path, 'rb'),
                size=self._fs.size(target_path),
            )

    def to_disk(
        self, destination_dir: str, move_files: bool = False, overwrite: bool = False
    ):
        destination_path = UPath(destination_dir)
        if self.relative_path:
            destination_path /= self.relative_path

        self._fs.mkdirs(destination_path.parent, exist_ok=True)

        if self._fs.exists(destination_path):
            assert overwrite, (
                f'Target {destination_path} already exists and `overwrite` is False'
            )

        self._fs.put(self.full_path, destination_path, recursive=True)

        if move_files:
            self._fs.rm(self.full_path, recursive=True)

    def open(self, path, mode='rb') -> IO:
        assert is_safe_relative_path(path)
        return self._fs.open(os.path.join(self.base_path, path), mode)

    def find(self, path: str) -> list[str]:
        assert is_safe_relative_path(path)
        return self._fs.find(os.path.join(self.base_path, path))

    def child(self, path: str) -> DiskFileSource:
        assert is_safe_relative_path(path)
        return DiskFileSource(self.base_path, path)


class ZipFileSource(BrowsableFileSource):
    """
    Allows us to "wrap" a :class:`zipfile.ZipFile` object and use it as a :class:`BrowsableFileSource`,
    i.e. it denotes a resource (single file or folder) stored in a ZipFile.
    """

    def __init__(
        self,
        zip_file: str,
        sub_path: str = '',
        fs: AbstractFileSystem | None = None,
    ):
        super().__init__(fs)
        assert is_safe_relative_path(sub_path)
        self.sub_path = sub_path
        self._zip_fs = ZipFileSystem(zip_file)

    def to_streamed_files(self) -> Iterable[StreamedFile]:
        for target_path in self._zip_fs.find(self.sub_path):
            yield StreamedFile(
                path=target_path,
                src=self._zip_fs.open(target_path),
                size=self._zip_fs.size(target_path),
            )

    def open(self, path, mode='rb') -> IO:
        assert 'r' in mode, 'Mode must be a read mode'
        assert all(c in 'rbt' for c in mode), f'Invalid mode for open command: {mode}'
        f = self._zip_fs.open(path)
        return io.TextIOWrapper(f) if 't' in mode else f

    def find(self, path: str) -> list[str]:
        return self._zip_fs.find(path)

    def child(self, path: str) -> ZipFileSource:
        assert is_safe_relative_path(path), 'Unsafe path provided'
        if self.sub_path:
            assert path.startswith(self.sub_path + os.path.sep), (
                'Provided `path` is not a sub path.'
            )
        return ZipFileSource(self._zip_fs.fo, path)

    def close(self):
        self._zip_fs.close()


class CombinedFileSource(FileSource):
    """
    Class for defining a :class:`FileSource` by combining multiple "subsources" into one.
    """

    def __init__(
        self, file_sources: Iterable[FileSource], fs: AbstractFileSystem | None = None
    ):
        """file_sources: an Iterable for getting FileSources."""
        super().__init__(fs)
        self._files = file_sources

    def to_streamed_files(self) -> Iterable[StreamedFile]:
        for file in self._files:
            yield from file.to_streamed_files()

    def to_disk(
        self, destination_dir: str, move_files: bool = False, overwrite: bool = False
    ):
        for file in self._files:
            file.to_disk(destination_dir, move_files, overwrite)


class StandardJSONDecoder(json.JSONDecoder):
    """Our standard JSONDecoder, with support for marshaling of datetime objects"""

    def __init__(self, *args, **kwargs):
        def dict_to_object(d: dict):
            if len(d) == 1 and (v := d.get('$datetime')) is not None:
                return datetime.fromisoformat(v)
            return d

        kwargs['object_hook'] = dict_to_object
        super().__init__(**kwargs)


def json_to_streamed_file(json_dict: dict[str, Any], path: str) -> StreamedFile:
    """Converts a json dictionary structure to a :class:`StreamedFile`."""

    class StandardJSONEncoder(json.JSONEncoder):
        """Our standard JSONEncoder with support for marshaling of datetime objects"""

        def default(self, obj):
            if isinstance(obj, datetime):
                return {'$datetime': obj.isoformat()}
            return super().default(obj)

    json_bytes = json.dumps(json_dict, cls=StandardJSONEncoder).encode()
    return StreamedFile(path=path, src=io.BytesIO(json_bytes), size=len(json_bytes))


def create_zipstream(streamed_files: Iterable[StreamedFile], compress: bool = False):
    """
    Creates a zip stream, i.e. a streamed zip file.
    """
    zs = zipstream.ZipStream(
        compress_type=zipfile.ZIP_DEFLATED if compress else zipfile.ZIP_STORED,
        compress_level=9,
    )

    def content_generator(file):
        with file.src as f:
            while data := f.read(1024 * 1024):
                yield data

    for streamed_file in streamed_files:
        zs.add(content_generator(streamed_file), streamed_file.path)

    yield from zs


async def create_zipstream_async(
    streamed_files: Iterable[StreamedFile], compress: bool = False
):
    for x in create_zipstream(streamed_files, compress):
        yield x


def _versioned_archive_file_object(
    target_dir: DirectoryObject, file_name: Callable[[str], str], fallback: bool
) -> PathObject:
    """
    Creates a file object for an archive file depending on the directory it is or
    will be created in, the recipe to construct the name from a version suffix, and
    a bool that denotes if alternative version suffixes should be considered.
    """
    suffixes = config.fs.archive_version_suffix

    if not isinstance(suffixes, list):
        suffixes = [suffixes]

    if len(suffixes) <= 1:
        return target_dir.join_file(file_name(f'-{suffixes[0]}' if suffixes[0] else ''))

    if not fallback:
        return target_dir.join_file(file_name(f'-{suffixes[0]}'))

    for suffix in suffixes:
        current_file = target_dir.join_file(file_name(f'-{suffix}'))
        if os.path.exists(current_file.os_path):
            return current_file

    return target_dir.join_file(file_name(f'-{suffixes[0]}'))


class UploadFiles(DirectoryObject):
    """Abstract base class for upload files."""

    def __init__(
        self,
        upload_id: str,
        create: bool = False,
        *,
        fs: AbstractFileSystem | None = None,
    ):
        self.logger = utils.get_logger(__name__, upload_id=upload_id)

        super().__init__(os_path=self.base_folder_for(upload_id), create=create, fs=fs)

        if not create and not self.exists():
            raise KeyError(upload_id)

        self.upload_id = upload_id

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    @classmethod
    def _file_area(cls) -> UPath:
        """
        Full path to where the upload files of this class are stored (i.e. either
        staging or public file area).
        """
        raise NotImplementedError()

    @property
    def external_os_path(self):
        """
        Full path to where the upload files of this class are stored on the server.
        This is equal to `self.os_path` if no external path substitutes for staging
        and public area are configured. This is helpful, when nomad is run in a container
        and the mounted path used by nomad are different from the actual paths on the
        host server.
        """
        raise NotImplementedError()

    @classmethod
    def base_folder_for(cls, upload_id: str) -> str:
        """
        Full path to the base folder for the upload files (of this class) for the
        specified upload_id.
        """
        full_path = cls._file_area() / upload_id[: config.fs.prefix_size] / upload_id
        return full_path.as_posix()

    @classmethod
    def exists_for(cls, upload_id: str) -> bool:
        """
        If an UploadFiles object (of this class) has been created for this upload_id.
        """
        return os.path.exists(cls.base_folder_for(upload_id))

    def to_staging(
        self, create: bool = False, include_archive: bool = False
    ) -> StagingUploadFiles | None:
        """Casts to or creates corresponding staging upload files or returns None."""
        raise NotImplementedError()

    @staticmethod
    def get(upload_id: str, create: bool = False) -> UploadFiles | None:
        for class_type in (PublicUploadFiles, StagingUploadFiles):
            if class_type.exists_for(upload_id):
                return class_type(upload_id, create)

        return None

    def is_empty(self) -> bool:
        """If this upload has no content yet."""
        raise NotImplementedError()

    def raw_exists(self, path: str) -> bool:
        """
        Returns True if the specified path is a valid raw path (either file or directory)
        """
        raise NotImplementedError()

    def raw_isfile(self, path: str) -> bool:
        """
        Returns True if the specified path points to a file (rather than a directory).
        """
        raise NotImplementedError()

    def raw_listdir(
        self,
        path: str = '',
        recursive=False,
        files_only=False,
        depth: int = -1,
    ) -> Iterable[RawPathInfo]:
        """
        Returns an iterable of RawPathInfo, one for each element (file or folder) in
        the directory specified by `path`. If `recursive` is set to True, subdirectories are
        also crawled. If `files_only` is set, only the file objects found are returned.
        If path is not a valid directory, the result will be empty. Selecting empty string
        as path (which is the default value) gives the content of the whole raw directory.
        The `path_prefix` argument can be used to filter out elements where the path starts
        with a specific prefix.

        The `depth` argument can be used to limit the depth of the recursion.
        """
        raise NotImplementedError()

    def raw_listdir_page(
        self,
        path: str = '',
        *,
        start: int = 0,
        end: int | None = None,
        recursive=False,
        files_only=False,
        depth: int = -1,
        order: Literal['asc', 'desc'] = 'asc',
        group_directories_first: bool = False,
        include_total_size: bool = False,
    ) -> RawDirPage:
        """
        Returns a paginated slice of raw directory metadata.

        Subclasses can override this to avoid materializing full metadata for items outside
        the requested page.
        """
        items = list(self.raw_listdir(path, recursive, files_only, depth))
        total_size = sum(item.size for item in items) if include_total_size else None

        if group_directories_first:
            folders = [item for item in items if not item.is_file]
            files = [item for item in items if item.is_file]
            ordered = folders + files
            if order != 'asc':
                ordered = list(reversed(files)) + list(reversed(folders))
        else:
            ordered = items if order == 'asc' else list(reversed(items))

        if end is None:
            end = len(ordered)

        return RawDirPage(
            content=ordered[start:end], total=len(items), total_size=total_size
        )

    @deprecated(details='Use raw_exists() instead.')
    def raw_path_exists(self, path: str) -> bool:
        return self.raw_exists(path)

    @deprecated(details='Use raw_isfile() instead.')
    def raw_path_is_file(self, path: str) -> bool:
        return self.raw_isfile(path)

    @deprecated(details='Use raw_listdir() instead.')
    def raw_directory_list(
        self,
        path: str = '',
        recursive=False,
        files_only=False,
        depth: int = -1,
    ) -> Iterable[RawPathInfo]:
        return self.raw_listdir(path, recursive, files_only, depth)

    @contextmanager
    def raw_file(self, file_path: str, *args, **kwargs):
        """
        Opens a raw file and returns a file-like object. Additional args, kwargs are
        delegated to the respective `open` call.
        Arguments:
            file_path: The path to the file relative to the upload.
        Raises:
            KeyError: If the file does not exist.
        """
        raise NotImplementedError()

    def raw_file_size(self, file_path: str) -> int:
        """
        Returns:
            The size of the given raw file.
        """
        raise NotImplementedError()

    def raw_file_mime_type(self, file_path: str) -> str:
        assert self.raw_isfile(file_path), (
            'Provided path does not specify a file, or is invalid.'
        )
        with self.raw_file(file_path, 'br') as raw_file:
            return (
                magic.from_buffer(raw_file.read(2048), mime=True)
                or 'application/octet-stream'
            )

    @contextmanager
    def read_archive(self, entry_id: str) -> Iterator[ArchiveReader]:
        """
        Returns an :class:`nomad.archive.ArchiveReader` that contains the
        given entry_id.
        """
        raise NotImplementedError()

    def close(self):
        """Release possibly held system resources (e.g. file handles)."""
        pass

    def delete(self) -> None:
        super().delete()
        if config.fs.prefix_size > 0 and not self._fs.ls(
            parent := os.path.dirname(self.os_path), False
        ):
            self._fs.rm(parent, recursive=True)

    def files_to_bundle(
        self, export_settings: BundleExportSettings
    ) -> Iterable[FileSource]:
        """
        A generator of :class:`FileSource` objects, defining the files/folders to be included in an
        upload bundle when *exporting*. The arguments allows for further filtering of what to include.

        Note, this only yields files to copy from the regular upload directory, not "special" files,
        like the bundle_info.json file, which is created by the :class:`BundleExporter`.
        """
        raise NotImplementedError()

    @classmethod
    def files_from_bundle(
        cls,
        bundle_file_source: BrowsableFileSource,
        import_settings: BundleImportSettings,
    ) -> Iterable[FileSource]:
        """
        Returns an Iterable of :class:`FileSource`, defining the files/folders to be included in an
        upload bundle when *importing*. Only the files specified by the import_settings are included.
        """
        raise NotImplementedError()

    def archive_hdf5_location(self, entry_id: str) -> str:
        """
        Returns the OS path to the target HDF5 file.
        The str will be passed to h5py module for reading and writing.
        We do not provide a raw IO object here, since later this file may be a web resource.
        """
        raise NotImplementedError()


class StagingUploadFiles(UploadFiles):
    def __init__(self, upload_id: str, create: bool = False):
        super().__init__(upload_id, create)

        self._raw_dir = self.join_dir('raw', create)
        self._archive_dir = self.join_dir('archive', create)

        self._size = 0

    @classmethod
    def _file_area(cls):
        return UPath(config.fs.staging)

    @property
    def _frozen_file(self):
        return self.join_file('.frozen')

    @property
    def external_os_path(self):
        if not config.fs.staging_external:
            return self.os_path

        return self.os_path.replace(config.fs.staging, config.fs.staging_external)

    def to_staging(
        self, create: bool = False, include_archive: bool = False
    ) -> StagingUploadFiles | None:
        return self

    @property
    def size(self) -> int:
        return self._fs.du(self._raw_dir.os_path)

    def _full_path(self, path: str):
        return UPath(self._raw_dir.os_path) / path

    def is_empty(self) -> bool:
        return not self._fs.ls(self._raw_dir.os_path, False)

    def raw_exists(self, path: str) -> bool:
        if not is_safe_relative_path(path):
            return False
        return self._fs.exists(self._full_path(path))

    def raw_isfile(self, path: str) -> bool:
        if not is_safe_relative_path(path):
            return False
        return self._fs.isfile(self._full_path(path))

    def raw_create_directory(self, path: str):
        assert path and is_safe_relative_path(path), 'Bad path provided'
        self._fs.makedirs(self._full_path(path).as_posix(), True)

    def raw_listdir(
        self,
        path: str = '',
        recursive=False,
        files_only=False,
        depth: int = -1,
    ) -> Iterable[RawPathInfo]:
        if not is_safe_relative_path(path) or depth == 0:
            return

        fs = self._fs
        for target in fs.find(
            os.path.join(self._raw_dir.os_path, path),
            (depth if depth > 0 else None) if recursive else 1,
            not files_only,
        ):
            relpath = UPath(os.path.relpath(target, self._raw_dir.os_path))
            if not (isfile := fs.isfile(target)) and relpath == UPath(path):
                # skip folder itself
                continue
            yield RawPathInfo(
                path=relpath.as_posix(),
                is_file=isfile,
                size=fs.size(target) if isfile else fs.du(target),
                access='unpublished',
            )

    def raw_listdir_page(
        self,
        path: str = '',
        *,
        start: int = 0,
        end: int | None = None,
        recursive=False,
        files_only=False,
        depth: int = -1,
        order: Literal['asc', 'desc'] = 'asc',
        group_directories_first: bool = False,
        include_total_size: bool = False,
    ) -> RawDirPage:
        """List a raw directory page using local filesystem primitives for staging uploads."""
        if not isinstance(self._fs, LocalFileSystem):
            return super().raw_listdir_page(
                path,
                start=start,
                end=end,
                recursive=recursive,
                files_only=files_only,
                depth=depth,
                order=order,
                group_directories_first=group_directories_first,
                include_total_size=include_total_size,
            )

        if not is_safe_relative_path(path) or depth == 0:
            return RawDirPage(
                content=[], total=0, total_size=0 if include_total_size else None
            )

        os_path = self._full_path(path).as_posix()
        if not self._fs.exists(os_path):
            return RawDirPage(
                content=[], total=0, total_size=0 if include_total_size else None
            )

        items: list[tuple[str, bool]] = []
        total_size = 0 if include_total_size else None
        normalised_path = path.rstrip('/')

        if self._fs.isfile(os_path):
            items.append((normalised_path, True))
            if include_total_size:
                total_size = os.stat(os_path).st_size
        else:
            total_size = self._collect_local_raw_paths(
                os_path,
                normalised_path,
                items,
                recursive=recursive,
                files_only=files_only,
                depth=depth,
                include_total_size=include_total_size,
            )

        if group_directories_first:
            folders = sorted(
                (item for item in items if not item[1]), key=lambda item: item[0]
            )
            files = sorted(
                (item for item in items if item[1]), key=lambda item: item[0]
            )
            ordered = folders + files
            if order != 'asc':
                ordered = list(reversed(files)) + list(reversed(folders))
        else:
            ordered = sorted(items, key=lambda item: item[0], reverse=order != 'asc')

        if end is None:
            end = len(ordered)

        content = [
            self._raw_path_info(relative_path, is_file)
            for relative_path, is_file in ordered[start:end]
        ]

        return RawDirPage(content=content, total=len(items), total_size=total_size)

    def _collect_local_raw_paths(
        self,
        os_path: str,
        relative_path: str,
        items: list[tuple[str, bool]],
        *,
        recursive: bool,
        files_only: bool,
        depth: int,
        include_total_size: bool,
    ) -> int | None:
        """Collect relative raw paths for a local directory walk and optionally accumulate size."""
        total_size = 0 if include_total_size else None
        remaining_depth = (
            depth if recursive and depth > 0 else (None if recursive else 1)
        )

        def _walk(current_os_path: str, current_relative_path: str, current_depth):
            nonlocal total_size

            with os.scandir(current_os_path) as iterator:
                children = sorted(iterator, key=lambda entry: entry.name)

            for child in children:
                child_relative_path = (
                    child.name
                    if not current_relative_path
                    else f'{current_relative_path}/{child.name}'
                )
                is_file = child.is_file(follow_symlinks=False)

                if is_file:
                    items.append((child_relative_path, True))
                    if include_total_size:
                        total_size += child.stat(follow_symlinks=False).st_size
                    continue

                if not files_only:
                    items.append((child_relative_path, False))
                    if include_total_size:
                        total_size += self._fs.du(child.path)

                if current_depth is None or current_depth > 1:
                    _walk(
                        child.path,
                        child_relative_path,
                        None if current_depth is None else current_depth - 1,
                    )

        _walk(os_path, relative_path, remaining_depth)
        return total_size

    def _raw_path_info(self, relative_path: str, is_file: bool) -> RawPathInfo:
        """Build a RawPathInfo for a known local raw path."""
        os_path = self._full_path(relative_path).as_posix()
        size = os.stat(os_path).st_size if is_file else self._fs.du(os_path)
        return RawPathInfo(
            path=relative_path, is_file=is_file, size=size, access='unpublished'
        )

    @contextmanager
    def raw_file(self, file_path: str, *args, **kwargs):
        assert is_safe_relative_path(file_path)

        full_path = self.raw_file_object(file_path).os_path

        try:
            with self._fs.open(full_path, *args, **kwargs) as f:
                yield f
        except (FileNotFoundError, IsADirectoryError) as e:
            raise KeyError(full_path) from e

    def raw_file_size(self, file_path: str) -> int:
        assert is_safe_relative_path(file_path)
        return self._fs.size(self.raw_file_object(file_path).os_path)

    def raw_file_object(self, file_path: str) -> PathObject:
        assert is_safe_relative_path(file_path)
        return self._raw_dir.join_file(file_path)

    def archive_hdf5_location(self, entry_id: str) -> str:
        return self.join_dir('archive').join_file(f'{entry_id}.h5').os_path

    def write_archive(self, entry_id: str, data: Any) -> int:
        """Writes the data as archive file and returns the archive file size."""
        archive_file_object = self._archive_file_object(entry_id)
        try:
            write_archive(archive_file_object.os_path, 1, data=[(entry_id, data)])
        except Exception:
            # in case of failure, remove the possible corrupted archive file
            archive_file_object.delete()

            raise

        return archive_file_object.size

    @contextmanager
    def read_archive(self, entry_id: str) -> Iterator[ArchiveReader]:
        try:
            with read_archive(
                self._archive_file_object(entry_id, True).os_path
            ) as archive:
                yield archive
        except FileNotFoundError as e:
            raise KeyError(entry_id) from e

    def _archive_file_object(self, entry_id: str, fallback: bool = False) -> PathObject:
        def versioned_file_name(version_suffix):
            return f'{entry_id}{version_suffix}.msg'

        return _versioned_archive_file_object(
            self._archive_dir, versioned_file_name, fallback=fallback
        )

    def add_rawfiles(
        self,
        path: str,
        target_dir: str = '',
        cleanup_source_file_and_dir: bool = False,
        updated_files: set[str] | None = None,
        auto_decompress: bool = True,
    ) -> None:
        """Adds files or directories to the upload, optionally decompressing archives.

        If `path` refers to an archive (ZIP, TAR) and `auto_decompress` is True,
        the archive is extracted before merging. Otherwise, archives are treated as single files.

        Args:
            path (str): Path to the file or directory to add.
            target_dir (str, optional): Relative path within the upload's raw directory.
                Defaults to "".
            cleanup_source_file_and_dir (bool, optional): If True, deletes the source path
                and its parent directory after processing. Defaults to False.
            updated_files (set[str], optional): Set to track paths of files updated or added.
            auto_decompress (bool, optional): If True, automatically decompress archives.
                Defaults to True.

        Raises:
            AssertionError: If file format is unrecognized or merge conflicts occur.
        """
        assert not self.is_frozen
        assert self._fs.exists(path), f'{path} does not exist'
        assert is_safe_relative_path(target_dir)

        archive_format = get_compression_format(path) if auto_decompress else None
        if archive_format == 'error':
            raise ValueError('Bad archive.')

        src_fs: AbstractFileSystem
        if archive_format == 'zip':
            src_fs = ZipFileSystem(path)
            src_root = src_parent = ''
        elif archive_format == 'tar':
            src_fs = TarFileSystem(path)
            src_root = src_parent = ''
        else:
            src_fs = self._fs
            src_root = path
            src_parent = os.path.dirname(path) if self._fs.isfile(path) else path

        dst_root = os.path.join(self._raw_dir.os_path, target_dir)

        try:
            for item, info in src_fs.find(src_root, None, True, True).items():
                rel_path = os.path.relpath(item, src_parent)
                dst_path = os.path.join(dst_root, rel_path)
                if info['type'] == 'file':
                    if self._fs.exists(dst_path) and not self._fs.isfile(dst_path):
                        raise ValueError(
                            f'Cannot merge a file with a directory or vice versa: {rel_path}.'
                        )
                    if src_fs is not self._fs or item != dst_path:
                        self._fs.mkdirs(os.path.dirname(dst_path), exist_ok=True)
                        with (
                            src_fs.open(item) as src_f,
                            self._fs.open(dst_path, 'wb') as dst_f,
                        ):
                            shutil.copyfileobj(src_f, dst_f)

                    if updated_files is not None:
                        updated_files.add(os.path.join(target_dir, rel_path))
                elif info['type'] == 'directory':
                    if self._fs.exists(dst_path) and not self._fs.isdir(dst_path):
                        raise ValueError(
                            f'Cannot merge a file with a directory or vice versa: {rel_path}.'
                        )
                    self._fs.mkdirs(dst_path, True)
        finally:
            if hasattr(src_fs, 'close'):
                src_fs.close()

            if cleanup_source_file_and_dir:
                self._fs.rm(path, recursive=True)
                if self._fs.exists(parent := os.path.dirname(path)) and not self._fs.ls(
                    parent, False
                ):
                    self._fs.rm(parent, recursive=True)

            self._size = self._fs.du(dst_root)

    def delete_rawfiles(self, path, updated_files: set[str] | None = None):
        assert is_safe_relative_path(path)
        raw_os_path = UPath(self.os_path) / 'raw'
        os_path = raw_os_path / path
        if not self._fs.exists(os_path):
            return
        if updated_files is not None:
            updated_files.update(
                os.path.relpath(target, raw_os_path)
                for target in self._fs.find(os_path)
            )
        self._fs.rm(os_path, recursive=True)
        if raw_os_path == os_path:
            # Special case - deleting everything, i.e. the entire raw folder. Need to recreate.
            self._fs.makedirs(os_path)

    def copy_or_move_rawfile(
        self, src: str, dest: str, copy_or_move, updated_files: set[str] | None = None
    ):
        assert is_safe_relative_path(src)
        assert is_safe_relative_path(dest)
        src_full_path = os.path.join(self._raw_dir.os_path, src)
        dest_full_path = os.path.join(self._raw_dir.os_path, dest)
        if not self._fs.exists(src_full_path):
            return
        if not self._fs.isfile(src_full_path):
            raise ValueError('Copying a directory is not possible.')
        if self._fs.exists(dest_full_path):
            raise ValueError('A file with the same name already exists.')

        if copy_or_move.lower() == 'copy':
            self._fs.cp_file(src_full_path, dest_full_path)
        elif copy_or_move.lower() == 'move':
            self._fs.mv(src_full_path, dest_full_path)

        if updated_files is not None:
            updated_files.add(dest)
            # if both the new and old name are the same then no new entry will be
            # added to the set. but if different, we add the old one so that later on
            # when self.matchall is called in data.py, the old filename is removed
            # from mongo database
            updated_files.add(src)

    def metadata_file_cached(self, path_dir: str = ''):
        """
        Gets the content of the metadata file located in the directory defined by `path_dir`.
        The `path_dir` should be relative to the `raw` folder.
        """

        def json_load(_f):
            return json.load(f)

        def yaml_load(_f):
            return yaml.safe_load(_f)

        def dummy_load(_):
            return {}

        loader = {'json': json_load, 'yaml': yaml_load, 'yml': yaml_load}

        base = UPath(self._raw_dir.os_path) / path_dir
        for ext in config.process.metadata_file_extensions:
            if not self._fs.isfile(
                full_path := base / f'{config.process.metadata_file_name}.{ext}'
            ):
                continue
            try:
                with self._fs.open(full_path.as_posix()) as f:
                    return loader.get(ext, dummy_load)(f)
            except Exception as e:
                # ignore the file contents if the file is not parsable, just warn.
                self.logger.warn(
                    'could not parse nomad.yaml/json', path=path_dir, exc_info=e
                )
        return {}

    @property
    def is_frozen(self) -> bool:
        """Returns True if this upload is already *bagged*."""
        return self._frozen_file.exists()

    def pack(
        self,
        entries: list[datamodel.EntryMetadata],
        with_embargo: bool,
        create: bool = True,
        include_raw: bool = True,
        include_archive: bool = True,
    ) -> None:
        """
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
        """
        self.logger.info('started to pack upload')

        # freeze the upload
        assert not self.is_frozen, 'Cannot pack an upload that is packed, or packing.'
        with self._fs.open(self._frozen_file.os_path, 'w') as f:
            f.write('frozen')

        # Check embargo flag consistency
        for entry in entries:
            assert entry.with_embargo == with_embargo

        access = 'restricted' if with_embargo else 'public'
        other_access = (
            'public' if with_embargo else 'restricted'
        )  # The "inverted" access

        # Get or create a target dir in the public area
        target_dir = DirectoryObject(
            PublicUploadFiles.base_folder_for(self.upload_id), create=create
        )
        if os.listdir(target_dir.os_path):
            # Target dir contains files. Check that the target access is identical
            assert PublicUploadFiles(self.upload_id).access == access, (
                'Inconsistent access'
            )

        # zip archives
        if include_archive:
            with utils.timer(self.logger, 'packed msgpack archive') as log_data:
                log_data.update(
                    number_of_entries=self._pack_archive_files(
                        target_dir, list(entry.entry_id for entry in entries), access
                    )
                )
                target_dir.msg_fp(other_access).delete()
                target_dir.h5_fp(other_access).delete()

        # zip raw files
        if include_raw:
            with utils.timer(self.logger, 'packed raw files'):
                self._pack_raw_files(target_dir, access)
                target_dir.zip_fp(other_access).delete()

    def _pack_archive_files(
        self, target_dir: DirectoryObject, entries: list[str], access: str
    ):
        number_of_entries = len(entries)

        def create_iterator():
            for item in entries:
                fo = self._archive_file_object(item)
                yield item, fo.os_path if fo.exists() else None

        try:
            combine_archive(
                target_dir.msg_fp(access).os_path, number_of_entries, create_iterator()
            )

            import h5py

            with h5py.File(target_dir.h5_fp(access).os_path, 'w') as hdf5_target:
                for entry_id in entries:
                    with h5py.File(
                        self.archive_hdf5_location(entry_id), 'a'
                    ) as hdf5_source:
                        group = hdf5_target.create_group(entry_id)
                        for key in hdf5_source.keys():
                            hdf5_source.copy(key, group)
        except Exception as e:
            self.logger.error('exception during packing archives', exc_info=e)
            raise

        return number_of_entries

    def _pack_raw_files(self, target_dir: DirectoryObject, access: str):
        try:
            with zipfile.ZipFile(
                target_dir.zip_fp(access).os_path, mode='w'
            ) as raw_zip:
                for path_info in self.raw_listdir(recursive=True):
                    basename = os.path.basename(path_info.path)
                    # TODO remove extra handling of POTCAR files once processed uploads are published.
                    if basename.startswith('POTCAR'):
                        if not basename.endswith('.stripped'):
                            continue  # Skip the unstripped POTCAR files when publishing
                        if basename.endswith('.stripped.stripped'):
                            continue  # Skip redundantly stripped POTCAR files (created due to bug #979) when publishing
                    raw_zip.write(
                        self._raw_dir.join_file(path_info.path).os_path, path_info.path
                    )
        except Exception as e:
            self.logger.error('exception during packing raw files', exc_info=e)
            raise

    def entry_files(
        self, mainfile: str, with_mainfile: bool = True, with_cutoff: bool = True
    ) -> Iterable[str]:
        """
        Returns all the auxfiles and mainfile for a given mainfile. This implements
        nomad's logic about what is part of an entry and what not. The mainfile
        is the first element, the rest is sorted.
        Arguments:
            mainfile: The mainfile path relative to upload
            with_mainfile: Do include the mainfile, default is True
        """
        mainfile_object = self._raw_dir.join_file(mainfile)
        if not mainfile_object.exists():
            raise KeyError(mainfile)

        mainfile_basename = os.path.basename(mainfile)
        entry_dir = os.path.dirname(mainfile_object.os_path)
        entry_relative_dir = entry_dir[len(self._raw_dir.os_path) + 1 :]

        file_count = 0
        aux_files: list[str] = []
        dir_elements = os.listdir(entry_dir)
        dir_elements.sort()
        for dir_element in dir_elements:
            if dir_element != mainfile_basename and os.path.isfile(
                os.path.join(entry_dir, dir_element)
            ):
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
        """
        Calculates a hash for the given entry based on file contents and aux file contents.
        Arguments:
            mainfile: The mainfile path relative to the upload that identifies the entry in
                the folder structure.
            mainfile_key: The mainfile_key of the entry (if any)
        Returns:
            The calculated hash
        Raises:
            KeyError: If the mainfile does not exist.
        """
        hash = hashlib.sha512()
        for filepath in self.entry_files(mainfile):
            with self._fs.open(self._raw_dir.join_file(filepath).os_path) as f:
                for data in iter(lambda: f.read(65536), b''):
                    hash.update(data)
        if mainfile_key:
            hash.update(mainfile_key.encode('utf8'))
        return utils.make_websave(hash)

    def files_to_bundle(
        self, export_settings: BundleExportSettings
    ) -> Iterable[FileSource]:
        # Defines files for upload bundles of staging uploads.
        if export_settings.include_raw_files:
            yield DiskFileSource(self.os_path, 'raw')
        if export_settings.include_archive_files:
            yield DiskFileSource(self.os_path, 'archive')

    @classmethod
    def files_from_bundle(
        cls,
        bundle_file_source: BrowsableFileSource,
        import_settings: BundleImportSettings,
    ) -> Iterable[FileSource]:
        # Files to import for a staging upload
        if import_settings.include_raw_files:
            yield bundle_file_source.child('raw')
        if import_settings.include_archive_files:
            yield bundle_file_source.child('archive')
        if import_settings.include_bundle_info:
            yield bundle_file_source.child(bundle_info_filename)


class PublicUploadFiles(UploadFiles):
    @classmethod
    def _file_area(cls):
        return UPath(config.fs.public)

    @property
    def external_os_path(self):
        if not config.fs.public_external:
            return self.os_path

        return self.os_path.replace(config.fs.public, config.fs.public_external)

    @cached_property
    def access(self):
        """
        Which "access" is used, either 'public' (uploads without embargo) or 'restricted'
        (uploads with embargo). This is reflected in the names of the files holding the
        raw data and the archive data. The reason for this is so that it should be easy to
        see, by just looking at the files, if a published upload is embargoed or not.

        The access is determined by inspecting which files exist/contain data. If both
        public and restricted files exist/contain data, or if neither exists/contain data,
        a KeyError will be thrown (this should not happen if the upload is correctly packed).
        The inspection of the files is only done on the first call, and the cached result
        is used in subsequent calls. The only way to change the access is to call :func:`re_pack`.
        """
        # Determine access by inspecting the files
        files_found = False
        sole_access = None
        for access in ('public', 'restricted'):
            raw_zip_file_object = self.zip_fp(access)
            archive_msg_file_object = self.msg_fp(access)
            archive_hdf5_file_object = self.h5_fp(access)
            found = (
                (
                    raw_zip_file_object.exists()
                    and raw_zip_file_object.size > empty_zip_file_size
                )
                or (
                    archive_msg_file_object.exists()
                    and archive_msg_file_object.size > empty_archive_file_size
                )
                or (
                    archive_hdf5_file_object.exists()
                    and archive_hdf5_file_object.size > empty_hdf5_file_size
                )
            )
            if found:
                if files_found:
                    raise KeyError(
                        'Inconsistency: both public and restricted files found'
                    )
                files_found = True
                sole_access = access

        if not files_found:
            raise KeyError('Neither public nor restricted files found')

        return sole_access

    def raw_zip_file_object(self) -> PathObject:
        """
        Gets the raw zip file, either public or restricted, depending on which one is used.
        If both public and restricted files exist, or if none of them exist, a KeyError will
        be thrown.
        """
        return self.zip_fp(self.access)

    @contextmanager
    def _zip_fs(self):
        zip_fs = ZipFileSystem(self.raw_zip_file_object().os_path)
        yield zip_fs
        zip_fs.close()

    def archive_hdf5_location(self, entry_id: str) -> str:
        hdf5_file_object = self.h5_fp(self.access)
        if not hdf5_file_object.exists():
            raise FileNotFoundError()

        return hdf5_file_object.os_path

    @property
    def _missing_raw_files(self):
        return not self._fs.exists(self.raw_zip_file_object().os_path)

    @contextmanager
    def _open_msg_file(self) -> Iterator[ArchiveReader]:
        msg_file_object = self.msg_fp(self.access, fallback=True)

        if not msg_file_object.exists():
            raise FileNotFoundError()

        with read_archive(msg_file_object.os_path) as archive:
            yield archive

    def to_staging(
        self, create: bool = False, include_archive: bool = False
    ) -> StagingUploadFiles | None:
        if StagingUploadFiles.exists_for(self.upload_id):
            if create:
                raise FileExistsError('Staging upload does already exist')
            return StagingUploadFiles(self.upload_id)

        if not create:
            return None

        staging_upload_files = StagingUploadFiles(self.upload_id, create=True)
        if (raw_zip_file := self.raw_zip_file_object()).exists():
            staging_upload_files.add_rawfiles(raw_zip_file.os_path)

        if include_archive:
            with suppress(FileNotFoundError):
                with self._open_msg_file() as archive:
                    for entry_id, data in archive.items():
                        staging_upload_files.write_archive(
                            entry_id.strip(), to_json(data)
                        )

                import h5py

                with h5py.File(self.archive_hdf5_location('')) as hdf5_source:
                    for entry_id, data in hdf5_source.items():
                        with h5py.File(
                            staging_upload_files.archive_hdf5_location(entry_id), 'w'
                        ) as hdf5_target:
                            for key in data.keys():
                                data.copy(key, hdf5_target)

        return staging_upload_files

    def is_empty(self) -> bool:
        with self._zip_fs() as zip_fs:
            return not zip_fs.ls('', False)

    def raw_exists(self, path: str) -> bool:
        if not is_safe_relative_path(path):
            return False
        if self._missing_raw_files:
            # We consider the empty path (i.e. root) to always "exists".
            return not path
        with self._zip_fs() as zip_fs:
            return zip_fs.exists(path)

    def raw_isfile(self, path: str) -> bool:
        if not is_safe_relative_path(path) or self._missing_raw_files:
            return False
        with self._zip_fs() as zip_fs:
            return zip_fs.isfile(path)

    def raw_listdir(
        self,
        path: str = '',
        recursive=False,
        files_only=False,
        depth: int = -1,
    ) -> Iterable[RawPathInfo]:
        if not is_safe_relative_path(path) or depth == 0:
            return
        if not path and self._missing_raw_files:
            return

        with self._zip_fs() as zip_fs:
            for target in zip_fs.find(
                path, (depth if depth > 0 else None) if recursive else 1, not files_only
            ):
                if not (isfile := zip_fs.isfile(target)) and UPath(target) == UPath(
                    path
                ):
                    # skip folder itself
                    continue
                yield RawPathInfo(
                    path=target,
                    is_file=isfile,
                    size=zip_fs.size(target) if isfile else zip_fs.du(target),
                    access=self.access,
                )

    @contextmanager
    def raw_file(self, file_path: str, *args, **kwargs):
        assert is_safe_relative_path(file_path)
        mode = kwargs.pop('mode', None)
        if len(args) > 0:
            mode = args[0]
        mode = mode or 'rb'
        encoding = kwargs.pop('encoding', None)

        try:
            with self._zip_fs() as zip_fs, zip_fs.open(file_path, **kwargs) as f:
                yield io.TextIOWrapper(f, encoding=encoding) if 't' in mode else f
        except (FileNotFoundError, IsADirectoryError, KeyError) as e:
            raise KeyError(file_path) from e

    def raw_file_size(self, file_path: str) -> int:
        assert is_safe_relative_path(file_path)

        with suppress(FileNotFoundError):
            with self._zip_fs() as zip_fs:
                if file_size := zip_fs.size(file_path):
                    return file_size

        raise KeyError(file_path)

    @contextmanager
    def read_archive(self, entry_id: str) -> Iterator[ArchiveReader]:
        try:
            with self._open_msg_file() as archive:
                if entry_id not in archive:
                    raise KeyError(entry_id)
                yield archive
        except FileNotFoundError as e:
            raise KeyError(entry_id) from e

    def re_pack(self, with_embargo: bool) -> None:
        """
        Repacks the files when changing the embargo flag on the upload. That is: when lifting the
        embargo the file names of the raw zip file and the archive file change from containing
        the keyword "restricted" to "public". Adding embargo to a non-embargoed published
        upload is also supported, but only admins should be allowed to do this. The existing
        files are just renamed, so this should be a rather quick operation. The upload must
        be correctly packed (i.e. there cannot be non-empty public and restricted files
        at the same time).
        """
        if (self.access == 'restricted') == with_embargo:
            return

        self.close()

        new_access = 'restricted' if with_embargo else 'public'

        def _move_to(src, dest):
            if src.exists():
                dest.delete()
                os.rename(src.os_path, dest.os_path)

        _move_to(self.msg_fp(self.access), self.msg_fp(new_access))
        _move_to(self.raw_zip_file_object(), self.zip_fp(new_access))
        _move_to(self.h5_fp(self.access), self.h5_fp(new_access))

        self.__dict__.pop('access', None)  # clear cached_property

    def files_to_bundle(
        self, export_settings: BundleExportSettings
    ) -> Iterable[FileSource]:
        # Defines files for upload bundles of published uploads.
        for filename in sorted(os.listdir(self.os_path)):
            if filename.startswith('raw-') and export_settings.include_raw_files:
                yield DiskFileSource(self.os_path, filename)
            if (
                filename.startswith('archive-')
                and export_settings.include_archive_files
            ):
                yield DiskFileSource(self.os_path, filename)

    @classmethod
    def files_from_bundle(
        cls,
        bundle_file_source: BrowsableFileSource,
        import_settings: BundleImportSettings,
    ) -> Iterable[FileSource]:
        for filename in bundle_file_source.find(''):
            if filename.startswith('raw-') and import_settings.include_raw_files:
                yield bundle_file_source.child(filename)
            if (
                filename.startswith('archive-')
                and import_settings.include_archive_files
            ):
                yield bundle_file_source.child(filename)
            if filename == bundle_info_filename and import_settings.include_bundle_info:
                yield bundle_file_source.child(filename)
