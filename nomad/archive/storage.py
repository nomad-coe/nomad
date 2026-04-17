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
from __future__ import annotations

from collections.abc import Generator, Iterable, Mapping, Sequence
from io import BufferedReader, BytesIO
from typing import Any, cast

import msgspec

from nomad import utils
from nomad.config import config

_toc_uuid_size = utils.default_hash_len + 1
_toc_item_size = _toc_uuid_size + 25  # packed(uuid + [10-byte-pos, 10-byte-pos])
_entries_per_block = config.archive.block_size // _toc_item_size
_bytes_per_block = _entries_per_block * _toc_item_size


def unpackb(o):
    return msgspec.msgpack.decode(o)  # type: ignore


def _decode(position: bytes) -> tuple[int, int]:
    return int.from_bytes(
        position[:5], byteorder='little', signed=False
    ), int.from_bytes(position[5:], byteorder='little', signed=False)


def _unpack_entry(data: bytes) -> tuple[Any, tuple[Any, Any]]:
    entry_uuid = unpackb(data[:_toc_uuid_size])
    positions_encoded = unpackb(data[_toc_uuid_size:])
    return entry_uuid, (_decode(positions_encoded[0]), _decode(positions_encoded[1]))


def to_json(data):
    if hasattr(data, 'to_json'):
        return data.to_json()

    # no need to convert build-in types
    return data


class ArchiveError(Exception):
    """An error that indicates a broken archive."""

    pass


class ArchiveItem:
    def __init__(self, f: BytesIO, offset: int = 0):
        self._f = f
        self._offset = offset

    def _direct_read(self, size: int, offset: int):
        self._f.seek(offset)
        return self._f.read(size)

    def _read(self, position: tuple[int, int]):
        start, end = position
        raw_data = self._direct_read(end - start, start + self._offset)
        return unpackb(raw_data)

    def _child(self, child_toc_entry):
        if isinstance(child_toc_entry, dict):
            if child_toc_entry.get('toc', None):
                return ArchiveDict(child_toc_entry, self._f, self._offset)

            return self._read(child_toc_entry['pos'])

        if isinstance(child_toc_entry, list):
            return ArchiveList(child_toc_entry, self._f, self._offset)

        assert False, 'unreachable'


class ArchiveList(ArchiveItem, Sequence):
    def __init__(self, toc_entry: list, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._toc_entry = toc_entry
        self._full_list = None

    def __getitem__(self, index):
        if self._full_list is None:
            child_entry = self._toc_entry[index]
            return self._child(child_entry)

        return self._full_list[index]

    def __len__(self):
        return self._toc_entry.__len__()

    def to_json(self):
        if self._full_list is None:
            self._full_list = [
                to_json(self._child(child_entry)) for child_entry in self._toc_entry
            ]

        return self._full_list


class ArchiveDict(ArchiveItem, Mapping):
    def __init__(self, toc_entry: dict, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._toc_entry = toc_entry

    def __getitem__(self, key):
        try:
            child_toc_entry = self._toc_entry['toc'][key]
            return self._child(child_toc_entry)
        except KeyError:
            return self.to_json().__getitem__(key)

    def __iter__(self):
        return self.to_json().__iter__()

    def __len__(self):
        return self.to_json().__len__()

    def to_json(self):
        # todo: potential bug here, what if children are ArchiveList or ArchiveDict?
        return self._read(self._toc_entry['pos'])


class ArchiveReader(ArchiveDict):
    def __init__(self, file_or_path: str | BytesIO):
        self._file_or_path = file_or_path

        if isinstance(self._file_or_path, str):
            f: BytesIO = cast(
                BytesIO,
                open(
                    self._file_or_path, 'rb', buffering=config.archive.read_buffer_size
                ),
            )
        elif isinstance(self._file_or_path, BytesIO | BufferedReader):
            f = cast(BytesIO, self._file_or_path)
        else:
            raise ValueError('not a file or path')

        # noinspection PyTypeChecker
        super().__init__(None, f)

        # this number is determined by the msgpack encoding of the file beginning:
        # { 'toc_pos': <...>
        #              ^11
        self._toc_entry = self._read(_decode(self._direct_read(10, 11)))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _locate_position(self, key: str) -> tuple:
        positions = self._toc_entry[key]
        return _decode(positions[0]), _decode(positions[1])

    def __getitem__(self, key):
        toc_position, data_position = self._locate_position(utils.adjust_uuid_size(key))

        return ArchiveDict(self._read(toc_position), self._f, data_position[0])

    def get_raw(self, key: str) -> tuple[dict, Generator]:
        """
        Get raw bytes of the data and the TOC of the entry.
        This is used to read the data without decoding it.
        This is used in combining individual entries into a single archive file.
        """
        toc_position, data_position = self._locate_position(utils.adjust_uuid_size(key))

        def _iter(position: tuple[int, int]):
            start, end = position
            total_size = end - start
            while total_size > 0:
                size = min(total_size, config.archive.copy_chunk_size)
                yield self._direct_read(size, start + self._offset)
                start += size
                total_size -= size

        return self._read(toc_position), _iter(data_position)

    def __iter__(self):
        return self._toc_entry.__iter__()

    def __len__(self):
        return self._toc_entry.__len__()

    def close(self):
        if isinstance(self._file_or_path, str):
            self._f.close()

    def is_closed(self):
        return self._f.closed if isinstance(self._file_or_path, str) else True


def read_archive(file_or_path: str | BytesIO, **kwargs) -> ArchiveReader:
    """
    Allows to read a msgpack-based archive.

    Arguments:
        file_or_path: A file path or file-like to the archive file that should be read. The
            respective file has to be closed by the user. The returned obj supports the
            'with' statement and has a 'close' method.

    Returns:
        A mapping (dict-like) that can be used to access the archive data. The mapping
        will lazily load data as it is used. The mapping needs to be closed or used within
        a 'with' statement to free the underlying file resource after use.
    """
    from .storage_v2 import ArchiveReader as ArchiveReaderNew
    from .storage_v2 import ArchiveWriter as ArchiveWriterNew

    # todo: replace implementation to enable automatic conversion
    # if isinstance(file_or_path, str):
    #     from nomad.archive.converter import convert_archive
    #
    #     convert_archive(file_or_path, overwrite=True)
    #
    # return ArchiveReaderNew(file_or_path, **kwargs)

    if isinstance(file_or_path, str):
        with open(file_or_path, 'rb') as f:
            magic = f.read(ArchiveWriterNew.magic_len)
    else:
        file_or_path.seek(0)
        magic = file_or_path.read(ArchiveWriterNew.magic_len)
        file_or_path.seek(0)

    if magic == ArchiveWriterNew.magic:
        return ArchiveReaderNew(file_or_path, **kwargs)  # type: ignore

    return ArchiveReader(file_or_path)


def combine_archive(path: str, n_entries: int, data: Iterable[tuple]):
    from .storage_v2 import ArchiveReader as ArchiveReaderNew
    from .storage_v2 import ArchiveWriter as ArchiveWriterNew

    with ArchiveWriterNew(
        path, n_entries, toc_depth=config.archive.toc_depth
    ) as writer:
        for uuid, src_path in data:
            if not src_path:
                writer.add(uuid, {})
                continue

            with read_archive(src_path) as reader:
                if isinstance(reader, ArchiveReaderNew):
                    toc, data = reader.get_raw(uuid)
                    writer.add_raw(uuid, toc, data)
                else:
                    # rare case, old reader new writer, toc is not compatible, has to repack
                    writer.add(uuid, to_json(reader[uuid]))


if __name__ == '__main__':
    pass
