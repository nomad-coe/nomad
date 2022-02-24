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

from collections.abc import Mapping, Sequence
from io import BufferedReader, BytesIO
import struct
from typing import Any, BinaryIO, Dict, Iterable, List, Tuple, Union, cast

import msgpack
from msgpack.fallback import Packer, StringIO

from nomad import utils

__packer = msgpack.Packer(autoreset=True, use_bin_type=True)

_block_size = 262144
_toc_uuid_size = utils.default_hash_len + 1
_toc_item_size = _toc_uuid_size + 25  # packed(uuid + [10-byte-pos, 10-byte-pos])
_toc_entries_per_block = _block_size // _toc_item_size
_toc_bytes_per_block = _toc_entries_per_block * _toc_item_size


def packb(o):
    return __packer.pack(o)


def unpackb(o):
    return msgpack.unpackb(o, raw=False)


def _encode(start: int, end: int) -> bytes:
    return start.to_bytes(5, byteorder='little', signed=False) + end.to_bytes(
        5, byteorder='little', signed=False)


def _decode(position: bytes) -> Tuple[int, int]:
    return int.from_bytes(position[0:5], byteorder='little', signed=False), int.from_bytes(
        position[5:], byteorder='little', signed=False)


def _unpack_entry(data: bytes) -> Tuple[Any, Tuple[Any, Any]]:
    entry_uuid = unpackb(data[: _toc_uuid_size])
    positions_encoded = unpackb(data[_toc_uuid_size: _toc_item_size])
    return entry_uuid, (_decode(positions_encoded[0]), _decode(positions_encoded[1]))


class ArchiveError(Exception):
    ''' An error that indicates a broken archive. '''
    pass


class TOCPacker(Packer):
    '''
    A special msgpack packer that records a TOC while packing.

    Uses a combination of the pure python msgpack fallback packer and the "real"
    c-based packing.
    '''

    def __init__(self, toc_depth: int, *args, **kwargs):
        self.toc_depth: int = toc_depth
        # noinspection PyTypeChecker
        self.toc: Dict[str, Any] = None

        self._depth = 0

        # Because we cannot change msgpack interface of _pack, this _stack is used to
        # transfer the result of _pack calls in terms of the TOC.
        self._stack: List[Any] = []

        super().__init__(*args, **kwargs)

    def pack(self, obj):
        assert isinstance(obj, dict), f'TOC packer can only pack dicts, {obj.__class__}'
        self._depth = 0
        self._buffer = StringIO()
        result = super().pack(obj)
        self.toc = self._stack.pop()
        assert len(self._stack) == 0
        return result

    def _pos(self):
        return self._buffer.getbuffer().nbytes

    def _pack_dict(self, obj, *args, **kwargs):
        toc_result = {}
        start = self._pos()

        if self._depth >= self.toc_depth:
            pack_result = self._buffer.write(packb(obj))
        else:
            self._depth += 1
            pack_result = super()._pack(obj, *args, **kwargs)
            self._depth -= 1

            toc = {}
            # TODO: upgrade to 3.8 to avoid reversed copy
            for key, value in reversed(list(obj.items())):
                if (isinstance(value, (list, tuple)) and len(value) > 0 and isinstance(
                        value[0], dict)) or isinstance(value, dict):
                    # assume non emptiness and uniformity of array items
                    toc[key] = self._stack.pop()

            toc_result['toc'] = {key: value for key, value in reversed(list(toc.items()))}

        end = self._pos()
        toc_result['pos'] = [start, end]

        self._stack.append(toc_result)

        return pack_result

    def _pack_list(self, obj, *args, **kwargs):
        pack_result = super()._pack(obj, *args, **kwargs)

        len_obj = len(obj)
        # same assumption and condition as above
        if len_obj > 0 and isinstance(obj[0], dict):
            # replace last len_obj items as a list, for example,
            # [1, 2, 3, 4, 5] -> [1, 2, 3, [4, 5]]
            # implementation uses list copies
            # toc_result = self._stack[-len_obj:]
            # self._stack = self._stack[:-len_obj]
            # self._stack.append(toc_result)
            # implementation minimised memory usage
            len_stack: int = len(self._stack)
            start: int = len_stack - len_obj
            end: int = len_stack + 1
            self._stack[start] = self._stack[start:end]
            del self._stack[start + 1:end]

        return pack_result

    def _pack(self, obj, *args, **kwargs):
        if isinstance(obj, dict):
            return self._pack_dict(obj, *args, **kwargs)

        if isinstance(obj, list):
            return self._pack_list(obj, *args, **kwargs)

        return self._buffer.write(packb(obj))


class ArchiveWriter:
    def __init__(self, file_or_path: Union[str, BytesIO], n_entries: int, entry_toc_depth: int):
        self._file_or_path = file_or_path
        self._n_entries = n_entries

        self._pos = 0
        # noinspection PyTypeChecker
        self._toc_position: Tuple[int, int] = None
        self._toc: Dict[str, Tuple[Tuple[int, int], Tuple[int, int]]] = {}
        # noinspection PyTypeChecker
        self._f: BinaryIO = None
        self._toc_packer = TOCPacker(toc_depth=entry_toc_depth)

    def __enter__(self):
        if isinstance(self._file_or_path, str):
            self._f = open(self._file_or_path, 'wb', buffering=_block_size)
        elif isinstance(self._file_or_path, BytesIO):
            self._f = self._file_or_path
            self._f.seek(0)
        else:
            raise ValueError('not a file or path')

        # write empty placeholder header
        self._write_map_header(3)
        self._write('toc_pos')
        self._write(_encode(0, 0))

        self._write('toc')
        toc_start, _ = self._write_map_header(self._n_entries)
        _, toc_end = self._writeb(b'0' * _toc_item_size * self._n_entries)
        self._toc_position = toc_start, toc_end

        self._write('data')
        self._write_map_header(self._n_entries)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_val is not None:
            raise exc_val

        # go back and write the real TOC to the header
        self._f.seek(0)
        self._pos = 0

        assert len(self._toc) == self._n_entries

        self._write_map_header(3)
        self._write('toc_pos')
        self._write(_encode(*self._toc_position))

        self._write('toc')
        toc = {uuid: [_encode(*positions[0]), _encode(*positions[1])] for uuid, positions in
               sorted(self._toc.items(), key=lambda item: item[0])}
        toc_position = self._write(toc)

        assert toc_position == self._toc_position, f'{toc_position} - {self._toc_position}'

        if isinstance(self._file_or_path, str):
            self._f.close()

    # noinspection SpellCheckingInspection
    def _writeb(self, b: bytes) -> Tuple[int, int]:
        start = self._pos
        self._pos += self._f.write(b)
        return start, self._pos

    def _write(self, b: Any) -> Tuple[int, int]:
        return self._writeb(packb(b))

    def _write_map_header(self, n) -> Tuple[int, int]:
        if n <= 0x0f:
            return self._writeb(struct.pack('B', 0x80 + n))
        if n <= 0xffff:
            return self._writeb(struct.pack(">BH", 0xde, n))
        if n <= 0xffffffff:
            return self._writeb(struct.pack(">BI", 0xdf, n))

        raise ValueError("Dict is too large")

    def add(self, uuid: str, data: Any) -> None:
        uuid = utils.adjust_uuid_size(uuid)

        self._toc_packer.reset()
        packed = self._toc_packer.pack(data)
        toc = self._toc_packer.toc

        self._write(uuid)
        self._write_map_header(2)
        self._write('toc')
        toc_pos = self._write(toc)
        self._write('data')
        data_pos = self._writeb(packed)

        self._toc[uuid] = toc_pos, data_pos


class ArchiveItem:
    def __init__(self, toc_entry: Union[dict, list], f: BytesIO, offset: int = 0):
        self.toc_entry = toc_entry
        self._f = f
        self._offset = offset

    def _seek(self, offset: int) -> int:
        return self._f.seek(offset)

    def _direct_read(self, size: int) -> bytes:
        return self._f.read(size)

    def _read(self, position: Tuple[int, int]):
        start, end = position
        self._seek(start + self._offset)
        return unpackb(self._direct_read(end - start))

    def _child(self, child_toc_entry):
        if isinstance(child_toc_entry, dict):
            if 'toc' in child_toc_entry:
                return ArchiveObject(child_toc_entry, self._f, self._offset)

            return self._read(child_toc_entry['pos'])

        if isinstance(child_toc_entry, list):
            return ArchiveList(child_toc_entry, self._f, self._offset)

        assert False, 'unreachable'


class ArchiveList(ArchiveItem, Sequence):
    def __getitem__(self, index):
        child_toc_entry = self.toc_entry[index]
        return self._child(child_toc_entry)

    def __len__(self):
        return self.toc_entry.__len__()


class ArchiveObject(ArchiveItem, Mapping):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._data = None

    def __getitem__(self, key):
        if self._data is None:
            try:
                child_toc_entry = self.toc_entry['toc'][key]
            except KeyError:
                return self.to_dict().__getitem__(key)

            return self._child(child_toc_entry)

        return self._data.__getitem__(key)

    def __iter__(self):
        if self._data is None:
            return self.toc_entry['toc'].__iter__()

        return self._data.__iter__()

    def __len__(self):
        if self._data is None:
            return self.toc_entry['toc'].__len__()

        return self._data.__len__()

    def to_dict(self):
        return self._read(self.toc_entry['pos'])


class ArchiveReader(ArchiveObject):
    def __init__(self, file_or_path: Union[str, BytesIO], use_blocked_toc: bool = True):
        self.file_or_path = file_or_path
        self.use_blocked_toc: bool = use_blocked_toc

        if isinstance(self.file_or_path, str):
            f = cast(BytesIO, open(self.file_or_path, 'rb', _block_size))
        elif isinstance(self.file_or_path, (BytesIO, BufferedReader)):
            f = self.file_or_path
        else:
            raise ValueError('not a file or path')

        super().__init__(None, f)

        # noinspection PyTypeChecker
        self.toc_position: Tuple[int, int] = None

        # noinspection PyTypeChecker
        self._toc: Dict[str, Any] = None
        # noinspection PyTypeChecker
        self._toc_block_info: List[Tuple[str, str]] = None
        # noinspection PyTypeChecker
        self._toc_offset: int = None
        # noinspection PyTypeChecker
        self._n_toc: int = None

        self.open()

    def __enter__(self):
        self.open()

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def open(self):
        # this number is determined by the msgpack encoding of the file beginning:
        # { 'toc_pos': <...>
        #              ^11
        self._seek(11)
        self.toc_position = _decode(self._direct_read(10))

        if not self.use_blocked_toc:
            self.toc_entry = self._read(self.toc_position)
            return

        self._toc = {}

        toc_start = self.toc_position[0]
        self._seek(toc_start)

        b = self._direct_read(1)[0]
        if b & 0b11110000 == 0b10000000:
            self._n_toc = b & 0b00001111
            self._toc_offset = toc_start + 1
        elif b == 0xde:
            self._n_toc, = struct.unpack_from(">H", self._direct_read(2))
            self._toc_offset = toc_start + 3
        elif b == 0xdf:
            self._n_toc, = struct.unpack_from(">I", self._direct_read(4))
            self._toc_offset = toc_start + 5
        else:
            raise ArchiveError('Archive top-level TOC is not a msgpack map (dictionary).')

        self._toc_block_info = [None] * (self._n_toc // _toc_entries_per_block + 1)

    def is_closed(self):
        if isinstance(self.file_or_path, str):
            return self._f.closed

        raise RuntimeError('not a file')

    def close(self):
        if isinstance(self.file_or_path, str):
            self._f.close()

    def _load_toc_block(self, i_entry: int):
        '''
        Load the i_entry-th block into self.toc.
        '''
        i_block = i_entry // _toc_entries_per_block

        if self._toc_block_info[i_block]:
            return self._toc_block_info[i_block]

        self._seek(i_block * _toc_bytes_per_block + self._toc_offset)
        block_data = self._direct_read(_toc_bytes_per_block)

        first, last = None, None

        total_entries = min(_toc_entries_per_block, self._n_toc - i_block * _toc_entries_per_block)

        offset = 0
        for i in range(total_entries):
            entry_uuid, positions = _unpack_entry(block_data[offset:offset + _toc_item_size])
            self._toc[entry_uuid] = positions
            offset += _toc_item_size

            if i == 0:
                first = entry_uuid

            if i + 1 == total_entries:
                last = entry_uuid

        self._toc_block_info[i_block] = (first, last)

        return self._toc_block_info[i_block]

    def __getitem__(self, key):
        key = utils.adjust_uuid_size(key)

        if self.use_blocked_toc and self.toc_entry is None:
            if self._n_toc == 0:
                raise KeyError(key)

            positions = self._toc.get(key)
            # TODO use hash algorithm instead of binary search
            if positions is None:
                r_start = 0
                r_end = self._n_toc
                i_entry = None
                while r_start <= r_end:
                    mid_entry = r_start + (r_end - r_start) // 2
                    if i_entry == mid_entry:
                        break

                    i_entry = mid_entry

                    first, last = self._load_toc_block(i_entry)

                    if key < first:
                        r_end = i_entry - 1
                    elif key > last:
                        r_start = i_entry + 1
                    else:
                        break

                positions = self._toc.get(key)
                if positions is None:
                    raise KeyError(key)

            toc_position, data_position = positions

        else:
            positions_encoded = self.toc_entry[key]
            toc_position = _decode(positions_encoded[0])
            data_position = _decode(positions_encoded[1])

        toc = self._read(toc_position)

        return ArchiveObject(toc, self._f, data_position[0])

    def __iter__(self):
        if self.toc_entry is None:
            # is not necessarily read when using blocked toc
            self.toc_entry = self._read(self.toc_position)

        return self.toc_entry.__iter__()

    def __len__(self):
        if self.toc_entry is None:
            # is not necessarily read when using blocked toc
            self.toc_entry = self._read(self.toc_position)

        return self.toc_entry.__len__()


def write_archive(path_or_file: Union[str, BytesIO], n_entries: int,
                  data: Iterable[Tuple[str, Any]], entry_toc_depth: int = 2) -> None:
    '''
    Writes a msgpack-based archive file. The file contents will be a valid msgpack-object.
    The data will contain extra table-of-contents (TOC) objects that map some keys to
    positions in the file. Data can be partially read from these positions and deserialized
    with msgpack.

    The data in the archive file will have the following layout:

    .. code-block:: python

        {
            'toc_pos': b[start, end],
            'toc': {
                entry_uuid: [b[start, end], b[start, end]], ...
            },
            'data': {
                entry_uuid: {
                    'toc': {
                        key: {
                            'pos': [start, end],
                            'toc': ...
                        },
                        key: [
                            {
                                'pos': [start, end]
                                'toc': ...
                            }, ...
                        ],
                        ...
                    },
                    'data': ...
                }, ...
            }
        }


    The top-level TOC will map entry_uuids to positions. The key 'toc_pos' holds the
    position of the entry TOC, the second ('toc') the position of each entry. These positions
    will be absolute positions in the file. The top-level TOC will be ordered by entry_uuid.
    The top-level TOC positions are 2*5byte encoded integers. This will give the top-level TOC a
    predictable layout and will allow to partially read this TOC.

    The TOC of each entry will have the same structure as the data up to a certain
    TOC depth. A TOC object will hold the position of the object it refers to (key 'pos')
    and further deeper TOC data (key 'toc'). Only data objects (dict instances) will
    have TOC objects and only object count towards the TOC depth. Positions in the entry
    TOCs are regular msgpack encoded integers.

    Arguments:
        path_or_file: A file path or file-like to the archive file that should be written.
        n_entries: The number of entries that will be added to the file.
        data: The file contents as an iterator of entry id, data tuples.
        entry_toc_depth: The depth of the table of contents in each entry. Only objects will
            count for calculating the depth.
    '''
    with ArchiveWriter(path_or_file, n_entries, entry_toc_depth=entry_toc_depth) as writer:
        for uuid, entry in data:
            writer.add(uuid, entry)


def read_archive(file_or_path: str, **kwargs) -> ArchiveReader:
    '''
    Allows to read a msgpack-based archive.

    Arguments:
        file_or_path: A file path or file-like to the archive file that should be read. The
            respective file has to be closed by the user. The returned obj supports the
            'with' statement and has a 'close' method.

    Returns:
        A mapping (dict-like) that can be used to access the archive data. The mapping
        will lazily load data as it is used. The mapping needs to be closed or used within
        a 'with' statement to free the underlying file resource after use.
    '''

    return ArchiveReader(file_or_path, **kwargs)


if __name__ == '__main__':
    pass
