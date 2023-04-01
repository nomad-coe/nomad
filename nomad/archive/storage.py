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

from typing import Iterable, Any, Tuple, Dict, BinaryIO, Union, List, cast
from io import BytesIO, BufferedReader
from collections.abc import Mapping, Sequence

from memoization import cached
import msgpack
from msgpack.fallback import Packer, StringIO
import struct
import json

from nomad import utils
from nomad.config import archive

__packer = msgpack.Packer(autoreset=True, use_bin_type=True)

_toc_uuid_size = utils.default_hash_len + 1
_toc_item_size = _toc_uuid_size + 25  # packed(uuid + [10-byte-pos, 10-byte-pos])
_entries_per_block = archive.block_size // _toc_item_size
_bytes_per_block = _entries_per_block * _toc_item_size


def packb(o):
    return __packer.pack(o)


def unpackb(o):
    return msgpack.unpackb(o, raw=False)


def _encode(start: int, end: int) -> bytes:
    return start.to_bytes(5, byteorder='little', signed=False) + end.to_bytes(
        5, byteorder='little', signed=False)


def _decode(position: bytes) -> Tuple[int, int]:
    return int.from_bytes(position[:5], byteorder='little', signed=False), int.from_bytes(
        position[5:], byteorder='little', signed=False)


def _unpack_entry(data: bytes) -> Tuple[Any, Tuple[Any, Any]]:
    entry_uuid = unpackb(data[: _toc_uuid_size])
    positions_encoded = unpackb(data[_toc_uuid_size:])
    return entry_uuid, (_decode(positions_encoded[0]), _decode(positions_encoded[1]))


def _to_son(data):
    if isinstance(data, ArchiveList):
        data = data.to_list()

    elif isinstance(data, ArchiveDict):
        data = data.to_dict()

    # no need to convert build-in types
    return data


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
        self.toc_depth = toc_depth
        # noinspection PyTypeChecker
        self.toc: Dict[str, Any] = None
        self._depth = 0

        # Because we cannot change msgpacks interface of _pack, this _stack is used to
        # transfer the result of _pack calls in terms of the TOC.
        self._stack: List[Any] = []

        super().__init__(*args, **kwargs)

    def _pos(self):
        return self._buffer.getbuffer().nbytes

    def _pack_list(self, obj, *args, **kwargs):
        pack_result = super()._pack(obj, *args, **kwargs)

        toc_result = []
        # same assumption and condition as above
        if len(obj) > 0 and isinstance(obj[0], dict):
            for _ in obj:
                toc_result.append(self._stack.pop())

            self._stack.append(list(reversed(toc_result)))

        return pack_result

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
            for key, value in reversed(list(obj.items())):
                if isinstance(value, dict) or (isinstance(value, (list, tuple)) and len(
                        value) > 0 and isinstance(value[0], dict)):
                    # assumes non emptiness and uniformity of array items
                    toc[key] = self._stack.pop()

            toc_result['toc'] = {key: value for key, value in reversed(list(toc.items()))}

        end = self._pos()
        toc_result['pos'] = [start, end]

        self._stack.append(toc_result)

        return pack_result

    def _pack(self, obj, *args, **kwargs):
        if isinstance(obj, dict):
            return self._pack_dict(obj, *args, **kwargs)

        if isinstance(obj, list):
            return self._pack_list(obj, *args, **kwargs)

        return self._buffer.write(packb(obj))

    def pack(self, obj):
        assert isinstance(obj, dict), f'TOC packer can only pack dicts, {obj.__class__}'
        self._depth = 0
        self._buffer = StringIO()
        result = super().pack(obj)
        self.toc = self._stack.pop()
        assert len(self._stack) == 0
        return result


class ArchiveWriter:
    def __init__(self, file_or_path: Union[str, BytesIO], n_entries: int, entry_toc_depth: int):
        self.file_or_path = file_or_path
        self.n_entries = n_entries

        self._pos = 0
        # noinspection PyTypeChecker
        self._toc_position: Tuple[int, int] = None
        self._toc: Dict[str, Tuple[Tuple[int, int], Tuple[int, int]]] = {}
        # noinspection PyTypeChecker
        self._f: BinaryIO = None
        self._toc_packer = TOCPacker(toc_depth=entry_toc_depth)

    def __enter__(self):
        if isinstance(self.file_or_path, str):
            self._f = open(self.file_or_path, 'wb')
        elif isinstance(self.file_or_path, BytesIO):
            self._f = self.file_or_path
            self._f.seek(0)
        else:
            raise ValueError('not a file or path')

        # write empty placeholder header
        self._write_map_header(3)
        self._writeb('toc_pos')
        self._writeb(_encode(0, 0))

        self._writeb('toc')
        toc_start, _ = self._write_map_header(self.n_entries)
        _, toc_end = self._write(b'0' * _toc_item_size * self.n_entries)
        self._toc_position = toc_start, toc_end

        self._writeb('data')
        self._write_map_header(self.n_entries)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_val is not None:
            raise exc_val

        # go back and write the real TOC to the header
        self._f.seek(0)
        self._pos = 0

        assert len(self._toc) == self.n_entries
        toc_items = sorted(self._toc.items(), key=lambda item: item[0])
        toc = {uuid: [_encode(*positions[0]), _encode(
            *positions[1])] for uuid, positions in toc_items}

        self._write_map_header(3)
        self._writeb('toc_pos')
        self._writeb(_encode(*self._toc_position))

        self._writeb('toc')
        toc_position = self._writeb(toc)
        assert toc_position == self._toc_position, f'{toc_position} - {self._toc_position}'

        if isinstance(self.file_or_path, str):
            self._f.close()

    def _write_map_header(self, n):
        if n <= 0x0f:
            return self._write(struct.pack('B', 0x80 + n))
        if n <= 0xffff:
            return self._write(struct.pack(">BH", 0xde, n))
        if n <= 0xffffffff:
            return self._write(struct.pack(">BI", 0xdf, n))
        raise ValueError("Dict is too large")

    def _write(self, b: bytes) -> Tuple[int, int]:
        start = self._pos
        self._pos += self._f.write(b)
        return start, self._pos

    def _writeb(self, obj):
        return self._write(packb(obj))

    def add(self, uuid: str, data: Any) -> None:
        uuid = utils.adjust_uuid_size(uuid)

        self._toc_packer.reset()
        packed = self._toc_packer.pack(data)
        toc = self._toc_packer.toc

        self._writeb(uuid)
        self._write_map_header(2)
        self._writeb('toc')
        toc_pos = self._writeb(toc)
        self._writeb('data')
        data_pos = self._write(packed)

        self._toc[uuid] = (toc_pos, data_pos)


class ArchiveItem:
    def __init__(self, f: BytesIO, offset: int = 0):
        self._f = f
        self._offset = offset

    def _direct_read(self, size: int, offset: int):
        self._f.seek(offset)
        return self._f.read(size)

    def _read(self, position: Tuple[int, int]):
        start, end = position
        raw_data = self._direct_read(end - start, start + self._offset)
        return unpackb(raw_data)

    @cached(thread_safe=False, max_size=512)
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

    def to_list(self):
        if self._full_list is None:
            self._full_list = [_to_son(self._child(child_entry)) for child_entry in self._toc_entry]

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
            return self.to_dict().__getitem__(key)

    def __iter__(self):
        return self._toc_entry['toc'].__iter__()

    def __len__(self):
        return self._toc_entry['toc'].__len__()

    @cached(thread_safe=False)
    def to_dict(self):
        # todo: potential bug here, what if children are ArchiveList or ArchiveDict?
        return self._read(self._toc_entry['pos'])


class ArchiveReader(ArchiveDict):
    def __init__(self, file_or_path: Union[str, BytesIO], use_blocked_toc=True):
        self._file_or_path = file_or_path

        if isinstance(self._file_or_path, str):
            f: BytesIO = cast(
                BytesIO, open(self._file_or_path, 'rb', buffering=archive.read_buffer_size))
        elif isinstance(self._file_or_path, (BytesIO, BufferedReader)):
            f = cast(BytesIO, self._file_or_path)
        else:
            raise ValueError('not a file or path')

        # noinspection PyTypeChecker
        super().__init__(None, f)

        # this number is determined by the msgpack encoding of the file beginning:
        # { 'toc_pos': <...>
        #              ^11
        self._toc_position = _decode(self._direct_read(10, 11))

        self._use_blocked_toc = use_blocked_toc

        if not self._use_blocked_toc:
            self._toc_entry = self._read(self._toc_position)
            return

        toc_start = self._toc_position[0]
        b = self._direct_read(1, toc_start)[0]
        if b & 0b11110000 == 0b10000000:
            self._toc_number = b & 0b00001111
            self._toc_offset = toc_start + 1
        elif b == 0xde:
            self._toc_number, = struct.unpack_from(">H", self._direct_read(2, toc_start + 1))
            self._toc_offset = toc_start + 3
        elif b == 0xdf:
            self._toc_number, = struct.unpack_from(">I", self._direct_read(4, toc_start + 1))
            self._toc_offset = toc_start + 5
        else:
            raise ArchiveError('Archive top-level TOC is not a msgpack map (dictionary).')

        self._toc: Dict[str, Any] = {}
        self._toc_block_info = [None] * (self._toc_number // _entries_per_block + 1)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _load_toc_block(self, i_entry: int):
        i_block = i_entry // _entries_per_block

        if self._toc_block_info[i_block]:
            return self._toc_block_info[i_block]

        i_offset = i_block * _bytes_per_block + self._toc_offset
        block_data = self._direct_read(_bytes_per_block, i_offset)

        first, last = None, None

        entries_current_block = min(
            _entries_per_block, self._toc_number - i_block * _entries_per_block)

        offset: int = 0
        for i in range(entries_current_block):
            entry_uuid, positions = _unpack_entry(block_data[offset:offset + _toc_item_size])
            self._toc[entry_uuid] = positions
            offset += _toc_item_size

            if i == 0:
                first = entry_uuid

            if i + 1 == entries_current_block:
                last = entry_uuid

        self._toc_block_info[i_block] = (first, last)  # type: ignore

        return self._toc_block_info[i_block]

    def __getitem__(self, key):
        key = utils.adjust_uuid_size(key)

        if self._use_blocked_toc and self._toc_entry is None:
            if self._toc_number == 0:
                raise KeyError(key)

            positions = self._toc.get(key)
            # TODO use hash algorithm instead of binary search
            if positions is None:
                r_start = 0
                r_end = self._toc_number
                i_entry = None
                while r_start <= r_end:
                    m_entry = r_start + (r_end - r_start) // 2
                    if i_entry == m_entry:
                        break

                    i_entry = m_entry

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
            positions = self._toc_entry[key]
            toc_position = _decode(positions[0])
            data_position = _decode(positions[1])

        return ArchiveDict(self._read(toc_position), self._f, data_position[0])

    def __iter__(self):
        if self._toc_entry is None:
            # is not necessarily read when using blocked toc
            self._toc_entry = self._read(self._toc_position)

        return self._toc_entry.__iter__()

    def __len__(self):
        if self._toc_entry is None:
            # is not necessarily read when using blocked toc
            self._toc_entry = self._read(self._toc_position)

        return self._toc_entry.__len__()

    def close(self):
        if isinstance(self._file_or_path, str):
            self._f.close()

    def is_closed(self):
        return self._f.closed


def serialise_container(v):
    if isinstance(v, ArchiveList):
        return v.to_list()
    if isinstance(v, ArchiveDict):
        return v.to_dict()

    return v


def write_archive(
        path_or_file: Union[str, BytesIO], n_entries: int,
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

    The TOC of each entry will have the same structure than the data up to a certain
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


def read_archive(file_or_path: Union[str, BytesIO], **kwargs) -> ArchiveReader:
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
    def benchmark():
        from time import time
        import sys

        with open('archive_test.json') as f:
            example_data = json.load(f)

        size = 5000 if len(sys.argv) == 1 else int(sys.argv[1])
        access_every = 2
        example_archive = [(utils.create_uuid(), example_data) for _ in range(0, size)]
        example_uuid = example_archive[int(size / 2)][0]

        # this impl
        # create archive
        start = time()
        buffer = BytesIO()
        write_archive(buffer, len(example_archive), example_archive, entry_toc_depth=2)
        print('archive.py: create archive (1): ', time() - start)

        # read single entry from archive
        buffer = BytesIO(buffer.getbuffer())
        for use_blocked_toc in [False, True]:
            start = time()
            for _ in range(0, 23):
                read_archive(buffer, use_blocked_toc=use_blocked_toc)[example_uuid]['run']['system']
            print(
                f'archive.py: access single entry system (23), blocked {use_blocked_toc:d}: ',
                (time() - start) / 23)

        # read every n-ed entry from archive
        buffer = BytesIO(buffer.getbuffer())
        for use_blocked_toc in [False, True]:
            start = time()
            for _ in range(0, 23):
                with read_archive(buffer, use_blocked_toc=use_blocked_toc) as data:
                    for i, entry in enumerate(example_archive):
                        if i % access_every == 0:
                            data[entry[0]]['run']['system']
            print(
                f'archive.py: access every {access_every:d}-ed entry single entry system (23), '
                f'blocked {use_blocked_toc:d}: ', (time() - start) / 23)

        # just msgpack
        start = time()
        packb(example_archive)
        print('msgpack: create archive (1): ', time() - start)

    benchmark()
