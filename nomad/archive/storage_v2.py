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

import struct
from io import BytesIO

import msgpack

from nomad import utils, config
from nomad.archive import ArchiveError

_packer = msgpack.Packer(autoreset=True, use_bin_type=True)


class Utility:
    # https://github.com/msgpack/msgpack/blob/master/spec.md#str-format-family
    toc_uuid_size = utils.default_hash_len + 1  # 1 byte for 0b101XXXXX
    # https://github.com/msgpack/msgpack/blob/master/spec.md#array-format-family
    # 0b1001XXXX for the array ==> 1 byte
    # 0xc4 0bXXXXXXXX for bin 8 ==> 2 bytes
    # 10 bytes for positions ==> 10 bytes
    # 25 ==> 1+(2+10)*2
    toc_item_size = toc_uuid_size + 25  # packed(uuid + [10-byte-pos, 10-byte-pos])
    entries_per_block = config.archive.block_size // toc_item_size
    bytes_per_block = entries_per_block * toc_item_size

    @staticmethod
    def encode(start: int, end: int) -> bytes:
        '''
        Encode the start and end offsets into byte string.
        '''
        return start.to_bytes(5, byteorder='little', signed=False) + end.to_bytes(
            5, byteorder='little', signed=False)

    @staticmethod
    def decode(position: bytes) -> tuple[int, int]:
        '''
        Decode the start and end offsets from a byte string.
        '''
        return int.from_bytes(position[:5], byteorder='little', signed=False), int.from_bytes(
            position[5:], byteorder='little', signed=False)

    # noinspection SpellCheckingInspection
    @staticmethod
    def packb(o) -> bytes:
        return _packer.pack(o)

    # noinspection SpellCheckingInspection
    @staticmethod
    def unpackb(o):
        return msgpack.unpackb(o, raw=False)

    @staticmethod
    def unpack_entry(data: bytes) -> tuple[str, tuple]:
        '''
        Unpack a single entry into the following format:
            entry_id, ((toc_start, toc_end), (data_start, data_end))
        '''
        entry_uuid: str = Utility.unpackb(data[: Utility.toc_uuid_size])
        positions_encoded: list = Utility.unpackb(data[Utility.toc_uuid_size:])
        return entry_uuid, (Utility.decode(positions_encoded[0]), Utility.decode(positions_encoded[1]))


class TOCPacker:
    '''
    A special msgpack packer that records a TOC while packing.
    '''

    def __init__(self, toc_depth: int, transform=None):
        self._toc_depth: int = toc_depth
        self._depth: int = 0
        self._buffer: BytesIO = None  # type: ignore
        self._transform = transform or (lambda x: x)

    @property
    def _pos(self):
        return self._buffer.getbuffer().nbytes

    def _pack(self, obj) -> dict:
        '''
        Pack a given object and record its TOC position in the buffer.
        '''

        def _pack_direct(_obj: bytes) -> None:
            '''
            Write bytes directly to the buffer.
            '''
            self._buffer.write(_obj)

        def _pack_raw(_obj) -> None:
            '''
            Pack a given object and write it to the buffer.
            '''
            _pack_direct(Utility.packb(_obj))

        if self._depth >= self._toc_depth or not isinstance(obj, (dict, list)):
            start_pos = self._pos
            _pack_raw(obj)
            return {'pos': [start_pos, self._pos]}

        start_pos = self._pos

        self._depth += 1

        obj_toc: dict | list
        all_small_obj: bool = False
        threshold: int = config.archive.small_obj_optimization_threshold
        if isinstance(obj, dict):
            _pack_direct(_packer.pack_map_header(len(obj)))
            obj_toc = {}
            for k, v in self._transform(obj.items()):
                _pack_raw(k)
                obj_toc[k] = self._pack(v)
            if all(v['pos'][1] < v['pos'][0] + threshold for v in obj_toc.values()):
                all_small_obj = True
        elif isinstance(obj, list):
            _pack_direct(_packer.pack_array_header(len(obj)))
            obj_toc = [self._pack(v) for v in obj]
            if all(v['pos'][1] < v['pos'][0] + threshold for v in obj_toc):
                all_small_obj = True
        else:
            raise ArchiveError(f'Expecting dict or list, got {obj.__class__}.')

        self._depth -= 1

        if self._pos < start_pos + threshold or all_small_obj:
            return {'pos': [start_pos, self._pos]}

        return {'toc': obj_toc, 'pos': [start_pos, self._pos]}

    def pack(self, obj):
        if not isinstance(obj, dict):
            raise ArchiveError(f'TOC packer can only pack dicts, {obj.__class__}')

        self._depth = 0
        self._buffer = BytesIO()
        toc: dict = self._pack(obj)
        return self._buffer.getvalue(), toc


class ArchiveWriter:
    '''
    Writes a msgpack-based archive file. The file contents will be a valid msgpack object.
    The resulting file contains extra table of contents (TOC) sections that map keys/indices to
    positions in the file. Data can be partially read from these positions.

    The data in the archive file will have the following layout:

    .. code-block:: python

        {
            'toc_pos': b[toc_start, toc_end],
            'toc': {
                entry_uuid_1: [b[toc_start, toc_end], b[data_start, data_end]],
                entry_uuid_2: [b[toc_start, toc_end], b[data_start, data_end]], ...
            },
            'data': {
                entry_uuid_1: {
                    'toc': {
                        key_1: { # for dicts
                            'pos': [start, end],
                            'toc': ...
                        },
                        key_2: { # for lists
                            'pos': [start, end],
                            'toc': [
                                {
                                    'pos': [start, end]
                                    'toc': ...
                                }, ...
                            ],
                        }, ...
                    },
                    'data': ...
                }, ...
            }
        }


    The top-level TOC will map entry_uuids to positions.
    The first key 'toc_pos' holds the position of the top-level TOC.
    The second key 'toc' holds TOC and data positions of each entry.
    These positions will be absolute positions in the file.
    The third key 'data' holds TOC and data for each entry.
    The top-level TOC will be ordered by entry_uuid.
    The top-level TOC positions are 2*5 bytes encoded integers.
    This will give the top-level TOC a predictable layout and will allow to partially read this TOC.

    The TOC of each entry will have the same structure for the data up to a certain TOC depth.
    A TOC section will hold the position of the object it refers to (key 'pos')
    and further deeper TOC data (key 'toc').
    Positions in the entry TOCs are regular msgpack encoded integers that represent offsets relative to the start.
    '''

    magic: bytes = b'nomad-archive-v2023'
    magic_len: int = len(magic)

    def __init__(self, file_or_path: str | BytesIO, n_entries: int, toc_depth: int):
        self.file_or_path: str | BytesIO = file_or_path
        self.n_entries: int = n_entries

        self._pos: int = 0
        self._toc_position: tuple[int, int] = None  # type: ignore
        self._toc: dict[str, tuple[tuple[int, int], tuple[int, int]]] = {}
        self._f: BytesIO = None  # type: ignore
        self._toc_packer: TOCPacker = TOCPacker(toc_depth=toc_depth)

    def __enter__(self):
        if isinstance(self.file_or_path, str):
            self._f = open(self.file_or_path, 'wb', buffering=config.archive.block_size)
        elif isinstance(self.file_or_path, BytesIO):
            self._f = self.file_or_path
            self._f.seek(0)
        else:
            raise ValueError('not a file or path')

        self._write(self.magic)
        # write empty placeholder header
        self._write_map_header(3)
        self._writeb('toc_pos')
        self._writeb(Utility.encode(0, 0))

        self._writeb('toc')
        start = self._write_map_header(self.n_entries)[0]
        end = self._write(b'0' * Utility.toc_item_size * self.n_entries)[1]
        self._toc_position = start, end

        self._writeb('data')
        self._write_map_header(self.n_entries)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_val is not None:
            self.close()
            raise exc_val

        assert len(self._toc) == self.n_entries

        # go back and write the real TOC to the header
        self._f.seek(self.magic_len)
        self._pos = self.magic_len

        self._write_map_header(3)
        self._writeb('toc_pos')
        self._writeb(Utility.encode(*self._toc_position))

        self._writeb('toc')
        toc_position = self._writeb({uuid: [
            Utility.encode(*pos[0]), Utility.encode(*pos[1])
        ] for uuid, pos in sorted(self._toc.items())})
        assert toc_position == self._toc_position, f'{toc_position} - {self._toc_position}'

        self.close()

    def close(self):
        if isinstance(self.file_or_path, str):
            self._f.close()

    def _write_map_header(self, n) -> tuple[int, int]:
        return self._write(_packer.pack_map_header(n))

    def _write(self, b: bytes) -> tuple[int, int]:
        start = self._pos
        self._pos += self._f.write(b)
        return start, self._pos

    # noinspection SpellCheckingInspection
    def _writeb(self, obj) -> tuple[int, int]:
        return self._write(Utility.packb(obj))

    def add(self, uuid: str, data):
        uuid = utils.adjust_uuid_size(uuid)

        packed, toc = self._toc_packer.pack(data)

        self._writeb(uuid)
        self._write_map_header(2)
        self._writeb('toc')
        toc_pos = self._writeb(toc)
        self._writeb('data')
        data_pos = self._write(packed)

        self._toc[uuid] = toc_pos, data_pos


def to_json(v):
    return v.to_json() if isinstance(v, ArchiveItem) else v


class ArchiveReadCounter:
    def __init__(self):
        self._counter = 0
        self._call_counter = 0

    def __iadd__(self, other):
        self._counter += other
        self._call_counter += 1
        return self

    def __call__(self, *args, **kwargs):
        return self._counter

    def bytes_per_call(self):
        return self._counter / self._call_counter

    def clear(self):
        self._counter = 0
        self._call_counter = 0


class ArchiveItem:
    def __init__(self, f: BytesIO, offset: int = 0, *, counter: ArchiveReadCounter = None):
        self._f: BytesIO = f
        self._offset: int = offset
        # to record how many bytes have been read
        self._counter: ArchiveReadCounter = counter
        # to record how many items have been accessed
        self._accessed_items: int = 0

    def __eq__(self, other):
        return self.to_json() == to_json(other)

    def __str__(self):
        return self.to_json().__str__()

    def _direct_read(self, size: int, offset: int):
        if self._f.closed:
            raise ArchiveError('Archive is closed')
        if self._counter:
            self._counter += size
        self._f.seek(offset)
        return self._f.read(size)

    def _read(self, position: tuple[int, int]):
        start, end = position
        return Utility.unpackb(self._direct_read(end - start, start + self._offset))

    def _child(self, toc: dict, offset: int = None):
        child_offset: int = offset or self._offset

        self._accessed_items += 1

        if (child_toc := toc.get('toc', None)) is None:
            start, end = toc['pos']
            if self._offset == 0:
                start += offset
                end += offset
            return self._read((start, end))

        if isinstance(child_toc, list):
            return ArchiveList(toc, self._f, child_offset, counter=self._counter)

        if isinstance(child_toc, dict):
            return ArchiveDict(toc, self._f, child_offset, counter=self._counter)

        raise ArchiveError(f'Invalid TOC: {toc}')

    def all(self):
        '''
        Suitable for read all only once.
        It will be overriden by reading and decoding the whole data at once.
        '''
        return self.to_json()

    def to_json(self):
        '''
        Ensure the result is JSON serializable.
        Suitable for multiple accesses.
        '''
        raise NotImplementedError


class ArchiveList(ArchiveItem):
    def __init__(self, toc: dict, f: BytesIO, offset: int = 0, *, counter: ArchiveReadCounter = None):
        super().__init__(f, offset, counter=counter)
        self._toc: list = toc['toc']
        self._pos: tuple = tuple(toc['pos'])
        self._cache = [None] * self._toc.__len__()
        self._full_cache: list = None  # type: ignore
        self._index: int = 0

    def __getitem__(self, index):
        if self._full_cache is not None:
            return self._full_cache[index]

        if isinstance(index, slice):
            index_range = range(*index.indices(self._toc.__len__()))
        elif isinstance(index, int):
            index_range = [index]
        else:
            raise TypeError(f'Invalid type: {type(index)} for index {index}')

        for i in index_range:
            if self._cache[i] is None:
                self._cache[i] = self._child(self._toc[i])

        return self._cache[index]

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index >= self._toc.__len__():
            raise StopIteration

        item = self[self._index]
        self._index += 1
        return item

    def __len__(self):
        return self._toc.__len__()

    def all(self):
        if self._full_cache is None:
            self._full_cache = self._read(self._pos)

        return self._full_cache

    def to_json(self):
        if self._full_cache is None:
            if config.archive.fast_loading and 2 * self._accessed_items < len(self._toc):
                self._full_cache = self._read(self._pos)
            else:
                self._full_cache = [to_json(v) for v in self]

        return self._full_cache


class ArchiveDict(ArchiveItem):
    def __init__(self, toc: dict, f: BytesIO, offset: int = 0, *, counter: ArchiveReadCounter = None):
        super().__init__(f, offset, counter=counter)
        self._toc: dict = toc['toc']
        self._pos: tuple = tuple(toc['pos'])
        self._cache: dict = {}
        self._full_cache: dict = None  # type: ignore

    def __getitem__(self, key):
        if key not in self._toc:
            raise KeyError(key)

        if self._full_cache is not None:
            return self._full_cache[key]

        if key not in self._cache:
            self._cache[key] = self._child(self._toc[key])

        return self._cache[key]

    def __contains__(self, item):
        return item in self._toc

    def __iter__(self):
        return self._toc.__iter__()

    def __len__(self):
        return self._toc.__len__()

    def get(self, key, default=None):
        return self[key] if key in self._toc else default

    def items(self):
        for k in self._toc:
            yield k, self[k]

    def keys(self):
        return self._toc.keys()

    def values(self):
        for k in self._toc:
            yield self[k]

    def all(self):
        if self._full_cache is None:
            self._full_cache = self._read(self._pos)

        return self._full_cache

    def to_json(self):
        if self._full_cache is None:
            if config.archive.fast_loading and 2 * self._accessed_items < len(self._toc):
                self._full_cache = self._read(self._pos)
            else:
                self._full_cache = {k: to_json(v) for k, v in self.items()}

        return self._full_cache


class ArchiveReader(ArchiveItem):
    def __init__(self, file_or_path: str | BytesIO, use_blocked_toc: bool = True, counter: ArchiveReadCounter = None):
        self._file_or_path: str | BytesIO = file_or_path

        if isinstance(self._file_or_path, str):
            f = open(self._file_or_path, 'rb', buffering=config.archive.read_buffer_size)
        elif isinstance(self._file_or_path, BytesIO):
            f = self._file_or_path
        else:
            raise ValueError('not a file or path')

        super().__init__(f, counter=counter)

        self._cache: dict = {}
        self._full_cache: dict = None  # type: ignore

        # this number is determined by the msgpack encoding of the file beginning:
        # { 'toc_pos': <...>
        #              ^11
        # 11 ==> 1 (0b0000XXXX for map) + 1 (0b101XXXXX for str key) + 7 ('toc_pos') + 2 (0xc4 0bXXXXXXXX for bin 8)
        self._toc_position = tuple(Utility.decode(self._direct_read(10, 11 + ArchiveWriter.magic_len)))
        self._toc_entry: dict | None = None

        self._use_blocked_toc: bool = use_blocked_toc

        # if not using blocked TOC, load everything now
        if not self._use_blocked_toc:
            self._ensure_toc()
            return

        # determine the map storing the TOC
        # https://github.com/msgpack/msgpack/blob/master/spec.md#map-format-family
        toc_start: int = self._toc_position[0]
        if (b := self._direct_read(1, toc_start)[0]) & 0b11110000 == 0b10000000:
            self._toc_number: int = b & 0b00001111
            self._toc_offset: int = toc_start + 1
        elif b == 0xde:
            self._toc_number, = struct.unpack_from(">H", self._direct_read(2, toc_start + 1))
            self._toc_offset = toc_start + 3
        elif b == 0xdf:
            self._toc_number, = struct.unpack_from(">I", self._direct_read(4, toc_start + 1))
            self._toc_offset = toc_start + 5
        else:
            raise ArchiveError('Top level TOC is not a msgpack map (dictionary).')

        self._toc: dict = {}
        self._toc_block_info: list = [None] * (self._toc_number // Utility.entries_per_block + 1)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        if exc_val:
            raise exc_val

    def _load_toc_block(self, i_entry: int) -> tuple:
        i_block: int = i_entry // Utility.entries_per_block

        if self._toc_block_info[i_block]:
            return self._toc_block_info[i_block]

        i_offset: int = i_block * Utility.bytes_per_block + self._toc_offset
        block_data = self._direct_read(Utility.bytes_per_block, i_offset)

        first, last = None, None

        offset: int = 0
        for i in range(entries_current_block := min(
                Utility.entries_per_block, self._toc_number - i_block * Utility.entries_per_block)):
            entry_uuid, positions = Utility.unpack_entry(block_data[offset:offset + Utility.toc_item_size])
            self._toc[entry_uuid] = positions
            offset += Utility.toc_item_size

            if i == 0:
                first = entry_uuid

            if i + 1 == entries_current_block:
                last = entry_uuid

        self._toc_block_info[i_block] = first, last

        return self._toc_block_info[i_block]

    def __getitem__(self, key: str) -> ArchiveDict:
        key = utils.adjust_uuid_size(key)

        if self._full_cache is not None:
            return self._full_cache[key]

        if key in self._cache:
            return self._cache[key]

        if self._use_blocked_toc and self._toc_entry is None:
            if self._toc_number == 0:
                raise KeyError(key)

            if (positions := self._toc.get(key)) is None:
                r_start, r_end, i_entry = 0, self._toc_number, -9999
                while r_start <= r_end:
                    if i_entry == (m_entry := r_start + (r_end - r_start) // 2):
                        break

                    i_entry = m_entry

                    first, last = self._load_toc_block(i_entry)

                    if key < first:
                        r_end = i_entry - 1
                    elif key > last:
                        r_start = i_entry + 1
                    else:
                        break

                if (positions := self._toc.get(key)) is None:
                    raise KeyError(key)

            toc_position, data_position = positions
        else:
            positions = self._toc_entry[key]
            toc_position = Utility.decode(positions[0])
            data_position = Utility.decode(positions[1])

        self._cache[key] = self._child(self._read(toc_position), data_position[0])  # type: ignore
        return self._cache[key]

    def __contains__(self, item):
        try:
            _ = self[item]
            return True
        except KeyError:
            return False

    def __iter__(self):
        self._ensure_toc()
        return self._toc_entry.__iter__()

    def __len__(self):
        self._ensure_toc()
        return self._toc_entry.__len__()

    def get(self, key: str, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def items(self):
        self._ensure_toc()
        for k in self._toc_entry:
            yield k, self[k]

    def keys(self):
        self._ensure_toc()
        return self._toc_entry.keys()

    def values(self):
        self._ensure_toc()
        for k in self._toc_entry:
            yield self[k]

    def _ensure_toc(self):
        if self._toc_entry is None:
            self._toc_entry = self._read(self._toc_position)

    def close(self, close_unowned: bool = False):
        if close_unowned or isinstance(self._file_or_path, str):
            self._f.close()

    def is_closed(self):
        # If the input is a BytesIO, it is assumed that the file is always closed
        # If the input is a path, need to check if the file is closed
        return self._f.closed if isinstance(self._file_or_path, str) else True

    def to_json(self):
        if self._full_cache is None:
            self._full_cache = {k: to_json(v) for k, v in self.items()}

        return self._full_cache


def write_archive(path_or_file: str | BytesIO, data: list, toc_depth: int = config.archive.toc_depth):
    with ArchiveWriter(path_or_file, len(data), toc_depth=toc_depth) as writer:
        for uuid, entry in data:
            writer.add(uuid, entry)
