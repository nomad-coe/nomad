# Copyright 2018 Markus Scheidgen, Alvin Noe Ladines
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

from typing import Iterable, Any, Tuple, Dict, BinaryIO, Union, List, cast
from io import BytesIO, BufferedReader
from collections.abc import Mapping, Sequence
import msgpack
from msgpack.fallback import Packer, StringIO
import struct
import json
import math
import re

from nomad import utils


__packer = msgpack.Packer(autoreset=True, use_bin_type=True)


def packb(o, **kwargs):
    return __packer.pack(o)


def unpackb(o, **kwargs):
    return msgpack.unpackb(o, raw=False)


def adjust_uuid_size(uuid):
    uuid = uuid.rjust(utils.default_hash_len, ' ')
    assert len(uuid) == utils.default_hash_len, 'uuids must have the right fixed size'
    return uuid


class ArchiveError(Exception):
    ''' An error that indicates a broken archive. '''
    pass


class ArchiveQueryError(Exception):
    '''
    An error that indicates that an archive query is either not valid or does not fit to
    the queried archive.
    '''
    pass


class TOCPacker(Packer):
    '''
    A special msgpack packer that records a TOC while packing.

    Uses a combination of the pure python msgpack fallback packer and the "real"
    c-based packing.
    '''
    def __init__(self, toc_depth: int, *args, **kwargs):
        self.toc_depth = toc_depth
        self.toc: Dict[str, Any] = None
        self._depth = 0

        # Because we cannot change msgpacks interface of _pack, this _stack is used to
        # tranfer the result of _pack calls in terms of the TOC.
        self._stack: List[Any] = []

        super().__init__(*args, **kwargs)

    def pack(self, obj, *args, **kwargs):
        assert isinstance(obj, dict), 'TOC packer can only pack dicts, %s' % obj.__class__
        self._depth = 0
        self._buffer = StringIO()
        result = super().pack(obj, *args, **kwargs)
        self.toc = self._stack.pop()
        assert len(self._stack) == 0
        return result

    def _pos(self):
        return self._buffer.getbuffer().nbytes

    def _pack(self, obj, *args, **kwargs):
        if isinstance(obj, dict):
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
                    if isinstance(value, (list, tuple)):
                        # this makes some assumptions about the non emptyness and uniformness
                        # of array items
                        if len(value) > 0 and isinstance(value[0], dict):
                            toc[key] = self._stack.pop()

                    elif isinstance(value, (dict, list, tuple)):
                        toc[key] = self._stack.pop()

                toc_result['toc'] = {key: value for key, value in reversed(list(toc.items()))}

            end = self._pos()
            toc_result['pos'] = [start, end]

            self._stack.append(toc_result)

        elif isinstance(obj, list):
            toc_result = []
            pack_result = super()._pack(obj, *args, **kwargs)

            # same assumption and condition as above
            if len(obj) > 0 and isinstance(obj[0], dict):
                for _ in obj:
                    toc_result.append(self._stack.pop())

                self._stack.append(list(reversed(toc_result)))

        else:
            pack_result = self._buffer.write(packb(obj))

        return pack_result


_toc_uuid_size = utils.default_hash_len + 1
_toc_item_size = _toc_uuid_size + 25  # packed(uuid + [10-byte-pos, 10-byte-pos])


class ArchiveWriter:
    def __init__(self, file_or_path: Union[str, BytesIO], n_entries: int, entry_toc_depth: int):
        self.file_or_path = file_or_path
        self.n_entries = n_entries

        self._pos = 0
        self._toc_position: Tuple[int, int] = None
        self._toc: Dict[str, Tuple[Tuple[int, int], Tuple[int, int]]] = {}
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
        self.write(packb('toc_pos'))
        self.write(packb(ArchiveWriter._encode_position(0, 0)))

        self.write(packb('toc'))
        toc_start, _ = self._write_map_header(self.n_entries)
        _, toc_end = self.write(b'0' * _toc_item_size * self.n_entries)
        self._toc_position = toc_start, toc_end

        self.write(packb('data'))
        self._write_map_header(self.n_entries)

        return self

    def write(self, b: bytes) -> Tuple[int, int]:
        start = self._pos
        self._pos += self._f.write(b)
        return start, self._pos

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_val is not None:
            raise exc_val

        # go back and write the real TOC to the header
        self._f.seek(0)
        self._pos = 0

        assert len(self._toc) == self.n_entries
        toc_items = sorted(self._toc.items(), key=lambda item: item[0])
        toc = {
            uuid: [
                ArchiveWriter._encode_position(*positions[0]),
                ArchiveWriter._encode_position(*positions[1])]
            for uuid, positions in toc_items}

        self._write_map_header(3)
        self.write(packb('toc_pos'))
        self.write(
            packb(ArchiveWriter._encode_position(*self._toc_position), use_bin_type=True))

        self.write(packb('toc'))
        toc_position = self.write(packb(toc, use_bin_type=True))
        assert toc_position == self._toc_position, '%s - %s' % (toc_position, self._toc_position)

        if isinstance(self.file_or_path, str):
            self._f.close()

    @staticmethod
    def _encode_position(start: int, end: int) -> bytes:
        return start.to_bytes(5, byteorder='little', signed=False) + \
            end.to_bytes(5, byteorder='little', signed=False)

    def add(self, uuid: str, data: Any) -> None:
        uuid = adjust_uuid_size(uuid)

        self._toc_packer.reset()
        packed = self._toc_packer.pack(data)
        toc = self._toc_packer.toc

        self.write(packb(uuid))
        self._write_map_header(2)
        self.write(packb('toc'))
        toc_pos = self.write(packb(toc, use_bin_type=True))
        self.write(packb('data'))
        data_pos = self.write(packed)

        self._toc[uuid] = (toc_pos, data_pos)

    def _write_map_header(self, n):
        if n <= 0x0f:
            return self.write(struct.pack('B', 0x80 + n))
        if n <= 0xffff:
            return self.write(struct.pack(">BH", 0xde, n))
        if n <= 0xffffffff:
            return self.write(struct.pack(">BI", 0xdf, n))
        raise ValueError("Dict is too large")


class ArchiveItem:
    def __init__(self, toc_entry: list, f: BytesIO, offset: int = 0):
        self.toc_entry = toc_entry
        self._f = f
        self._offset = offset

    def _read(self, position: Tuple[int, int]):
        start, end = position
        self._f.seek(start + self._offset)
        return unpackb(self._f.read(end - start))

    def _getchild(self, child_toc_entry):
        if isinstance(child_toc_entry, dict):
            if 'toc' in child_toc_entry:
                return ArchiveObject(child_toc_entry, self._f, self._offset)
            else:
                return self._read(child_toc_entry['pos'])

        elif isinstance(child_toc_entry, list):
            return ArchiveList(child_toc_entry, self._f, self._offset)

        else:
            assert False, 'unreachable'


class ArchiveList(ArchiveItem, Sequence):

    def __getitem__(self, index):
        child_toc_entry = self.toc_entry[index]
        return self._getchild(child_toc_entry)

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
            return self._getchild(child_toc_entry)
        else:
            return self._data.__getitem__(key)

    def __iter__(self):
        if self._data is None:
            return self.toc_entry['toc'].__iter__()
        else:
            return self._data.__iter__()

    def __len__(self):
        if self._data is None:
            return self.toc_entry['toc'].__len__()
        else:
            return self._data.__len__()

    def to_dict(self):
        return self._read(self.toc_entry['pos'])


class ArchiveReader(ArchiveObject):
    def __init__(self, file_or_path: Union[str, BytesIO], use_blocked_toc=True):
        self.file_or_path = file_or_path

        f: BytesIO = None
        if isinstance(self.file_or_path, str):
            f = cast(BytesIO, open(self.file_or_path, 'rb'))
        elif isinstance(self.file_or_path, BytesIO):
            f = self.file_or_path
        elif isinstance(self.file_or_path, BufferedReader):
            f = self.file_or_path
        else:
            raise ValueError('not a file or path')

        super().__init__(None, f)

        self._toc_entry = None

        # this number is determined by the msgpack encoding of the file beginning:
        # { 'toc_pos': <...>
        #              ^11
        self._f.seek(11)
        self.toc_position = ArchiveReader._decode_position(self._f.read(10))

        self.use_blocked_toc = use_blocked_toc
        if use_blocked_toc:
            self._f.seek(11)
            self._toc: Dict[str, Any] = {}
            toc_start = self.toc_position[0]
            self._f.seek(toc_start)
            b = self._f.read(1)[0]
            if b & 0b11110000 == 0b10000000:
                self._n_toc = b & 0b00001111
                self._toc_offset = toc_start + 1
            elif b == 0xde:
                self._n_toc, = struct.unpack_from(">H", self._f.read(2))
                self._toc_offset = toc_start + 3
            elif b == 0xdf:
                self._n_toc, = struct.unpack_from(">I", self._f.read(4))
                self._toc_offset = toc_start + 5
            else:
                raise ArchiveError('Archive top-level TOC is not a msgpack map (dictionary).')

        else:
            self.toc_entry = self._read(self.toc_position)

    def __enter__(self):
        return self

    fs_block_size = 4096
    toc_block_size_entries = math.ceil(4096 / _toc_item_size)
    toc_block_size_bytes = math.ceil(4096 / 54) * _toc_item_size

    def _load_toc_block(self, i_entry: int):
        i_block = math.floor(i_entry / ArchiveReader.toc_block_size_entries)

        self._f.seek(i_block * ArchiveReader.toc_block_size_bytes + self._toc_offset)
        block_data = self._f.read(ArchiveReader.toc_block_size_bytes)
        first, last = None, None
        toc_block_size_entries = min(
            ArchiveReader.toc_block_size_entries,
            self._n_toc - i_block * ArchiveReader.toc_block_size_entries)

        for i in range(0, toc_block_size_entries):
            offset = i * _toc_item_size
            entry_uuid = unpackb(block_data[offset:offset + _toc_uuid_size])
            positions_encoded = unpackb(block_data[offset + _toc_uuid_size:offset + _toc_item_size])
            positions = (
                ArchiveReader._decode_position(positions_encoded[0]),
                ArchiveReader._decode_position(positions_encoded[1]))
            self._toc[entry_uuid] = positions

            if i == 0:
                first = entry_uuid

            if i + 1 == toc_block_size_entries:
                last = entry_uuid

        return first, last

    def __getitem__(self, key):
        key = adjust_uuid_size(key)

        if self.use_blocked_toc and self.toc_entry is None:
            if self._n_toc == 0:
                raise KeyError(key)

            positions = self._toc.get(key)
            # TODO use hash algo instead of binary search
            if positions is None:
                r_start = 0
                r_end = self._n_toc
                i_block = None
                while r_start <= r_end:
                    new_i_block = r_start + math.floor((r_end - r_start) / 2)
                    if i_block == new_i_block:
                        break
                    else:
                        i_block = new_i_block

                    first, last = self._load_toc_block(i_block)

                    if key < first:
                        r_end = i_block - 1
                    elif key > last:
                        r_start = i_block + 1
                    else:
                        break

                positions = self._toc.get(key)
                if positions is None:
                    raise KeyError(key)

            toc_position, data_position = positions

        else:
            positions_encoded = self.toc_entry[key]
            toc_position = ArchiveReader._decode_position(positions_encoded[0])
            data_position = ArchiveReader._decode_position(positions_encoded[1])

        toc = self._read(toc_position)
        return ArchiveObject(toc, self._f, data_position[0])

    def __iter__(self):
        if self.toc_entry is None:
            # is not necessarely read when using blocked toc
            self.toc_entry = self._read(self.toc_position)

        return self.toc_entry.__iter__()

    def __len__(self):
        if self.toc_entry is None:
            # is not necessarely read when using blocked toc
            self.toc_entry = self._read(self.toc_position)

        return self.toc_entry.__len__()

    def close(self):
        if isinstance(self.file_or_path, str):
            self._f.close()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @staticmethod
    def _decode_position(position: bytes) -> Tuple[int, int]:
        return int.from_bytes(position[0:5], byteorder='little', signed=False), \
            int.from_bytes(position[5:], byteorder='little', signed=False)

    def is_closed(self):
        return self._f.closed


def write_archive(
        path_or_file: Union[str, BytesIO], n_entries: int, data: Iterable[Tuple[str, Any]],
        entry_toc_depth: int = 2) -> None:
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

    The TOC of each entry will have the same structure than the data upto a certain
    TOC depth. A TOC object will hold the position of the object it refers to (key 'pos')
    and further deeper TOC data (key 'toc'). Only data objects (dict instances) will
    have TOC objects and only object count towards the TOC depth. Positions in the entry
    TOCs are regular msgpack encoded integers.

    Arguments:
        file_or_path: A file path or file-like to the archive file that should be written.
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
        will lazyly load data as it is used. The mapping needs to be closed or used within
        a 'with' statement to free the underlying file resource after use.
    '''

    return ArchiveReader(file_or_path, **kwargs)


__query_archive_key_pattern = re.compile(r'^([\s\w\-]+)(\[([-?0-9]*)(:([-?0-9]*))?\])?$')


def query_archive(f_or_archive_reader: Union[str, ArchiveReader, BytesIO], query_dict: dict, **kwargs):
    def _fix_index(index, length):
        if index is None:
            return index
        if index < 0:
            return max(-(length), index)
        else:
            return min(length, index)

    def _to_son(data):
        if isinstance(data, (ArchiveList, List)):
            data = [_to_son(item) for item in data]

        elif isinstance(data, ArchiveObject):
            data = data.to_dict()

        return data

    def _load_data(query_dict: Dict[str, Any], archive_item: ArchiveObject, main_section: bool = False):
        if not isinstance(query_dict, dict):
            return _to_son(archive_item)

        result = {}
        for key, val in query_dict.items():
            key = key.strip()

            # process array indices
            match = __query_archive_key_pattern.match(key)
            index: Union[Tuple[int, int], int] = None
            if match:
                key = match.group(1)

                # check if we have indices
                if match.group(2) is not None:
                    first_index, last_index = None, None
                    group = match.group(3)
                    first_index = None if group == '' else int(group)

                    if match.group(4) is not None:
                        group = match.group(5)
                        last_index = None if group == '' else int(group)
                        index = (0 if first_index is None else first_index, last_index)

                    else:
                        index = first_index  # one item

                else:
                    index = None

            else:
                raise ArchiveQueryError('invalid key format: %s' % key)

            # support for shorter uuids
            if main_section:
                archive_key = adjust_uuid_size(key)
            else:
                archive_key = key

            try:
                archive_child = archive_item[archive_key]
                is_list = isinstance(archive_child, (ArchiveList, list))

                if index is None and is_list:
                    index = (0, None)
                elif index is not None and not is_list:
                    raise ArchiveQueryError('cannot use list key on none list %s' % key)

                if index is None:
                    pass
                else:
                    length = len(archive_child)
                    if isinstance(index, tuple):
                        index = (_fix_index(index[0], length), _fix_index(index[1], length))
                        if index[0] == index[1]:
                            archive_child = [archive_child[index[0]]]
                        else:
                            archive_child = archive_child[index[0]: index[1]]
                    else:
                        archive_child = [archive_child[_fix_index(index, length)]]

                if isinstance(archive_child, (ArchiveList, list)):
                    result[key] = [_load_data(val, item) for item in archive_child]
                else:
                    result[key] = _load_data(val, archive_child)

            except (KeyError, IndexError):
                continue

        return result

    if isinstance(f_or_archive_reader, ArchiveReader):
        return _load_data(query_dict, f_or_archive_reader, True)

    elif isinstance(f_or_archive_reader, (BytesIO, str)):
        with ArchiveReader(f_or_archive_reader, **kwargs) as archive:
            return _load_data(query_dict, archive, True)

    else:
        raise TypeError('%s is neither a file-like nor ArchiveReader' % f_or_archive_reader)


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
                read_archive(buffer, use_blocked_toc=use_blocked_toc)[example_uuid]['section_run']['section_system']
            print('archive.py: access single entry system (23), blocked %d: ' % use_blocked_toc, (time() - start) / 23)

        # read every n-ed entry from archive
        buffer = BytesIO(buffer.getbuffer())
        for use_blocked_toc in [False, True]:
            start = time()
            for _ in range(0, 23):
                with read_archive(buffer, use_blocked_toc=use_blocked_toc) as data:
                    for i, entry in enumerate(example_archive):
                        if i % access_every == 0:
                            data[entry[0]]['section_run']['section_system']
            print('archive.py: access every %d-ed entry single entry system (23), blocked %d: ' % (access_every, use_blocked_toc), (time() - start) / 23)

        # just msgpack
        start = time()
        packb(example_archive)
        print('msgpack: create archive (1): ', time() - start)

    benchmark()
