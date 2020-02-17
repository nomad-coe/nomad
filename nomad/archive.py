from typing import Iterable, Any, Tuple, Dict, BinaryIO, Union, List, cast
from io import BytesIO
from collections.abc import Mapping, Sequence
import msgpack
from msgpack.fallback import Packer, StringIO
import struct
import json

from nomad import utils


__packer = msgpack.Packer(autoreset=True, use_bin_type=True)


def packb(o, **kwargs):
    return __packer.pack(o)


def unpackb(o, **kwargs):
    return msgpack.unpackb(o, raw=False)


class TOCPacker(Packer):
    """
    A special msgpack packer that records a TOC while packing.

    Uses a combination of the pure python msgpack fallback packer and the "real"
    c-based packing.
    """
    def __init__(self, toc_depth: int, *args, **kwargs):
        self.toc_depth = toc_depth
        self.toc: Dict[str, Any] = None
        self._depth = 0

        # Because we cannot change msgpacks interface of _pack, this _stack is used to
        # tranfer the result of _pack calls in terms of the TOC.
        self._stack: List[Any] = []

        super().__init__(*args, **kwargs)

    def pack(self, obj, *args, **kwargs):
        assert isinstance(obj, dict), 'TOC packer can only pack dicts'
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
                toc_result['toc'] = toc

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
        toc_item_length = len(packb(utils.create_uuid()))
        toc_item_length += len(packb([
            ArchiveWriter._encode_position(0, 0), ArchiveWriter._encode_position(0, 0)]))

        self._write_map_header(3)
        self.write(packb('toc_pos'))
        self.write(packb(ArchiveWriter._encode_position(0, 0)))

        self.write(packb('toc'))
        toc_start, _ = self._write_map_header(self.n_entries)
        _, toc_end = self.write(b'0' * toc_item_length * self.n_entries)
        self._toc_position = toc_start, toc_end

        self.write(packb('data'))
        self._write_map_header(self.n_entries)

        return self

    def write(self, b: bytes) -> Tuple[int, int]:
        start = self._pos
        self._pos += self._f.write(b)
        return start, self._pos

    def __exit__(self, exc_type, exc_val, exc_tb):
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
        assert toc_position == self._toc_position

        if isinstance(self.file_or_path, str):
            self._f.close()

    @staticmethod
    def _encode_position(start: int, end: int) -> bytes:
        return start.to_bytes(5, byteorder='little', signed=False) + \
            end.to_bytes(5, byteorder='little', signed=False)

    def add(self, uuid: str, data: Any) -> None:
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
    def __init__(self, file_or_path: Union[str, BytesIO]):
        self.file_or_path = file_or_path

        f: BytesIO = None
        if isinstance(self.file_or_path, str):
            f = cast(BytesIO, open(self.file_or_path, 'rb'))
        elif isinstance(self.file_or_path, BytesIO):
            f = self.file_or_path
        else:
            raise ValueError('not a file or path')

        super().__init__(None, f)

        # TODO do not load the whole top-level TOC. It has a fixed layout based on
        # the msgpack spec and can be loaded block-by block. The calc_ids are uniformly
        # distributed and sorted!
        # this number is determined by the msgpack encoding of the file beginning:
        # { 'toc_pos': <...>
        #              ^11
        self._f.seek(11)
        toc_position = ArchiveReader._decode_position(self._f.read(10))
        self.toc_entry = self._read(toc_position)

    def __enter__(self):
        return self

    def __getitem__(self, key):
        toc_position = ArchiveReader._decode_position(self.toc_entry[key][0])
        data_position = ArchiveReader._decode_position(self.toc_entry[key][1])
        toc = self._read(toc_position)
        return ArchiveObject(toc, self._f, data_position[0])

    def __iter__(self):
        return self.toc_entry.__iter__()

    def __len__(self):
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


def write_archive(
        path_or_file: Union[str, BytesIO], n_entries: int, data: Iterable[Tuple[str, Any]],
        entry_toc_depth: int = 2) -> None:
    """
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
    """
    with ArchiveWriter(path_or_file, n_entries, entry_toc_depth=entry_toc_depth) as writer:
        for uuid, entry in data:
            writer.add(uuid, entry)


def read_archive(file_or_path: str) -> ArchiveReader:
    """
    Allows to read a msgpack-based archive.

    Arguments:
        file_or_path: A file path or file-like to the archive file that should be read. The
            respective file has to be closed by the user. The returned obj supports the
            'with' statement and has a 'close' method.

    Returns:
        A mapping (dict-like) that can be used to access the archive data. The mapping
        will lazyly load data as it is used. The mapping needs to be closed or used within
        a 'with' statement to free the underlying file resource after use.
    """

    return ArchiveReader(file_or_path)


class PositionEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, bytes):
            return 'position'

        return json.JSONEncoder.default(self, obj)


if __name__ == '__main__':

    def benchmark():
        from time import time

        with open('local/test_be.json') as f:
            example_data = json.load(f)

        size = 1000
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
        start = time()
        for _ in range(0, 23):
            read_archive(buffer)[example_uuid]['section_run']['section_system']
        print('archive.py: access single entry system (23): ', (time() - start) / 23)

        # read every n-ed entry from archive
        buffer = BytesIO(buffer.getbuffer())
        start = time()
        for _ in range(0, 23):
            with read_archive(buffer) as data:
                for i, entry in enumerate(example_archive):
                    if i % access_every == 0:
                        data[entry[0]]['section_run']['section_system']
        print('archive.py: access every %d-ed entry single entry system (23): ' % access_every, (time() - start) / 23)

        # just msgpack
        start = time()
        packb(example_archive)
        print('msgpack: create archive (1): ', time() - start)

        # v0.8.0 impl
        from nomad.archive_library import filedb
        start = time()
        buffer = BytesIO()
        db = filedb.ArchiveFileDB(buffer, mode='w', max_lfragment=3)
        db.add_data({
            uuid: data for uuid, data in example_archive})
        db.close(save=False)
        print('filedb.py: create archive (1): ', time() - start)

        buffer = BytesIO(buffer.getbuffer())
        start = time()
        for _ in range(0, 23):
            db = filedb.ArchiveFileDB(buffer, mode='r', max_lfragment=3)
            db.get_docs(db.ids[example_uuid + '/section_run/section_system'][0])
        print('filedb.py: access single entry system (23): ', (time() - start) / 23)

        buffer = BytesIO(buffer.getbuffer())
        start = time()
        db = filedb.ArchiveFileDB(buffer, mode='r', max_lfragment=3)
        for _ in range(0, 23):
            for i, entry in enumerate(example_archive):
                if i % access_every == 0:
                    db.get_docs(db.ids[entry[0] + '/section_run/section_system'][0])
        print('filedb.py: access every %d-ed entry single entry system (23): ' % access_every, (time() - start) / 23)

    benchmark()
