
import pytest
import msgpack
from io import BytesIO

from nomad import utils
from nomad.archive import TOCPacker, write_archive, read_archive, ArchiveReader, ArchiveFileDB


def create_example_uuid(index: int = 0):
    return ('{:%dd}' % utils.default_hash_len).format(index)


@pytest.fixture(scope='session')
def example_uuid():
    return create_example_uuid()


@pytest.fixture(scope='session')
def example_entry():
    return {
        'run': {
            'program_name': 'VASP',
            'method': {
                'basis_sets': 'plane waves'
            },
            'system': [
                {
                    'atom_labels': ['H', 'H'],
                },
                {
                    'atom_labels': ['H', 'H'],
                    'symmetry': {
                        'space_group': 4
                    }
                },
            ]
        },
        'repo_entry': {
            'chemical_formula': 'H2'
        }
    }


def _unpack(data, pos=None):
    f = BytesIO(data)
    if pos is None:
        return msgpack.unpackb(f.read(), raw=False)
    else:
        f.seek(pos[0])
        return msgpack.unpackb(f.read(pos[1] - pos[0]), raw=False)


def test_toc_packer(example_entry):
    toc_packer = TOCPacker(toc_depth=2)

    toc_packer.reset()
    data = toc_packer.pack(example_entry)
    toc = toc_packer.toc

    assert toc is not None
    assert 'pos' in toc
    assert _unpack(data, toc['pos']) == example_entry

    assert 'run' in toc['toc']
    toc = toc['toc']['run']
    assert _unpack(data, toc['pos']) == example_entry['run']

    assert 'program_name' not in toc
    assert 'system' in toc['toc']
    toc = toc['toc']['system']
    assert isinstance(toc, list)
    assert 'pos' in toc[0]
    assert 'toc' not in toc[0]
    assert _unpack(data, toc[0]['pos']) == example_entry['run']['system'][0]

    assert data is not None
    assert msgpack.unpackb(data, raw=False) == example_entry


def test_write_archive_single(example_uuid, example_entry):
    f = BytesIO()
    print('#', len(example_uuid))
    write_archive(f, 1, [(example_uuid, example_entry)])
    packed_archive = f.getbuffer()
    archive = _unpack(packed_archive)

    assert 'toc_pos' in archive
    assert 'toc' in archive
    assert 'data' in archive
    assert example_uuid in archive['data']
    assert 'data' in archive['data'][example_uuid]
    assert archive['data'][example_uuid]['data'] == example_entry

    toc_packer = TOCPacker(toc_depth=2)
    toc_packer.reset()
    toc_packer.pack(example_entry)
    assert archive['data'][example_uuid]['toc'] == toc_packer.toc

    toc = _unpack(packed_archive, ArchiveReader._decode_position(archive['toc_pos']))
    assert example_uuid in toc
    assert _unpack(
        packed_archive, ArchiveReader._decode_position(toc[example_uuid][0])) == toc_packer.toc
    assert _unpack(
        packed_archive, ArchiveReader._decode_position(toc[example_uuid][1])) == example_entry


def test_write_archive_multi(example_uuid, example_entry):
    f = BytesIO()
    example_uuids = create_example_uuid(0), create_example_uuid(1)
    write_archive(f, 2, [
        (example_uuids[0], example_entry),
        (example_uuids[1], example_entry)])
    packed_archive = f.getbuffer()
    archive = _unpack(packed_archive)

    example_uuid = example_uuids[1]
    assert 'toc_pos' in archive
    assert 'toc' in archive
    assert 'data' in archive
    assert example_uuid in archive['data']
    assert 'data' in archive['data'][example_uuid]
    assert archive['data'][example_uuid]['data'] == example_entry

    toc = archive['toc']
    assert len(toc) == 2
    assert example_uuid in toc


@pytest.mark.parametrize('use_blocked_toc', [False, True])
def test_read_archive_single(example_uuid, example_entry, use_blocked_toc):
    f = BytesIO()
    write_archive(f, 1, [(example_uuid, example_entry)])
    packed_archive = f.getbuffer()

    f = BytesIO(packed_archive)
    data = read_archive(f, use_blocked_toc=use_blocked_toc)

    assert example_uuid in data
    assert data[example_uuid]['run']['system'][1] == example_entry['run']['system'][1]
    assert data[example_uuid]['run'].to_dict() == example_entry['run']
    assert data[example_uuid].to_dict() == example_entry


@pytest.mark.parametrize('use_blocked_toc', [False, True])
def test_read_archive_multi(example_uuid, example_entry, use_blocked_toc):
    archive_size = ArchiveReader.toc_block_size_entries * 2 + 23
    f = BytesIO()
    write_archive(
        f, archive_size,
        [(create_example_uuid(i), example_entry) for i in range(0, archive_size)])
    packed_archive = f.getbuffer()

    f = BytesIO(packed_archive)
    with ArchiveReader(f, use_blocked_toc=use_blocked_toc) as reader:
        if use_blocked_toc:
            reader._load_toc_block(0)
            assert reader._toc.get(create_example_uuid(0)) is not None
            assert len(reader._toc) == ArchiveReader.toc_block_size_entries
            reader._load_toc_block(archive_size - 1)
            assert reader._toc.get(create_example_uuid(archive_size - 1)) is not None
            assert len(reader._toc) > ArchiveReader.toc_block_size_entries

        for i in range(0, archive_size):
            reader.get(create_example_uuid(i)) is not None


def test_archivefiledb():
    payload = [
        {'calc1': {
            'secA': {'subsecA1': [{'propA1a': 1.0}]}, 'secB': {'propB1a': ['a', 'b']}}},
        {'calc2': {
            'secA': {'subsecA1': [{'propA1a': 2.0}]}, 'secB': {'propB1a': ['c', 'd']}}}]

    f = BytesIO()
    msgdb = ArchiveFileDB(f, mode='w', entry_toc_depth=1)
    msgdb.add_data(payload)
    msgdb.close(save=False)

    f = BytesIO(f.getbuffer())
    msgdb = ArchiveFileDB(f, mode='r')

    assert msgdb.query({'calc1': '*'}) == payload[0]
    assert msgdb.query({'calc2': {'secA': {'subsecA1[0]': '*'}}}) == {'calc2': {'secA': {'subsecA1[0]': {'propA1a': 2.0}}}}

    msgdb.close()
