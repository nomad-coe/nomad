
import pytest
import msgpack
from io import BytesIO
import os.path

from nomad import utils, config
from nomad.archive import TOCPacker, write_archive, read_archive, ArchiveReader, query_archive


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


def test_write_archive_empty():
    f = BytesIO()
    write_archive(f, 0, [])


def test_short_uuids():
    f = BytesIO()
    write_archive(f, 1, [('0', {'archive': 'test'})])

    packed_archive = f.getbuffer()
    f = BytesIO(packed_archive)
    with read_archive(f) as archive:
        assert '0' in archive
        assert archive['0'].to_dict() == {'archive': 'test'}


def test_write_file(raw_files, example_uuid):
    path = os.path.join(config.fs.tmp, 'test.msg')
    write_archive(path, 1, [(example_uuid, {'archive': 'test'})])
    with read_archive(path) as archive:
        assert example_uuid in archive
        assert archive[example_uuid].to_dict() == {'archive': 'test'}


def test_write_archive_single(example_uuid, example_entry):
    f = BytesIO()
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

    with pytest.raises(KeyError):
        data['does not exist']

    with pytest.raises(KeyError):
        data[example_uuid]['does not exist']

    with pytest.raises(IndexError):
        data[example_uuid]['run']['system'][2]


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


def test_query():
    payload = {
        'c1': {
            's1': {
                'ss1': [{'p1': 1.0, 'p2': 'x'}, {'p1': 1.5, 'p2': 'y'}]
            },
            's2': {'p1': ['a', 'b']}
        },
        'c2': {
            's1': {'ss1': [{'p1': 2.0}]},
            's2': {'p1': ['c', 'd']}
        }
    }

    f = BytesIO()
    write_archive(f, 2, [(k, v) for k, v in payload.items()], entry_toc_depth=1)
    packed_archive = f.getbuffer()

    f = BytesIO(packed_archive)
    assert query_archive(f, {'c1': '*'}) == {'c1': payload['c1']}
    assert query_archive(f, {'c1': '*', 'c2': {'s1': '*'}}) == {'c1': payload['c1'], 'c2': {'s1': payload['c2']['s1']}}
    assert query_archive(f, {'c2': {'s1': {'ss1[0]': '*'}}}) == {'c2': {'s1': {'ss1': payload['c2']['s1']['ss1'][0]}}}
    assert query_archive(f, {'c1': {'s1': {'ss1[1:]': '*'}}}) == {'c1': {'s1': {'ss1': payload['c1']['s1']['ss1'][1:]}}}
    assert query_archive(f, {'c1': {'s1': {'ss1[:2]': '*'}}}) == {'c1': {'s1': {'ss1': payload['c1']['s1']['ss1'][:2]}}}
    assert query_archive(f, {'c1': {'s1': {'ss1[0:2]': '*'}}}) == {'c1': {'s1': {'ss1': payload['c1']['s1']['ss1'][0:2]}}}
    assert query_archive(f, {'c1': {'s1': {'ss1[-2]': '*'}}}) == {'c1': {'s1': {'ss1': payload['c1']['s1']['ss1'][-2]}}}
    assert query_archive(f, {'c1': {'s1': {'ss1[:-1]': '*'}}}) == {'c1': {'s1': {'ss1': payload['c1']['s1']['ss1'][:-1]}}}
    assert query_archive(f, {'c1': {'s1': {'ss1[1:-1]': '*'}}}) == {'c1': {'s1': {'ss1': payload['c1']['s1']['ss1'][1:-1]}}}
    assert query_archive(f, {'c2': {'s1': {'ss1[-3:-1]': '*'}}}) == {'c2': {'s1': {'ss1': payload['c2']['s1']['ss1'][-3:-1]}}}


def test_read_springer():
    springer = read_archive(config.normalize.springer_db_path)
    with pytest.raises(KeyError):
        springer['doesnotexist']


if __name__ == '__main__':
    import sys
    import pprint
    with open(sys.argv[1], 'rb') as f:
        data = msgpack.unpack(f)
        pprint.pprint(data)
