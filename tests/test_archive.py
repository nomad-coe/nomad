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

from typing import Dict, Any
import pytest
import msgpack
from io import BytesIO
import os.path
import json

from nomad import utils, config
from nomad.datamodel import EntryArchive
from nomad.archive import (
    TOCPacker, write_archive, read_archive, ArchiveReader, ArchiveQueryError, query_archive,
    write_partial_archive_to_mongo, read_partial_archive_from_mongo, read_partial_archives_from_mongo,
    create_partial_archive, compute_required_with_referenced)


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


test_query_example: Dict[Any, Any] = {
    'c1': {
        's1': {
            'ss1': [{'p1': 1.0, 'p2': 'x'}, {'p1': 1.5, 'p2': 'y'}]
        },
        's2': [{'p1': ['a', 'b'], 'p2': True}]
    },
    'c2': {
        's1': {
            'ss1': [{'p1': 2.0}]
        },
        's2': [{'p1': ['c', 'd']}]
    }
}


@pytest.mark.parametrize('query,ref', [
    ({'c1': '*'}, {'c1': test_query_example['c1']}),
    ({'c1': '*', 'c2': {'s1': '*'}}, {'c1': test_query_example['c1'], 'c2': {'s1': test_query_example['c2']['s1']}}),
    ({'c2': {'s1': {'ss1[0]': '*'}}}, {'c2': {'s1': {'ss1': test_query_example['c2']['s1']['ss1'][0:1]}}}),
    ({'c1': {'s1': {'ss1[1:]': '*'}}}, {'c1': {'s1': {'ss1': test_query_example['c1']['s1']['ss1'][1:]}}}),
    ({'c1': {'s1': {'ss1[:2]': '*'}}}, {'c1': {'s1': {'ss1': test_query_example['c1']['s1']['ss1'][:2]}}}),
    ({'c1': {'s1': {'ss1[0:2]': '*'}}}, {'c1': {'s1': {'ss1': test_query_example['c1']['s1']['ss1'][0:2]}}}),
    ({'c1': {'s1': {'ss1[-2]': '*'}}}, {'c1': {'s1': {'ss1': test_query_example['c1']['s1']['ss1'][-2:-1]}}}),
    ({'c1': {'s1': {'ss1[-10]': '*'}}}, {'c1': {'s1': {'ss1': test_query_example['c1']['s1']['ss1'][-2:-1]}}}),
    ({'c1': {'s1': {'ss1[:-1]': '*'}}}, {'c1': {'s1': {'ss1': test_query_example['c1']['s1']['ss1'][:-1]}}}),
    ({'c1': {'s1': {'ss1[1:-1]': '*'}}}, {'c1': {'s1': {'ss1': test_query_example['c1']['s1']['ss1'][1:-1]}}}),
    ({'c2': {'s1': {'ss1[-3:-1]': '*'}}}, {'c2': {'s1': {'ss1': [test_query_example['c2']['s1']['ss1'][-1]]}}}),
    ({'c1': {'s2[0]': {'p1': '*'}}}, {'c1': {'s2': [{'p1': test_query_example['c1']['s2'][0]['p1']}]}}),
    ({'c1': {'s3': '*'}}, {'c1': {}}),
    ({'c1': {'s1[0]': '*'}}, ArchiveQueryError())
])
def test_query(query, ref):
    f = BytesIO()
    write_archive(f, 2, [(k, v) for k, v in test_query_example.items()], entry_toc_depth=1)
    packed_archive = f.getbuffer()

    f = BytesIO(packed_archive)
    if isinstance(ref, Exception):
        with pytest.raises(ref.__class__):
            query_archive(f, query)

    else:
        assert query_archive(f, query) == ref


@pytest.mark.parametrize('key', ['simple', '  fixedsize', 'z6qp-VxV5uacug_1xTBhm5xxU2yZ'])
def test_keys(key):
    f = BytesIO()
    write_archive(f, 1, [(key, dict(example='content'))])
    packed_archive = f.getbuffer()
    f = BytesIO(packed_archive)
    assert key.strip() in query_archive(f, {key: '*'})


def test_read_springer():
    springer = read_archive(config.normalize.springer_db_path)
    with pytest.raises(KeyError):
        springer['doesnotexist']


@pytest.fixture(scope='session')
def archive():
    return EntryArchive.m_from_dict(json.loads('''
        {
            "section_metadata": {
                "calc_id": "test_id",
                "encyclopedia": {
                    "properties": {
                        "electronic_dos": "/section_run/0/section_single_configuration_calculation/1/section_dos/0"
                    }
                }
            },
            "section_run": [
                {
                    "section_single_configuration_calculation": [
                        {
                            "energy_total": 0.1
                        },
                        {
                            "energy_total": 0.2,
                            "section_dos": [
                                {
                                    "dos_kind": "test"
                                }
                            ],
                            "section_eigenvalues": [
                                {
                                    "eigenvalues_kind": "test"
                                }
                            ],
                            "single_configuration_calculation_to_system_ref": "/section_run/0/section_system/1"
                        },
                        {
                            "energy_total": 0.1
                        }
                    ],
                    "section_system": [
                        {
                            "atom_labels": ["He"]
                        },
                        {
                            "atom_labels": ["H"],
                            "section_symmetry": [
                                {
                                    "space_group_number": 221
                                }
                            ]
                        }
                    ]
                }
            ],
            "section_workflow": {
                "calculation_result_ref": "/section_run/0/section_single_configuration_calculation/1"
            }
        }
        '''))


def assert_partial_archive(archive: EntryArchive) -> EntryArchive:
    # test contents
    assert archive.section_workflow.calculation_result_ref is not None
    assert archive.section_metadata.encyclopedia is not None
    # test refs
    assert archive.section_workflow.calculation_result_ref.energy_total is not None
    assert len(archive.section_workflow.calculation_result_ref.section_eigenvalues) == 0
    # test refs of refs
    assert archive.section_workflow.calculation_result_ref.single_configuration_calculation_to_system_ref.atom_labels == ['H']
    assert archive.section_workflow.calculation_result_ref.single_configuration_calculation_to_system_ref.section_symmetry[0].space_group_number == 221

    return archive


def test_partial_archive(archive):
    partial_archive_dict = create_partial_archive(archive)
    partial_archive = EntryArchive.m_from_dict(partial_archive_dict)
    assert_partial_archive(partial_archive)


def test_parital_archive_read_write(archive, mongo):
    write_partial_archive_to_mongo(archive)
    assert_partial_archive(read_partial_archive_from_mongo('test_id'))


def test_partial_archive_re_write(archive, mongo):
    write_partial_archive_to_mongo(archive)
    archive.section_metadata.comment = 'changed'
    write_partial_archive_to_mongo(archive)
    archive = assert_partial_archive(read_partial_archive_from_mongo('test_id'))
    assert archive.section_metadata.comment == 'changed'


def test_read_partial_archives(archive, mongo):
    write_partial_archive_to_mongo(archive)
    assert_partial_archive(read_partial_archives_from_mongo(['test_id'])['test_id'])


def test_compute_required_with_referenced(archive):
    required = compute_required_with_referenced({
        'section_workflow': {
            'calculation_result_ref': {
                'energy_total': '*',
                'single_configuration_calculation_to_system_ref': '*'
            }
        }
    })

    assert required == {
        'section_workflow': {
            'calculation_result_ref': '*'
        },
        'section_run': {
            'section_single_configuration_calculation': {
                'energy_total': '*',
                'single_configuration_calculation_to_system_ref': '*'
            },
            'section_system': '*'
        }
    }


def test_compute_required_incomplete(archive):
    required = compute_required_with_referenced({
        'section_workflow': {
            'calculation_result_ref': {
                'energy_total': '*',
                'section_dos': '*'
            }
        }
    })

    assert required is None

    required = compute_required_with_referenced({
        'section_workflow': {
            'calculation_result_ref': {
                'energy_total': '*',
                'single_configuration_calculation_to_system_ref': {
                    'section_symmetry': '*'
                }
            }
        }
    })

    assert required is not None
