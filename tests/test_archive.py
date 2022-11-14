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

from typing import Dict, Any, Union
import pytest
import msgpack
from io import BytesIO
import os.path
import json

from nomad import utils, config
from nomad.metainfo import MSection, Quantity, Reference, SubSection, QuantityReference, MetainfoError, Context
from nomad.datamodel import EntryArchive
from nomad.archive.storage import TOCPacker, _decode, _entries_per_block
from nomad.archive import (
    write_archive, read_archive, ArchiveReader, ArchiveQueryError, query_archive,
    write_partial_archive_to_mongo, read_partial_archive_from_mongo, read_partial_archives_from_mongo,
    create_partial_archive, compute_required_with_referenced, RequiredReader,
    RequiredValidationError)
from nomad.utils.exampledata import ExampleData


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

    toc = _unpack(packed_archive, _decode(archive['toc_pos']))
    assert example_uuid in toc
    assert _unpack(packed_archive, _decode(toc[example_uuid][0])) == toc_packer.toc
    assert _unpack(packed_archive, _decode(toc[example_uuid][1])) == example_entry


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
    archive_size = _entries_per_block * 2 + 23
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
            assert len(reader._toc) == _entries_per_block
            reader._load_toc_block(archive_size - 1)
            assert reader._toc.get(create_example_uuid(archive_size - 1)) is not None
            assert len(reader._toc) > _entries_per_block

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
def json_dict():
    return json.loads(
        '''
{
    "metadata": {
        "entry_id": "test_id"
    },
    "results": {
        "properties": {
            "electronic": {
                "dos_electronic": [{
                    "energies": "/run/0/calculation/1/dos_electronic/0/energies"
                }]
            }
        }
    },
    "run": [
        {
            "system": [
                {
                    "atoms": {
                        "labels": [
                            "He"
                        ]
                    },
                    "symmetry": [
                        {
                            "space_group_number": 221
                        }
                    ]
                },
                {
                    "atoms": {
                        "labels": [
                            "H"
                        ]
                    },
                    "symmetry": [
                        {
                            "space_group_number": 221
                        }
                    ]
                }
            ],
            "calculation": [
                {
                    "system_ref": "/run/0/system/1",
                    "energy": {
                        "total": {
                            "value": 0.1
                        }
                    }
                },
                {
                    "system_ref": "/run/0/system/1",
                    "energy": {
                        "total": {
                            "value": 0.2
                        }
                    },
                    "dos_electronic": [
                        {
                            "energies": [0.0, 0.1]
                        }
                    ],
                    "eigenvalues": [
                    ]
                },
                {
                    "system_ref": "/run/0/system/1",
                    "energy": {
                        "total": {
                            "value": 0.1
                        }
                    }
                }
            ]
        }
    ],
    "workflow": [
        {
            "calculation_result_ref": "/run/0/calculation/1"
        }
    ]
}
''')


@pytest.fixture(scope='session')
def archive(json_dict):
    archive = EntryArchive.m_from_dict(json_dict)
    assert archive.run is not None
    assert len(archive.run) == 1
    return archive


@pytest.mark.parametrize('definition_id,context,exception_type', [
    pytest.param(EntryArchive.m_def.definition_id + 'a', None, MetainfoError, id='wrong_id_no_context'),
    pytest.param(EntryArchive.m_def.definition_id[::-1], Context(), NotImplementedError, id='wrong_id_with_context'), ])
def test_archive_with_wrong_id(json_dict, definition_id, context, exception_type):
    '''
    Test that the archive with wrong id raises the expected exception.
    '''
    json_dict['m_def_id'] = definition_id
    with pytest.raises(exception_type):
        EntryArchive.m_from_dict(json_dict, m_context=context)

    del json_dict['m_def_id']


@pytest.mark.parametrize('m_def,m_def_id', [
    pytest.param(None, EntryArchive.m_def.definition_id, id='plain-definition-id'),
    pytest.param('nomad.datamodel.EntryArchive', None, id='plain-definition-python-style'),
    pytest.param(
        'nomad.datamodel.EntryArchive@' + EntryArchive.m_def.definition_id, None,
        id='plain-definition-with-correct-id'),
    pytest.param(
        'nomad.datamodel.EntryArchive@' + EntryArchive.m_def.definition_id[::-1], None,
        id='plain-definition-with-wrong-id'),
    pytest.param(
        'http://my.domain#/placeholder@' + EntryArchive.m_def.definition_id, None,
        id='url-definition'),
])
def test_archive_with_id_in_reference(json_dict, m_def, m_def_id, monkeypatch):
    '''
    Patch Context to return proper section definition to test if the archive is correctly created.
    '''

    def resolve_section_definition(self, definition: str, definition_id: str):  # pylint: disable=unused-argument
        return EntryArchive

    monkeypatch.setattr('nomad.metainfo.Context.resolve_section_definition', resolve_section_definition)

    if m_def is not None:
        json_dict['m_def'] = m_def
    if m_def_id is not None:
        json_dict['m_def_id'] = m_def_id

    archive = MSection.m_from_dict(json_dict, m_context=Context())
    assert archive.run is not None
    assert len(archive.run) == 1

    if 'm_def' in json_dict:
        del json_dict['m_def']
    if 'm_def_id' in json_dict:
        del json_dict['m_def_id']


@pytest.mark.parametrize('required, error', [
    pytest.param('include', None, id='include-all'),
    pytest.param('*', None, id='include-all-alias'),
    pytest.param({'metadata': '*'}, None, id='include-sub-section'),
    pytest.param({'metadata': {
        'entry_id': '*'
    }}, None, id='include-quantity'),
    pytest.param({
        'workflow': {
            'calculation_result_ref': {
                'energy': {
                    'total': '*'
                }
            }
        }
    }, None, id='resolve-with-required'),
    pytest.param({
        'workflow': {
            'calculation_result_ref': 'include-resolved'
        }
    }, None, id='resolve-with-directive'),
    pytest.param({
        'workflow':
            'include-resolved',
        'results': 'include-resolved'
    }, None, id='include-resolved'),
    pytest.param({
        'results': {
            'properties': {
                'electronic': {
                    'dos_electronic': {
                        'energies': 'include-resolved'
                    }
                }
            }
        }
    }, None, id='resolve-quantity-ref'),
    pytest.param({
        'metadata': {
            'entry_id': {
                'doesnotexist': '*'
            }
        }
    }, ['metadata', 'entry_id'], id='not-a-section'),
    pytest.param({
        'metadata': 'bad-directive'
    }, ['metadata'], id='bad-directive')
])
@pytest.mark.parametrize('resolve_inplace', [
    pytest.param(True, id='inplace'),
    pytest.param(False, id='root'),
])
def test_required_reader(archive, required, error, resolve_inplace):
    f = BytesIO()
    write_archive(f, 1, [('entry_id', archive.m_to_dict())], entry_toc_depth=2)
    packed_archive = f.getbuffer()

    archive_reader = ArchiveReader(BytesIO(packed_archive))
    try:
        required_reader = RequiredReader(required, resolve_inplace=resolve_inplace)
    except RequiredValidationError as e:
        assert error is not None, f'{e.msg}, loc={e.loc}'
        assert e.loc == error
        return

    assert error is None
    results = required_reader.read(archive_reader, 'entry_id', None)

    assert_required_results(results, required_reader.required, archive)


@pytest.fixture(scope='module')
def example_data_with_reference(elastic_module, raw_files_module, mongo_module, test_user, json_dict):
    '''
    Provides a couple of entries with references.

    Only used in test_required_reader_with_remote_reference.
    '''
    data = ExampleData(main_author=test_user)

    data.create_upload(upload_id='id_published_with_ref', upload_name='name_published', published=True)

    ref_list = [
        {'calculation_result_ref': '/run/0/calculation/1'},  # plain direct reference
        {'calculation_result_ref': '#/run/0/calculation/1'},  # new-style reference
        {'workflows_ref': ['../entries/id_01/archive#/workflow/0']},  # reference to another archive
        {'workflows_ref': ['../entries/id_05/archive#/workflow/0']},  # circular reference
        {'workflows_ref': ['../entries/id_04/archive#/workflow/0']},  # circular reference
        {'workflows_ref': ['https://another.domain/entries/id_03/archive#/workflow/0']}  # remote reference
    ]

    del json_dict['results']

    for index, ref in enumerate(ref_list):
        json_dict['workflow'][0] = ref
        data.create_entry(
            upload_id='id_published_with_ref',
            entry_id=f'id_{index + 1:02d}',
            entry_archive=EntryArchive.m_from_dict(json_dict))

    for archive in data.archives.values():
        archive.metadata.apply_archive_metadata(archive)

    data.save(with_files=True, with_es=True, with_mongo=True)

    from nomad.search import search
    results = search().data
    assert len(results) == 6
    for i in {0, 1, 5}:
        assert 'entry_references' not in results[i]
    for i in {2, 3, 4}:
        assert 'entry_references' in results[i]

    yield data
    data.delete()


@pytest.fixture(scope='function')
def remote_reference_required():
    '''
    Only used in test_required_reader_with_remote_reference.
    '''
    return {'workflow': 'include-resolved'}


@pytest.mark.parametrize(
    'resolve_inplace', [
        pytest.param(True, id='inplace'),
        pytest.param(False, id='root'),
    ])
@pytest.mark.parametrize(
    'entry_id, inplace_result', [
        pytest.param(
            'id_01', {'calculation_result_ref': {
                'system_ref': {'atoms': {'labels': ['H']}, 'symmetry': [{'space_group_number': 221}]},
                'energy': {'total': {'value': 0.2}}, 'dos_electronic': [{'energies': [0.0, 0.1]}]}},
            id='plain-direct-reference'),
        pytest.param(
            'id_02', {'calculation_result_ref': {
                'system_ref': {'atoms': {'labels': ['H']}, 'symmetry': [{'space_group_number': 221}]},
                'energy': {'total': {'value': 0.2}}, 'dos_electronic': [{'energies': [0.0, 0.1]}]}},
            id='new-style-reference'),
        pytest.param(
            'id_03', {'calculation_result_ref': {
                'system_ref': {'atoms': {'labels': ['H']}, 'symmetry': [{'space_group_number': 221}]},
                'energy': {'total': {'value': 0.2}}, 'dos_electronic': [{'energies': [0.0, 0.1]}]}},
            id='reference-to-another-archive'),
        pytest.param(
            # circular reference detected thus untouched
            'id_04', '../entries/id_04/archive#/workflow/0',
            id='circular-reference-1'),
        pytest.param(
            # circular reference detected thus untouched
            'id_05', '../entries/id_05/archive#/workflow/0',
            id='circular-reference-2'),
        pytest.param(
            # remote reference detected thus untouched
            'id_06', 'https://another.domain/entries/id_03/archive#/workflow/0',
            id='remote-reference'),
        pytest.param(
            'id_07', '../entries/id_07/archive#/workflow/0',
            id='does-not-exist'),
    ])
def test_required_reader_with_remote_reference(
        json_dict, remote_reference_required, resolve_inplace,
        example_data_with_reference, test_user, entry_id, inplace_result):
    archive = {'workflow': json_dict['workflow']}

    archive['workflow'][0]['workflows_ref'] = [f'../entries/{entry_id}/archive#/workflow/0']

    f = BytesIO()
    write_archive(f, 1, [('entry_id', archive)], entry_toc_depth=2)
    packed_archive = f.getbuffer()

    archive_reader = ArchiveReader(BytesIO(packed_archive))
    required_reader = RequiredReader(
        remote_reference_required, resolve_inplace=resolve_inplace, user=test_user)
    results = required_reader.read(archive_reader, 'entry_id', 'id_published_with_ref')
    ref_result = results['workflow'][0]

    while 'workflows_ref' in ref_result:
        ref_result = ref_result['workflows_ref'][0]

    if resolve_inplace or entry_id == 'id_07':
        assert ref_result == inplace_result
    else:
        # print(results)
        # if not resolved inplace, the target archive is copied to the current archive,
        # so the reference is overwritten by the following pattern,
        # whether the target destination is another reference is not controlled by this reference.
        assert ref_result.endswith(f'/archive/{entry_id}#/workflow/0')

        from nomad.datamodel import ClientContext
        archive_obj = EntryArchive.m_from_dict(results, m_context=ClientContext())
        resolved_obj = archive_obj.workflow[0].workflows_ref[0].m_proxy_resolve()
        from nomad.metainfo import MProxy
        if entry_id in ['id_01', 'id_02']:
            # reference to calculation_result_ref
            assert isinstance(resolved_obj.calculation_result_ref, MProxy)
        else:
            # reference to another workflows_ref
            assert isinstance(resolved_obj.workflows_ref[0], MProxy)


def assert_required_results(
        results: dict, required: dict, archive: MSection,
        current_results: Union[dict, str] = None,
        current_archive_serialized: Union[str, dict] = None):
    '''
    Asserts if the resulting dict from a :class:`RequiredReader` contains everything that
    was requested and if this is consistent with the archive that was read from.
    '''
    # initialize recursion
    if current_archive_serialized is None:
        current_archive_serialized = archive.m_to_dict()
    if current_results is None:
        current_results = results

    definition = required['_def']
    directive = required.get('_directive')

    # assert quantity values
    if isinstance(definition, Quantity):
        if current_results != current_archive_serialized:
            if isinstance(current_archive_serialized, str):
                # assume its a quantity reference
                pass
            else:
                assert False, 'quantity values do not match'
        else:
            return

    # deal with references
    if isinstance(current_archive_serialized, str):
        # The current value must be a reference.
        if directive in ['*', 'include']:
            assert current_results == current_archive_serialized
            return

        if isinstance(current_results, str):
            # It is a reference string, it should be resolveable within an
            # results based archive. We should continue the assert from the resolved
            # results and resolved section in the archive.
            assert current_results == current_archive_serialized
            resolved: Any = archive.m_def.section_cls.m_from_dict(results).m_resolve(current_results)

            if isinstance(resolved, MSection):
                current_results = resolved.m_to_dict()
            else:
                # its a quantity reference
                # assertion only works for np typed quantities with unit
                current_results = list(resolved.m)  # type: ignore

        resolved = archive.m_resolve(current_archive_serialized)
        if isinstance(resolved, MSection):
            current_archive_serialized = resolved.m_to_dict()
        else:
            # its a quantity reference
            # assertion only works for np typed quantities with unit
            assert current_results == list(resolved.m)
            return

    # continue recursion on directives, by extending the required with all possible
    # decends
    def prop_def(prop, section):
        prop_def = section.all_properties[prop]
        if isinstance(prop_def, SubSection):
            return prop_def.sub_section.m_resolved()
        if isinstance(prop_def, Quantity):
            if isinstance(prop_def.type, Reference):
                if isinstance(prop_def.type, QuantityReference):
                    return prop_def.type.target_quantity_def.m_resolved()
                return prop_def.type.target_section_def.m_resolved()

        return prop_def

    if directive is not None:
        # we use a made up required that applies the directive to all values
        required = {
            key: dict(_directive=directive, _prop=key, _def=prop_def(key, definition))
            for key in current_archive_serialized}

    # decend into references and subsections
    for key, value in required.items():
        if key.startswith('_'): continue
        prop = value['_prop']
        assert prop in current_results
        assert prop in current_archive_serialized
        prop_value = current_results[prop]
        prop_definition = value['_def']
        if isinstance(prop_value, list) and not isinstance(prop_definition, Quantity):
            for i, _ in enumerate(prop_value):
                assert_required_results(
                    results, value, archive, prop_value[i], current_archive_serialized[prop][i])
        else:
            assert_required_results(
                results, value, archive, prop_value, current_archive_serialized[prop])


def assert_partial_archive(archive: EntryArchive) -> EntryArchive:
    # test contents
    assert archive.workflow[0].calculation_result_ref is not None
    # test refs
    assert archive.workflow[0].calculation_result_ref.energy.total is not None
    assert len(archive.workflow[0].calculation_result_ref.eigenvalues) == 0
    # test refs of refs
    system = archive.workflow[0].calculation_result_ref.system_ref
    assert system.atoms.labels == ['H']
    assert system.symmetry[0].space_group_number == 221

    return archive


def test_partial_archive(archive):
    partial_archive_dict = create_partial_archive(archive)
    partial_archive = EntryArchive.m_from_dict(partial_archive_dict)
    assert_partial_archive(partial_archive)


def test_partial_archive_read_write(archive, mongo):
    write_partial_archive_to_mongo(archive)
    assert_partial_archive(read_partial_archive_from_mongo('test_id'))


def test_partial_archive_re_write(archive, mongo):
    write_partial_archive_to_mongo(archive)
    archive.metadata.comment = 'changed'
    write_partial_archive_to_mongo(archive)
    archive = assert_partial_archive(read_partial_archive_from_mongo('test_id'))
    assert archive.metadata.comment == 'changed'


def test_read_partial_archives(archive, mongo):
    write_partial_archive_to_mongo(archive)
    assert_partial_archive(read_partial_archives_from_mongo(['test_id'])['test_id'])


def test_compute_required_with_referenced(archive):
    required = compute_required_with_referenced({
        'workflow': {
            'calculation_result_ref': {
                'energy': {
                    'total': '*'
                },
                'system_ref': '*'
            }
        }
    })

    assert required == {
        'workflow': {
            'calculation_result_ref': '*'
        },
        'run': {
            'calculation': {
                'energy': {
                    'total': '*'
                },
                'system_ref': '*'
            },
            'system': '*'
        }
    }


def test_compute_required_incomplete(archive):
    required = compute_required_with_referenced({
        'workflow': {
            'calculation_result_ref': {
                'energy': {
                    'total': '*'
                },
                'dos_electronic': '*'
            }
        }
    })

    assert required is None

    required = compute_required_with_referenced({
        'workflow': {
            'calculation_result_ref': {
                'energy': {
                    'total': '*'
                },
                'system_ref': {
                    'symmetry': '*'
                }
            }
        }
    })

    assert required is not None


def test_compute_required_full():
    assert compute_required_with_referenced('*') is None
