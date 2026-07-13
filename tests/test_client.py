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

import asyncio
import contextlib
import json
import os

import pytest
import pytest_asyncio
from pydantic import ValidationError

from nomad.client.api import APIError, Auth
from nomad.client.archive import ArchiveQuery
from nomad.datamodel import EntryArchive, User
from nomad.datamodel.metainfo import SCHEMA_IMPORT_ERROR, runschema
from nomad.datamodel.metainfo.annotations import Condition, RegexCondition, Rule, Rules
from nomad.metainfo import MSection, SubSection
from nomad.utils.json_transformer import Transformer
from tests.processing import test_data as test_processing

# TODO: more tests


def test_headers_empty_if_no_token():
    auth = Auth(user=None, password=None)
    auth._token = None
    assert auth.headers() == {}


def test_headers_with_token():
    auth = Auth(user='u', password='p')
    auth._token = {'access_token': 'abc'}
    assert auth.headers() == {'Authorization': 'Bearer abc'}


def test_get_access_token_from_api_success(monkeypatch):
    class FakeResponse:
        status_code = 200

        def json(self):
            return {'access_token': 'tok123'}

    monkeypatch.setattr(
        'nomad.client.api.requests.post', lambda *a, **k: FakeResponse()
    )

    auth = Auth(user='u', password='p', from_api=True)
    auth._token = None
    auth.get_access_token_from_api()
    assert auth._token['access_token'] == 'tok123'


def test_get_access_token_from_api_failure(monkeypatch):
    class FakeResponse:
        status_code = 401

        def json(self):
            return {'detail': 'bad creds', 'code': 401}

    monkeypatch.setattr(
        'nomad.client.api.requests.post', lambda *a, **k: FakeResponse()
    )

    auth = Auth(user='u', password='p', from_api=True)
    auth._token = None
    with pytest.raises(APIError) as e:
        auth.get_access_token_from_api()
    assert 'bad creds' in str(e.value)


@pytest.fixture(autouse=True)
def quiet_archivequery_in_ci(monkeypatch):
    """Silence print and progress bar in CI."""
    if os.getenv('CI') == 'true':
        monkeypatch.setattr('builtins.print', lambda *a, **k: None)

        @contextlib.contextmanager
        def dummy_progressbar(*args, **kwargs):
            class DummyBar:
                def __getattr__(self, name):
                    return lambda *a, **k: None

            yield DummyBar()

        from nomad.client import archive

        monkeypatch.setattr(archive, 'progressbar', dummy_progressbar)


def assert_results(
    results: list[MSection], sub_section_defs: list[SubSection] | None = None, total=1
):
    assert len(results) == total
    for result in results:
        assert result.m_def == EntryArchive.m_def
        if sub_section_defs:
            current = result
            for sub_section_def in sub_section_defs:
                assert current.m_def is not None
                for other_sub_section_def in current.m_def.all_sub_sections.values():
                    if other_sub_section_def != sub_section_def:
                        assert (
                            len(current.m_get_sub_sections(other_sub_section_def)) == 0
                        )

                sub_sections = current.m_get_sub_sections(sub_section_def)
                assert len(sub_sections) > 0
                current = sub_sections[0]


@pytest_asyncio.fixture(scope='function')
async def many_uploads(
    non_empty_uploaded: tuple[str, str], user1: User, temporal_worker
):
    _, upload_file = non_empty_uploaded
    async with temporal_worker() as client:
        for index in range(0, 4):
            upload = await asyncio.to_thread(
                test_processing.run_processing,
                (f'test_upload_{index}', upload_file),
                user1,
            )
            await upload._start_publish_upload_workflow()
    yield


@pytest.mark.asyncio
async def test_async_query_basic(
    elastic_function, async_api_v1, published_wo_user_metadata
):
    async_query = ArchiveQuery()

    assert_results(await async_query.async_download())

    async_query = ArchiveQuery(
        query=dict(upload_id=[published_wo_user_metadata.upload_id])
    )

    assert_results(await async_query.async_download())


@pytest.mark.asyncio
@pytest.mark.skipif(runschema is None, reason=SCHEMA_IMPORT_ERROR)
@pytest.mark.parametrize(
    'q_required,sub_sections',
    (
        [
            ({'run': '*'}, [EntryArchive.run]),
            ({'run': {'system': '*'}}, [EntryArchive.run, runschema.run.Run.system]),
            ({'run[0]': {'system': '*'}}, [EntryArchive.run, runschema.run.Run.system]),
        ]
        if runschema is not None
        else []
    ),
)
async def test_async_query_required(
    elastic_function, async_api_v1, published_wo_user_metadata, q_required, sub_sections
):
    async_query = ArchiveQuery(required=q_required)

    assert_results(await async_query.async_download(), sub_section_defs=sub_sections)


@pytest.mark.asyncio
async def test_async_query_auth(
    elastic_function, async_api_v1, published, user2, user1
):
    async_query = ArchiveQuery(username=user2.username, password='password')

    assert_results(await async_query.async_download(), total=0)

    async_query = ArchiveQuery(username=user1.username, password='password')

    assert_results(await async_query.async_download(), total=1)


@pytest.mark.asyncio
async def test_async_query_parallel(
    elastic_function, async_api_v1, many_uploads, monkeypatch
):
    async_query = ArchiveQuery(required=dict(run='*'))

    assert_results(await async_query.async_download(), total=4)
    assert_results(await async_query.async_download(), total=0)

    async_query = ArchiveQuery(required=dict(run='*'), page_size=1)

    assert_results(await async_query.async_download(), total=4)


def load_example(path: str):
    current_dir = os.getcwd()
    rules_path = os.path.join(
        current_dir, 'examples', 'data', 'json_transformer', path + '.json'
    )
    expected_path = os.path.join(
        current_dir, 'examples', 'data', 'json_transformer', 'expected.json'
    )
    with open(rules_path) as file:
        rules_data = json.load(file)
    with open(expected_path) as f:
        expected = json.load(f)

    transformation_dict = {}
    try:
        rules = Rules(**rules_data['schema'])
        transformation_dict[path] = rules
    except ValidationError as ve:
        pytest.fail(f"Validation error in transformation '{path}': {ve}")
    data = rules_data['data']
    expected = expected[path]
    return transformation_dict, data, expected


def load_transformer(transformation_rules):
    """
    Fixture to provide a Transformer instance with predefined rules.
    """
    return Transformer(mapping_dict=transformation_rules)


@pytest.mark.parametrize(
    'transformation_name,target',
    [
        ('basic_transformation', {}),
        ('list_transformation', {}),
        ('dict_to_list_transformation', []),
        ('list_to_list_transformation', []),
        ('dict_to_list_with_none_transformation', []),
        ('dict_with_regex_transformation', {}),
        ('dict_with_regex_transformation_no_match', {}),
        ('dict_with_default_transformation', {}),
        (
            'dict_with_default_transformation_missing_b',
            {},
        ),
        (
            'dict_with_default_transformation_missing_a_and_b',
            {},
        ),
        ('nested_transformation', {}),
        ('complex_transformation', {}),
        ('conditional_transformation_not_met', {}),
        ('conditional_transformation_met', {}),
    ],
)
def test_transform(transformation_name, target):
    """
    General test for Transformer.transform method.
    """

    transformation_rules, source, expected = load_example(transformation_name)
    transformer = load_transformer(transformation_rules)

    assert source is not None, (
        f"'source' key missing in test data for '{transformation_name}'"
    )
    assert target is not None, (
        f"'target' key missing in test data for '{transformation_name}'"
    )
    assert expected is not None, (
        f"'expected' key missing in test data for '{transformation_name}'"
    )

    result = transformer.transform(source, transformation_name, target)
    assert result == expected, f"Failed for transformation '{transformation_name}''"


def test_transform_with_invalid_mapping_name():
    """
    Test that transforming with a non-existent transformation name raises a ValueError.
    """
    transformation_rules, test_data, expected = load_example('basic_transformation')
    transformer = load_transformer(transformation_rules)
    with pytest.raises(ValueError) as exc_info:
        transformer.transform(test_data, 'non_existent_transformation')
    assert (
        "Mapping name 'non_existent_transformation' not found in the transformation dictionary"
        in str(exc_info.value)
    )


def test_transform_with_complex_nested_structure():
    """
    Test transforming a deeply nested structure.
    """
    _, test_data, _ = load_example('nested_transformation')
    mapping = Rules(
        rules={'rule_f_nested_key': Rule(source='f.nested.key', target='result')}
    )
    transformer = load_transformer({'deep_nested': mapping})
    result = transformer.transform(test_data, 'deep_nested', {})
    assert result == {'result': 'value'}, 'Failed for deep_nested transformation'


def test_transform_with_default_null_values():
    """
    Test that setting values in lists with indices beyond current length inserts nulls appropriately.
    """
    target = {'new_list': []}
    Transformer.set_value('new_list[2]', 3, target)
    assert target['new_list'] == [
        None,
        None,
        3,
    ], 'Failed to insert None values correctly'


def test_transform_with_default_value():
    """
    Test that default_value is correctly set when source paths are missing.
    """
    transformation_rules, source_default_with_a, expected_default_with_a = load_example(
        'default_test_transformation_with_a'
    )
    transformer = load_transformer(transformation_rules)
    target_default_with_a = {}
    result = transformer.transform(
        source_default_with_a,
        'default_test_transformation_with_a',
        target_default_with_a,
    )
    assert result == expected_default_with_a, (
        'Failed default_test_transformation_with_a'
    )

    transformation_rules, source_default_without_a, expected_default_without_a = (
        load_example('default_test_transformation_without_a')
    )
    transformer = load_transformer(transformation_rules)
    target_default_without_a = {}
    result = transformer.transform(
        source_default_without_a,
        'default_test_transformation_without_a',
        target_default_without_a,
    )
    assert result == expected_default_without_a, (
        'Failed default_test_transformation_without_a'
    )

    transformation_rules, source_default_without_a_b, expected_default_without_a_b = (
        load_example('default_test_transformation_without_a_and_b')
    )
    transformer = load_transformer(transformation_rules)
    target_default_without_a_b = {}
    result = transformer.transform(
        source_default_without_a_b,
        'default_test_transformation_without_a_and_b',
        target_default_without_a_b,
    )
    assert result == expected_default_without_a_b, (
        'Failed default_test_transformation_without_a_and_b'
    )


def test_transform_with_use_rule_reference():
    """
    A rule with only `use_rule` should inherit source and target from the
    referenced rule and apply the transformation.
    """
    library = Rules(rules={'copy_name': Rule(source='person.name', target='out.name')})
    main = Rules(
        rules={'ref_name': Rule(target='placeholder', use_rule='#lib.copy_name')}
    )
    transformer = load_transformer({'main': main, 'lib': library})

    result = transformer.transform({'person': {'name': 'Ada'}}, 'main', {})
    assert result == {'out': {'name': 'Ada'}}


def test_transform_use_rule_referenced_fields_win():
    """
    `override_fields` lets the referenced rule's non-empty fields overwrite
    the caller's, so the referenced source/target should be used.
    """
    library = Rules(rules={'real': Rule(source='b_src', target='b_out')})
    main = Rules(
        rules={'shadowed': Rule(source='a_src', target='a_out', use_rule='#lib.real')}
    )
    transformer = load_transformer({'main': main, 'lib': library})

    result = transformer.transform({'a_src': 'A', 'b_src': 'B'}, 'main', {})
    assert result == {'b_out': 'B'}


def test_transform_use_rule_propagates_default_value():
    """
    When the referenced rule supplies a default_value and the source path is
    absent, the default is written to the target.
    """
    library = Rules(
        rules={
            'with_default': Rule(
                source='missing.path', target='out.value', default_value=42
            )
        }
    )
    main = Rules(rules={'use_default': Rule(target='x', use_rule='#lib.with_default')})
    transformer = load_transformer({'main': main, 'lib': library})

    assert transformer.transform({}, 'main', {}) == {'out': {'value': 42}}


def test_transform_use_rule_propagates_conditions():
    """
    Conditions defined on the referenced rule must be evaluated and gate
    whether the value is written.
    """
    library = Rules(
        rules={
            'gated': Rule(
                source='payload.value',
                target='out.value',
                conditions=[
                    Condition(
                        regex_condition=RegexCondition(
                            regex_path='payload.kind', regex_pattern=r'^ok$'
                        )
                    )
                ],
            )
        }
    )
    main = Rules(rules={'ref': Rule(target='x', use_rule='#lib.gated')})
    transformer = load_transformer({'main': main, 'lib': library})

    met = transformer.transform({'payload': {'kind': 'ok', 'value': 7}}, 'main', {})
    assert met == {'out': {'value': 7}}

    not_met = transformer.transform(
        {'payload': {'kind': 'nope', 'value': 7}}, 'main', {}
    )
    assert not_met == {}


def test_transform_use_rule_chained_reference():
    """A -> B -> C should resolve through both hops."""
    leaf = Rules(rules={'c': Rule(source='deep.val', target='final.val')})
    mid = Rules(rules={'b': Rule(target='_', use_rule='#leaf.c')})
    top = Rules(rules={'a': Rule(target='_', use_rule='#mid.b')})
    transformer = load_transformer({'top': top, 'mid': mid, 'leaf': leaf})

    result = transformer.transform({'deep': {'val': 'hi'}}, 'top', {})
    assert result == {'final': {'val': 'hi'}}


def test_transform_use_rule_circular_reference_raises():
    a = Rules(rules={'a': Rule(target='ta', use_rule='#b.b')})
    b = Rules(rules={'b': Rule(target='tb', use_rule='#a.a')})
    transformer = load_transformer({'a': a, 'b': b})

    with pytest.raises(ValueError, match='Circular reference'):
        transformer.transform({}, 'a', {})


def test_transform_use_rule_invalid_format_raises():
    main = Rules(rules={'r': Rule(target='t', use_rule='#no_dot_here')})
    transformer = load_transformer({'main': main})
    with pytest.raises(ValueError, match='Invalid use_rule format'):
        transformer.transform({}, 'main', {})


def test_transform_use_rule_unknown_mapping_raises():
    main = Rules(rules={'r': Rule(target='t', use_rule='#missing.x')})
    transformer = load_transformer({'main': main})
    with pytest.raises(ValueError, match="Mapping name 'missing' not found"):
        transformer.transform({}, 'main', {})


def test_transform_use_rule_unknown_rule_name_raises():
    lib = Rules(rules={'real': Rule(source='s', target='t')})
    main = Rules(rules={'r': Rule(target='t', use_rule='#lib.missing')})
    transformer = load_transformer({'main': main, 'lib': lib})
    with pytest.raises(ValueError, match="Rule name 'missing' not found"):
        transformer.transform({}, 'main', {})


def test_transform_with_array_notation():
    """
    Test that array notation in source paths is correctly handled.
    """
    transformation_rules, source, expected = load_example(
        'list_with_array_notation_transformation'
    )
    transformer = load_transformer(transformation_rules)
    result = transformer.transform(
        source, 'list_with_array_notation_transformation', {}, array_rules=True
    )
    assert result == expected, 'Failed list_with_array_notation_transformation'


@pytest.mark.parametrize(
    'source,target',
    [
        pytest.param('a[n]', 'b', id='missing-target-array-notation'),
        pytest.param('a', 'b[n]', id='missing-source-array-notation'),
        pytest.param('a[n].b[n]', 'b[n]', id='number-of-arrays-mismatch'),
        pytest.param('a[n1]', 'b[n2]', id='different-array-indices'),
        pytest.param('a[n1].b[n2]', 'b[n2].c[n1]', id='wrong-array-index-positions'),
    ],
)
def test_transform_with_incorrect_array_notation(source, target):
    """
    Test that incorrect array notation raises an error.
    """
    transformation_rules = {
        'test': Rules(rules={'rule': Rule(source=source, target=target)})
    }
    transformer = load_transformer(transformation_rules)
    with pytest.raises(ValueError):
        transformer.transform({}, 'test', {}, array_rules=True)


@pytest.mark.parametrize(
    'transformation_name,target',
    [
        ('basic_transformation', {}),
        ('list_transformation', {}),
        ('dict_to_list_transformation', []),
        ('list_to_list_transformation', []),
        ('dict_to_list_with_none_transformation', []),
        ('dict_with_regex_transformation', {}),
        ('dict_with_regex_transformation_no_match', {}),
        ('dict_with_default_transformation', {}),
        (
            'dict_with_default_transformation_missing_b',
            {},
        ),
        (
            'dict_with_default_transformation_missing_a_and_b',
            {},
        ),
        ('nested_transformation', {}),
        ('complex_transformation', {}),
        ('conditional_transformation_not_met', {}),
        ('conditional_transformation_met', {}),
    ],
)
def test_transform_use_array_rules_without_array_notation(transformation_name, target):
    """
    Test for Transformer.transform method with array rules without array notation.
    To ensure default behavior is not affected.
    """

    transformation_rules, source, expected = load_example(transformation_name)
    transformer = load_transformer(transformation_rules)

    assert source is not None, (
        f"'source' key missing in test data for '{transformation_name}'"
    )
    assert target is not None, (
        f"'target' key missing in test data for '{transformation_name}'"
    )
    assert expected is not None, (
        f"'expected' key missing in test data for '{transformation_name}'"
    )

    result = transformer.transform(
        source, transformation_name, target, array_rules=True
    )
    assert result == expected, f"Failed for transformation '{transformation_name}''"


def test_transform_array_rules_preseves_other_params():
    transformation_rule = Rule(
        source='a[n]',
        target='b[n]',
        default_value=42,
        conditions=[
            Condition(
                regex_condition=RegexCondition(
                    regex_path='payload.kind', regex_pattern=r'^ok$'
                )
            )
        ],
        use_rule='#lib.some_rule',
    )
    array_rules = Transformer._resolve_array_rule(['a[0]'], transformation_rule, 'test')
    for rule_name, rule in array_rules['test'].rules.items():
        assert rule.default_value == transformation_rule.default_value, (
            'Default value not preserved in array rule resolution'
        )
        assert rule.conditions == transformation_rule.conditions, (
            'Conditions not preserved in array rule resolution'
        )
        assert rule.use_rule == transformation_rule.use_rule, (
            'use_rule not preserved in array rule resolution'
        )


def test_transform_with_array_rules_with_referenced_rule():
    """
    Test that array rules work correctly when referencing another rule.
    """
    library = Rules(rules={'base': Rule(source='a[n]', target='b[n]')})
    main = Rules(rules={'test': Rule(target='x', use_rule='#lib.base')})
    transformer = load_transformer({'main': main, 'lib': library})

    result = transformer.transform({'a': [1, 2]}, 'main', {}, array_rules=True)
    assert result == {'b': [1, 2]}, 'Failed for array rules with referenced rule'


def test_transform_with_inplace_transformation():
    """
    Test that inplace transformation works correctly.
    """
    transformation_rules, source, expected = load_example('inplace_transformation')
    transformer = load_transformer(transformation_rules)
    result = transformer.transform(source, 'inplace_transformation', inplace=True)
    result2 = transformer.transform(source, 'inplace_transformation', source)

    assert result == expected, (
        'Inplace transformation did not modify source as expected'
    )
    assert result2 == result, 'Standard way didnt match inplace transformation result'


def test_transform_with_inplace_deletion_transformation():
    """
    Test that inplace deletion transformation works correctly.
    """
    transformation_rules, source, expected = load_example(
        'inplace_deletion_transformation'
    )
    transformer = load_transformer(transformation_rules)
    result = transformer.transform(
        source, 'inplace_deletion_transformation', inplace=True, delete_sources=True
    )
    result2 = transformer.transform(
        source, 'inplace_deletion_transformation', source, delete_sources=True
    )
    assert result == expected, (
        'Inplace deletion transformation did not modify source as expected'
    )
    assert result2 == result, (
        'Standard way didnt match inplace deletion transformation result'
    )
