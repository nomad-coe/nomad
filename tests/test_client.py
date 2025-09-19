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
import contextlib
import json
import os

import pytest
from httpx import AsyncClient
from pydantic import ValidationError

from nomad.app.main import app
from nomad.client.archive import ArchiveQuery
from nomad.datamodel import EntryArchive, User
from nomad.datamodel.metainfo import SCHEMA_IMPORT_ERROR, runschema
from nomad.datamodel.metainfo.annotations import Rule, Rules
from nomad.metainfo import MSection, SubSection
from nomad.utils.json_transformer import Transformer
from tests.fixtures.users import users
from tests.processing import test_data as test_processing

# TODO: more tests


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
    results: list[MSection], sub_section_defs: list[SubSection] = None, total=1
):
    assert len(results) == total
    for result in results:
        assert result.m_def == EntryArchive.m_def
        if sub_section_defs:
            current = result
            for sub_section_def in sub_section_defs:
                for other_sub_section_def in current.m_def.all_sub_sections.values():
                    if other_sub_section_def != sub_section_def:
                        assert (
                            len(current.m_get_sub_sections(other_sub_section_def)) == 0
                        )

                sub_sections = current.m_get_sub_sections(sub_section_def)
                assert len(sub_sections) > 0
                current = sub_sections[0]


@pytest.fixture(scope='function')
def many_uploads(non_empty_uploaded: tuple[str, str], user1: User, proc_infra):
    _, upload_file = non_empty_uploaded
    for index in range(0, 4):
        upload = test_processing.run_processing(
            (f'test_upload_{index}', upload_file), user1
        )
        upload.publish_upload()  # pylint: disable=no-member
        try:
            upload.block_until_complete(interval=0.01)
        except Exception:
            pass


@pytest.fixture(scope='module')
def async_api_v1(monkeysession):
    """
    This fixture provides an HTTP client with AsyncClient that accesses
    the fast api. The patch will redirect all requests to the fast api under test.
    """
    test_client = AsyncClient(app=app)

    monkeysession.setattr(
        'nomad.client.archive.ArchiveQuery._fetch_url',
        'http://testserver/api/v1/entries/query',
    )
    monkeysession.setattr(
        'nomad.client.archive.ArchiveQuery._download_url',
        'http://testserver/api/v1/entries/archive/query',
    )

    monkeysession.setattr('httpx.AsyncClient.get', getattr(test_client, 'get'))
    monkeysession.setattr('httpx.AsyncClient.put', getattr(test_client, 'put'))
    monkeysession.setattr('httpx.AsyncClient.post', getattr(test_client, 'post'))
    monkeysession.setattr('httpx.AsyncClient.delete', getattr(test_client, 'delete'))

    def mocked_auth_headers(self) -> dict:
        for user in users.values():
            if user['username'] == self.user or user['email'] == self.user:
                return dict(Authorization=f'Bearer {user["user_id"]}')
        return {}

    monkeysession.setattr('nomad.client.api.Auth.headers', mocked_auth_headers)

    return test_client


def test_async_query_basic(async_api_v1, published_wo_user_metadata):
    async_query = ArchiveQuery()

    assert_results(async_query.download())

    async_query = ArchiveQuery(
        query=dict(upload_id=[published_wo_user_metadata.upload_id])
    )

    assert_results(async_query.download())


@pytest.mark.skipif(runschema is None, reason=SCHEMA_IMPORT_ERROR)
@pytest.mark.parametrize(
    'q_required,sub_sections',
    [
        ({'run': '*'}, [EntryArchive.run]),
        ({'run': {'system': '*'}}, [EntryArchive.run, runschema.run.Run.system]),
        ({'run[0]': {'system': '*'}}, [EntryArchive.run, runschema.run.Run.system]),
    ],
)
def test_async_query_required(
    async_api_v1, published_wo_user_metadata, q_required, sub_sections
):
    async_query = ArchiveQuery(required=q_required)

    assert_results(async_query.download(), sub_section_defs=sub_sections)


def test_async_query_auth(async_api_v1, published, user2, user1):
    async_query = ArchiveQuery(username=user2.username, password='password')

    assert_results(async_query.download(), total=0)

    async_query = ArchiveQuery(username=user1.username, password='password')

    assert_results(async_query.download(), total=1)


def test_async_query_parallel(async_api_v1, many_uploads, monkeypatch):
    async_query = ArchiveQuery(required=dict(run='*'))

    assert_results(async_query.download(), total=4)
    assert_results(async_query.download(), total=0)

    async_query = ArchiveQuery(required=dict(run='*'), page_size=1)

    assert_results(async_query.download(), total=4)


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
