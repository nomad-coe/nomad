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

import pytest

from nomad import utils
from nomad.search import search
from nomad.datamodel import Dataset
from nomad import processing as proc

from tests.utils import ExampleData


logger = utils.get_logger(__name__)


class TestEditRepo():

    def query(self, *uploads):
        return {'upload_id:any': uploads}

    @pytest.fixture(autouse=True)
    def set_api(self, client, elastic, mongo):
        self.api = client

    @pytest.fixture(autouse=True)
    def example_datasets(self, test_user, other_test_user, mongo):
        self.example_dataset = Dataset(
            dataset_id='example_ds', dataset_name='example_ds', user_id=test_user.user_id)
        self.example_dataset.a_mongo.create()

        self.other_example_dataset = Dataset(
            dataset_id='other_example_ds', dataset_name='other_example_ds',
            user_id=other_test_user.user_id)
        self.other_example_dataset.a_mongo.create()

    @pytest.fixture(autouse=True)
    def example_data(self, test_user, other_test_user, raw_files, elastic, mongo):
        # TODO
        example_data = ExampleData()

        example_data.create_upload('upload_1', main_author=test_user, published=True, embargo_length=0)
        example_data.create_entry(upload_id='upload_1')
        example_data.create_upload('upload_2', main_author=test_user, published=True, embargo_length=36)
        example_data.create_entry(upload_id='upload_2')
        example_data.create_entry(upload_id='upload_2')
        example_data.create_upload('upload_3', main_author=other_test_user, published=True, embargo_length=0)
        example_data.create_entry(upload_id='upload_3')

        example_data.save()

    @pytest.fixture(autouse=True)
    def auth(self, test_user_auth):
        self.test_user_auth = test_user_auth

    def perform_edit(self, query=None, verify=False, **kwargs):
        actions = {}
        for key, value in kwargs.items():
            if isinstance(value, list):
                actions[key] = [dict(value=i) for i in value]
            else:
                actions[key] = dict(value=value)

        data = dict(actions=actions)
        if query is not None:
            data.update(query=query)
        if verify:
            data.update(verify=verify)

        return self.api.post('entries/edit', headers=self.test_user_auth, json=data)

    def assert_edit(self, rv, quantity: str, success: bool, message: bool, status_code: int = 200):
        data = rv.json()
        assert rv.status_code == status_code, data
        actions = data.get('actions')
        assert actions is not None
        assert [quantity] == list(actions.keys())
        quantity_actions = actions[quantity]
        if not isinstance(quantity_actions, list):
            quantity_actions = [quantity_actions]
        has_failure = False
        has_message = False
        for action in quantity_actions:
            has_failure = has_failure or not action['success']
            has_message = has_message or ('message' in action)
        assert not has_failure == success
        assert has_message == message

    def mongo(self, *args, edited: bool = True, **kwargs):
        for calc_id in args:
            calc = proc.Calc.objects(calc_id='test_entry_id_%d' % calc_id).first()
            assert calc is not None
            if edited:
                assert calc.last_edit_time is not None
            metadata = calc.mongo_metadata(calc.upload).m_to_dict()
            for key, value in kwargs.items():
                if metadata.get(key) != value:
                    return False
        return True

    def assert_elastic(self, *args, invert: bool = False, **kwargs):
        def assert_entry(get_entries):
            for arg in args:
                entry_id = 'test_entry_id_%d' % arg
                entries = list(get_entries(entry_id))
                assert len(entries) > 0, entry_id
                for entry in entries:
                    for key, value in kwargs.items():
                        if key in ['authors', 'viewers']:
                            ids = [user['user_id'] for user in entry.get(key)]
                            if ids != value:
                                return False
                        else:
                            if entry.get(key) != value:
                                return False

            return True

        # test v1 data
        assert invert != assert_entry(
            lambda id: search(owner=None, query=dict(entry_id=id)).data)

    def test_edit_all_properties(self, test_user, other_test_user):
        edit_data = dict(
            comment='test_edit_props',
            references=['http://test', 'http://test2'],
            # reviewers=[other_test_user.user_id],  # TODO: need to set on upload level
            entry_coauthors=[other_test_user.user_id])
        rv = self.perform_edit(**edit_data, query=self.query('upload_1'))
        result = rv.json()
        assert rv.status_code == 200, result
        actions = result.get('actions')
        for key in edit_data:
            assert key in actions
            quantity_actions = actions.get(key)
            if not isinstance(quantity_actions, list):
                quantity_actions = [quantity_actions]
            for quantity_action in quantity_actions:
                assert quantity_action['success']

        assert self.mongo(1, comment='test_edit_props')
        assert self.mongo(1, references=['http://test', 'http://test2'])
        assert self.mongo(1, entry_coauthors=[other_test_user.user_id])
        # assert self.mongo(1, reviewers=[other_test_user.user_id])  TODO: need to be set on upload level

        self.assert_elastic(1, comment='test_edit_props')
        self.assert_elastic(1, references=['http://test', 'http://test2'])
        self.assert_elastic(1, authors=[test_user.user_id, other_test_user.user_id])
        # self.assert_elastic(1, viewers=[test_user.user_id, other_test_user.user_id])

        edit_data = dict(
            comment='',
            references=[],
            entry_coauthors=[])
        rv = self.perform_edit(**edit_data, query=self.query('upload_1'))
        result = rv.json()
        assert rv.status_code == 200
        actions = result.get('actions')
        for key in edit_data:
            assert key in actions
            quantity_actions = actions.get(key)
            if not isinstance(quantity_actions, list):
                quantity_actions = [quantity_actions]
            for quantity_action in quantity_actions:
                assert quantity_action['success']

        assert self.mongo(1, comment=None)
        assert self.mongo(1, references=[])
        assert self.mongo(1, entry_coauthors=[])
        assert self.mongo(1, reviewers=None)

        self.assert_elastic(1, comment=None)
        self.assert_elastic(1, references=[])
        self.assert_elastic(1, authors=[test_user.user_id])
        self.assert_elastic(1, viewers=[test_user.user_id])

    def test_edit_all(self):
        rv = self.perform_edit(comment='test_edit_all')
        self.assert_edit(rv, quantity='comment', success=True, message=False)
        assert self.mongo(1, 2, 3, comment='test_edit_all')
        self.assert_elastic(1, 2, 3, comment='test_edit_all')
        assert not self.mongo(4, comment='test_edit_all', edited=False)
        self.assert_elastic(4, comment='test_edit_all', edited=False, invert=True)

    def test_edit_multi(self):
        rv = self.perform_edit(comment='test_edit_multi', query=self.query('upload_1', 'upload_2'))
        self.assert_edit(rv, quantity='comment', success=True, message=False)
        assert self.mongo(1, 2, 3, comment='test_edit_multi')
        self.assert_elastic(1, 2, 3, comment='test_edit_multi')
        assert not self.mongo(4, comment='test_edit_multi', edited=False)
        self.assert_elastic(4, comment='test_edit_multi', edited=False, invert=True)

    def test_edit_some(self):
        rv = self.perform_edit(comment='test_edit_some', query=self.query('upload_1'))
        self.assert_edit(rv, quantity='comment', success=True, message=False)
        assert self.mongo(1, comment='test_edit_some')
        self.assert_elastic(1, comment='test_edit_some')
        assert not self.mongo(2, 3, 4, comment='test_edit_some', edited=False)
        self.assert_elastic(2, 3, 4, comment='test_edit_some', edited=False, invert=True)

    def test_edit_verify(self):
        rv = self.perform_edit(
            comment='test_edit_verify', verify=True, query=self.query('upload_1'))
        self.assert_edit(rv, quantity='comment', success=True, message=False)
        assert not self.mongo(1, comment='test_edit_verify', edited=False)

    def test_edit_empty_list(self, other_test_user):
        rv = self.perform_edit(entry_coauthors=[other_test_user.user_id], query=self.query('upload_1'))
        self.assert_edit(rv, quantity='entry_coauthors', success=True, message=False)
        rv = self.perform_edit(entry_coauthors=[], query=self.query('upload_1'))
        self.assert_edit(rv, quantity='entry_coauthors', success=True, message=False)
        assert self.mongo(1, entry_coauthors=[])

    def test_edit_duplicate_value(self, other_test_user):
        rv = self.perform_edit(entry_coauthors=[other_test_user.user_id, other_test_user.user_id], query=self.query('upload_1'))
        self.assert_edit(rv, status_code=400, quantity='entry_coauthors', success=False, message=True)

    def test_edit_main_author_as_coauthor(self, test_user):
        rv = self.perform_edit(entry_coauthors=[test_user.user_id], query=self.query('upload_1'))
        self.assert_edit(rv, status_code=400, quantity='entry_coauthors', success=False, message=True)

    def test_edit_ds(self):
        rv = self.perform_edit(
            datasets=[self.example_dataset.dataset_name], query=self.query('upload_1'))
        self.assert_edit(rv, quantity='datasets', success=True, message=False)
        assert self.mongo(1, datasets=[self.example_dataset.dataset_id])

    def test_edit_ds_remove_doi(self):
        rv = self.perform_edit(
            datasets=[self.example_dataset.dataset_name], query=self.query('upload_1'))

        assert rv.status_code == 200
        rv = self.api.post('datasets/%s/action/doi' % self.example_dataset.dataset_name, headers=self.test_user_auth)
        assert rv.status_code == 200
        rv = self.perform_edit(datasets=[], query=self.query('upload_1'))
        assert rv.status_code == 400
        data = rv.json()
        assert not data['success']
        assert self.example_dataset.dataset_name in data['message']
        assert Dataset.m_def.a_mongo.get(dataset_id=self.example_dataset.dataset_id) is not None

    def test_edit_ds_remove(self):
        rv = self.perform_edit(
            datasets=[self.example_dataset.dataset_name], query=self.query('upload_1'))
        assert rv.status_code == 200
        rv = self.perform_edit(datasets=[], query=self.query('upload_1'))
        assert rv.status_code == 200
        with pytest.raises(KeyError):
            assert Dataset.m_def.a_mongo.get(dataset_id=self.example_dataset.dataset_id) is None

    def test_edit_ds_user_namespace(self, test_user):
        assert Dataset.m_def.a_mongo.objects(
            dataset_name=self.other_example_dataset.dataset_name).first() is not None

        rv = self.perform_edit(
            datasets=[self.other_example_dataset.dataset_name], query=self.query('upload_1'))

        self.assert_edit(rv, quantity='datasets', success=True, message=True)
        new_dataset = Dataset.m_def.a_mongo.objects(
            dataset_name=self.other_example_dataset.dataset_name,
            user_id=test_user.user_id).first()
        assert new_dataset is not None
        assert self.mongo(1, datasets=[new_dataset.dataset_id])

    def test_edit_new_ds(self, test_user):
        rv = self.perform_edit(datasets=['new_dataset'], query=self.query('upload_1'))
        self.assert_edit(rv, quantity='datasets', success=True, message=True)
        new_dataset = Dataset.m_def.a_mongo.objects(dataset_name='new_dataset').first()
        assert new_dataset is not None
        assert new_dataset.user_id == test_user.user_id
        assert self.mongo(1, datasets=[new_dataset.dataset_id])

    def test_edit_bad_user(self):
        rv = self.perform_edit(entry_coauthors=['bad_user'], query=self.query('upload_1'))
        self.assert_edit(rv, status_code=400, quantity='entry_coauthors', success=False, message=True)

    def test_edit_user(self, other_test_user):
        rv = self.perform_edit(entry_coauthors=[other_test_user.user_id], query=self.query('upload_1'))
        self.assert_edit(rv, quantity='entry_coauthors', success=True, message=False)

    @pytest.mark.skip(reason='Not necessary during transition. Fails because main_author is not editable anyways.')
    def test_admin_only(self, other_test_user):
        rv = self.perform_edit(main_author=other_test_user.user_id)
        assert rv.status_code != 200
