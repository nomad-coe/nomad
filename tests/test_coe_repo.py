# Copyright 2018 Markus Scheidgen
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

import pytest
import datetime

from nomad.coe_repo import User, Calc, Upload

from tests.processing.test_data import processed_upload  # pylint: disable=unused-import
from tests.processing.test_data import uploaded_id  # pylint: disable=unused-import
from tests.processing.test_data import mocks_forall  # pylint: disable=unused-import
from tests.test_files import clear_files  # pylint: disable=unused-import


def assert_user(user, reference):
    assert user is not None
    assert user.user_id == reference.user_id
    assert user.email == reference.email


def test_token_authorize(test_user):
    user = User.verify_auth_token(test_user.email)
    assert_user(user, test_user)


def test_password_authorize(test_user):
    user = User.verify_user_password(test_user.email, 'password')
    assert_user(user, test_user)


def assert_coe_upload(upload_id, empty=False, metadata={}):
    coe_upload = Upload.from_upload_id(upload_id)

    if empty:
        assert coe_upload is None
    else:
        assert coe_upload is not None
        assert len(coe_upload.calcs) > 0
        for calc in coe_upload.calcs:
            assert_coe_calc(calc, metadata=metadata)

        if '_upload_time' in metadata:
            assert coe_upload.created.isoformat()[:26] == metadata['_upload_time']


def assert_coe_calc(calc: Calc, metadata={}):
    assert int(calc.pid) == int(metadata.get('_pid', calc.pid))
    assert calc.calc_id == metadata.get('_checksum', calc.calc_id)

    # calc data
    assert len(calc.filenames) == 5
    assert calc.chemical_formula is not None

    # user meta data
    assert calc.comment == metadata.get('comment', None)
    assert sorted(calc.references) == sorted(metadata.get('references', []))
    assert calc.uploader is not None
    assert calc.uploader.user_id == metadata.get('_uploader', calc.uploader.user_id)
    assert sorted(user.user_id for user in calc.coauthors) == sorted(metadata.get('coauthors', []))
    assert sorted(user.user_id for user in calc.shared_with) == sorted(metadata.get('shared_with', []))
    assert calc.with_embargo == metadata.get('with_embargo', False)


@pytest.mark.timeout(10)
def test_add_upload(clean_repository_db, processed_upload):
    empty = processed_upload.total_calcs == 0

    Upload.add(processed_upload)
    assert_coe_upload(processed_upload.upload_id, empty=empty)


@pytest.mark.timeout(10)
def test_add_upload_metadata(clean_repository_db, processed_upload, other_test_user, test_user):
    empty = processed_upload.total_calcs == 0

    metadata = {
        'comment': 'test comment',
        'with_embargo': True,
        'references': ['http://external.ref/one', 'http://external.ref/two'],
        '_uploader': other_test_user.user_id,
        'coauthors': [test_user.user_id],
        '_checksum': '1',
        '_upload_time': datetime.datetime.now().isoformat(),
        '_pid': 256
    }

    Upload.add(processed_upload, metadata=metadata)
    assert_coe_upload(processed_upload.upload_id, empty=empty, metadata=metadata)


class TestDataSets:

    @pytest.fixture(scope='function')
    def datasets(self, clean_repository_db):
        clean_repository_db.begin()
        one = Calc()
        two = Calc()
        three = Calc()
        clean_repository_db.add(one)
        clean_repository_db.add(two)
        clean_repository_db.add(three)
        one.children.append(two)
        two.children.append(three)
        clean_repository_db.commit()

        return one, two, three

    def assert_datasets(self, datasets, id_list):
        assert sorted([ds.id for ds in datasets]) == sorted(id_list)

    def test_all(self, datasets):
        one, two, three = datasets
        self.assert_datasets(one.all_datasets, [])
        self.assert_datasets(two.all_datasets, [one.coe_calc_id])
        self.assert_datasets(three.all_datasets, [one.coe_calc_id, two.coe_calc_id])

    def test_direct(self, datasets):
        one, two, three = datasets
        self.assert_datasets(one.direct_datasets, [])
        self.assert_datasets(two.direct_datasets, [one.coe_calc_id])
        self.assert_datasets(three.direct_datasets, [two.coe_calc_id])
