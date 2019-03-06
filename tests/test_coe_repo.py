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
from passlib.hash import bcrypt

from nomad.coe_repo import User, Calc, Upload
from nomad.coe_repo.calc import PublishContext
from nomad import processing, parsing, datamodel


def assert_user(user, reference):
    assert user is not None
    assert user.user_id == reference.user_id
    assert user.email == reference.email


def test_token_authorize(test_user):
    user = User.verify_auth_token(test_user.first_name.lower())
    assert_user(user, test_user)


def test_password_authorize(test_user):
    user = User.verify_user_password(test_user.email, 'password')
    assert_user(user, test_user)


def assert_coe_upload(upload_id, upload: datamodel.UploadWithMetadata = None, user_metadata: dict = None):
    coe_upload = Upload.from_upload_id(upload_id)

    if upload is not None:
        calcs = list(upload.calcs)
    elif coe_upload is None:
        calcs = []
    else:
        calcs = list(calc.to_calc_with_metadata() for calc in coe_upload.calcs)

    if len(calcs) == 0:
        assert coe_upload is None
    else:
        assert coe_upload is not None
        assert len(coe_upload.calcs) == len(calcs)
        for coe_calc, calc in zip(coe_upload.calcs, calcs):
            if user_metadata is not None:
                calc.apply_user_metadata(user_metadata)

            assert_coe_calc(coe_calc, calc)

        if upload is not None and upload.upload_time is not None:
            assert coe_upload.created.isoformat()[:26] == upload.upload_time.isoformat()


def assert_coe_calc(coe_calc: Calc, calc: datamodel.CalcWithMetadata):
    if calc.pid is not None:
        assert coe_calc.pid == calc.pid

    # calc data
    assert len(coe_calc.files) == len(calc.files)
    assert coe_calc.formula == calc.formula

    # user meta data
    assert coe_calc.comment == calc.comment
    assert len(coe_calc.references) == len(calc.references)
    assert coe_calc.uploader is not None
    if calc.uploader is not None:
        assert coe_calc.uploader.user_id == calc.uploader.id
    assert sorted(user.user_id for user in coe_calc.coauthors) == sorted(user.id for user in calc.coauthors)
    assert sorted(user.user_id for user in coe_calc.shared_with) == sorted(user.id for user in calc.shared_with)
    if calc.with_embargo is not None:
        assert coe_calc.with_embargo == calc.with_embargo
    else:
        assert not coe_calc.with_embargo


def test_add_normalized_calc(postgres, normalized: parsing.LocalBackend, test_user):
    calc_with_metadata = normalized.to_calc_with_metadata()
    calc_with_metadata.uploader = test_user.to_popo()
    calc_with_metadata.files = [calc_with_metadata.mainfile, '1', '2', '3', '4']
    coe_calc = Calc()
    coe_calc.apply_calc_with_metadata(calc_with_metadata, PublishContext())

    assert_coe_calc(coe_calc, calc_with_metadata)


def test_add_normalized_calc_with_metadata(
        postgres, normalized: parsing.LocalBackend, example_user_metadata: dict):

    calc_with_metadata = normalized.to_calc_with_metadata()
    calc_with_metadata.files = [calc_with_metadata.mainfile, '1', '2', '3', '4']
    calc_with_metadata.apply_user_metadata(example_user_metadata)
    coe_calc = Calc(coe_calc_id=calc_with_metadata.pid)
    coe_calc.apply_calc_with_metadata(calc_with_metadata, PublishContext())

    assert_coe_calc(coe_calc, calc_with_metadata)


def test_add_upload(processed: processing.Upload):
    upload_with_metadata = processed.to_upload_with_metadata()
    Upload.publish(upload_with_metadata)(True)
    assert_coe_upload(processed.upload_id, upload_with_metadata)


def test_rollback_upload(processed: processing.Upload):
    upload_with_metadata = processed.to_upload_with_metadata()
    Upload.publish(upload_with_metadata)(False)
    assert Upload.from_upload_id(processed.upload_id) is None


# def test_large_upload(processed: processing.Upload, example_user_metadata):
#     processed.metadata = example_user_metadata
#     upload_with_metadata = processed.to_upload_with_metadata()
#     calcs = list(upload_with_metadata.calcs)

#     if len(calcs) == 0:
#         return

#     def many_calcs():
#         count = 0
#         while True:
#             for calc in calcs:
#                 calc.pid = count + 10
#                 yield calc
#                 count += 1
#                 if count > 1000:
#                     return

#     import time
#     start = time.time()
#     upload_with_metadata.calcs = many_calcs()
#     Upload.publish(upload_with_metadata)(True)
#     print('########### %d' % (time.time() - start))


def test_add_upload_with_metadata(processed, example_user_metadata):
    processed.metadata = example_user_metadata
    upload_with_metadata = processed.to_upload_with_metadata()
    Upload.publish(upload_with_metadata)(True)
    assert_coe_upload(
        processed.upload_id, upload_with_metadata)


@pytest.mark.parametrize('crypted', [True, False])
def test_create_user(postgres, crypted):
    password = bcrypt.encrypt('test_password', ident='2y') if crypted else 'test_password'
    data = dict(
        email='test@email.com', last_name='Teser', first_name='testi', password=password)

    user = User.create_user(**data, crypted=crypted)

    authenticated_user = User.verify_user_password('test@email.com', 'test_password')
    assert authenticated_user is not None
    assert user.user_id == authenticated_user.user_id
    assert user.get_auth_token() is not None


class TestDataSets:

    @pytest.fixture(scope='function')
    def datasets(self, postgres):
        postgres.begin()
        one = Calc()
        two = Calc()
        three = Calc()
        postgres.add(one)
        postgres.add(two)
        postgres.add(three)
        one.children.append(two)
        two.children.append(three)
        postgres.commit()

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
