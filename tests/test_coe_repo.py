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
import json
import datetime

from nomad.coe_repo import User, Calc, CalcMetaData, StructRatio, Upload, add_upload, \
    UserMetaData, Citation, MetaDataCitation, Shareship, CoAuthorship, Ownership

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


def assert_coe_upload(upload_hash, repository_db, empty=False, meta_data={}):
    coe_uploads = repository_db.query(Upload).filter_by(upload_name=upload_hash)
    if empty:
        assert coe_uploads.count() == 0
    else:
        assert coe_uploads.count() == 1
        coe_upload = coe_uploads.first()
        coe_upload_id = coe_upload.upload_id
        one_calc_exist = False
        for calc in repository_db.query(Calc).filter_by(origin_id=coe_upload_id):
            one_calc_exist = True
            assert calc.origin_id == coe_upload_id
            assert_coe_calc(calc, repository_db, meta_data=meta_data)

        if '_upload_time' in meta_data:
            assert coe_upload.created.isoformat()[:26] == meta_data['_upload_time']

        assert one_calc_exist


def assert_coe_calc(calc, repository_db, meta_data={}):
    calc_id = calc.calc_id
    calc_meta_data = repository_db.query(CalcMetaData).filter_by(calc_id=calc_id).first()

    assert calc_meta_data is not None
    assert calc_meta_data.calc is not None
    assert calc_meta_data.chemical_formula is not None
    filenames = calc_meta_data.filenames.decode('utf-8')
    assert len(json.loads(filenames)) == 5

    # struct ratio
    struct_ratio = repository_db.query(StructRatio).filter_by(calc_id=calc_id).first()
    assert struct_ratio is not None
    assert struct_ratio.chemical_formula == calc_meta_data.chemical_formula
    assert struct_ratio.formula_units == 1

    # pid
    if '_pid' in meta_data:
        assert calc_id == int(meta_data['_pid'])

    # checksum
    if '_checksum' in meta_data:
        calc.checksum == meta_data['_checksum']

    # comments
    comment = repository_db.query(UserMetaData).filter_by(
        label=meta_data.get('comment', 'not existing comment'),
        calc_id=calc_id).first()
    if 'comment' in meta_data:
        assert comment is not None
    else:
        assert comment is None

    # references
    if 'references' in meta_data:
        for reference in meta_data['references']:
            citation = repository_db.query(Citation).filter_by(
                value=reference, kind='EXTERNAL').first()
            assert citation is not None
            assert repository_db.query(MetaDataCitation).filter_by(
                citation_id=citation.citation_id, calc_id=calc_id).first() is not None
    else:
        repository_db.query(MetaDataCitation).filter_by(calc_id=calc_id).first() is None

    # coauthors
    if 'coauthors' in meta_data:
        for coauthor in meta_data['coauthors']:
            assert repository_db.query(CoAuthorship).filter_by(
                user_id=coauthor, calc_id=calc_id).first() is not None
    else:
        assert repository_db.query(CoAuthorship).filter_by(calc_id=calc_id).first() is None

    # coauthors
    if 'shared_with' in meta_data:
        for coauthor in meta_data['shared_with']:
            assert repository_db.query(Shareship).filter_by(
                user_id=coauthor, calc_id=calc_id).first() is not None
    else:
        assert repository_db.query(Shareship).filter_by(calc_id=calc_id).first() is None

    # ownership
    owners = repository_db.query(Ownership).filter_by(calc_id=calc_id)
    assert owners.count() == 1
    if '_uploader' in meta_data:
        assert owners.first().user_id == meta_data['_uploader']

    # embargo/restriction/permission
    user_meta_data = repository_db.query(UserMetaData).filter_by(
        calc_id=calc_meta_data.calc_id).first()
    assert user_meta_data is not None
    assert user_meta_data.permission == (1 if meta_data.get('with_embargo', False) else 0)


@pytest.mark.timeout(10)
def test_add_upload(clean_repository_db, processed_upload):
    empty = processed_upload.total_calcs == 0

    processed_upload.upload_hash = str(1)
    add_upload(processed_upload)
    assert_coe_upload(processed_upload.upload_hash, clean_repository_db, empty=empty)

    processed_upload.upload_hash = str(2)
    add_upload(processed_upload)
    assert_coe_upload(processed_upload.upload_hash, clean_repository_db, empty=empty)


@pytest.mark.timeout(10)
def test_add_upload_metadata(clean_repository_db, processed_upload, other_test_user, test_user):
    empty = processed_upload.total_calcs == 0

    meta_data = {
        'comment': 'test comment',
        'with_embargo': True,
        'references': ['http://external.ref/one', 'http://external.ref/two'],
        '_uploader': other_test_user.user_id,
        'coauthors': [test_user.user_id],
        '_checksum': 1,
        '_upload_time': datetime.datetime.now().isoformat(),
        '_pid': 256
    }

    add_upload(processed_upload, meta_data=meta_data)
    assert_coe_upload(processed_upload.upload_hash, clean_repository_db, empty=empty, meta_data=meta_data)
