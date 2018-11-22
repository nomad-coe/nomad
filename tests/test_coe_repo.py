import pytest

from nomad.coe_repo import User, Calc, CalcMetaData, Upload, add_upload

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


def test_rollback(repository_db):
    repository_db.begin()

    calc = Calc(checksum='test')
    repository_db.add(calc)
    repository_db.flush()
    calc_id = calc.calc_id

    repository_db.rollback()

    assert repository_db.query(Calc).filter_by(calc_id=calc_id).first() is None


def assert_coe_upload(upload_hash, repository_db, empty=False):
    coe_upload = repository_db.query(Upload).filter_by(upload_name=upload_hash).first()
    if empty:
        assert coe_upload is None
    else:
        assert coe_upload is not None
        coe_upload_id = coe_upload.upload_id
        for calc in repository_db.query(Calc).filter_by(origin_id=coe_upload_id):
            assert calc.origin_id == coe_upload_id
            metadata = repository_db.query(CalcMetaData).filter_by(calc_id=calc.calc_id).first()
            assert metadata is not None
            assert metadata.chemical_formula is not None


@pytest.mark.timeout(10)
def test_add_upload(repository_db, processed_upload):
    coe_upload_id = add_upload(processed_upload, restricted=False)
    if coe_upload_id:
        assert_coe_upload(processed_upload.upload_hash, repository_db)

    coe_upload_id = add_upload(processed_upload, restricted=False)
    if coe_upload_id:
        assert_coe_upload(processed_upload.upload_hash, repository_db)
