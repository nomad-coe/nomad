from nomad.user import User


def assert_user(user, reference):
    assert user is not None
    assert user.user_id == reference.user_id


def test_token_authorize(test_user):
    user = User.verify_auth_token(test_user.email)
    assert_user(user, test_user)


def test_password_authorize(test_user):
    user = User.verify_user_password(test_user.email, 'password')
    assert_user(user, test_user)
