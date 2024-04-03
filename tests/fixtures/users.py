"""
User fixtures:
- user0: admin user
- user1: default user to use
- user2, user3: additional users for access or interaction tests
"""

import pytest

from nomad import infrastructure
from nomad.config import config
from nomad.datamodel import User
from tests.utils import fake_user_uuid

admin_user_id = fake_user_uuid(0)

users = {
    fake_user_uuid(0): dict(username='admin', email='admin', user_id=fake_user_uuid(0)),
    fake_user_uuid(1): dict(
        username='scooper',
        email='sheldon.cooper@nomad-coe.eu',
        first_name='Sheldon',
        last_name='Cooper',
        user_id=fake_user_uuid(1),
        is_oasis_admin=True,
    ),
    fake_user_uuid(2): dict(
        username='lhofstadter',
        email='leonard.hofstadter@nomad-fairdi.tests.de',
        first_name='Leonard',
        last_name='Hofstadter',
        user_id=fake_user_uuid(2),
    ),
    fake_user_uuid(3): dict(
        username='hwolowitz',
        email='howard.wolowitz@nomad-fairdi.tests.de',
        first_name='Howard',
        last_name='Wolowitz',
        user_id=fake_user_uuid(3),
    ),
}


@pytest.fixture(scope='session')
def user_molds():
    label_num = {f'user{i}': i for i in range(4)}
    return {label: users[fake_user_uuid(num)] for label, num in label_num.items()}


@pytest.fixture(scope='session')
def user0():
    return User(**users[fake_user_uuid(0)])


@pytest.fixture(scope='session')
def user1():
    return User(**users[fake_user_uuid(1)])


@pytest.fixture(scope='session')
def user2():
    return User(**users[fake_user_uuid(2)])


@pytest.fixture(scope='session')
def user3():
    return User(**users[fake_user_uuid(3)])


@pytest.fixture(scope='session')
def users_dict(user0, user1, user2, user3):
    return {
        'user0': user0,
        'user1': user1,
        'user2': user2,
        'user3': user3,
    }


@pytest.fixture(scope='session', autouse=True)
def configure_admin_user_id(monkeysession):
    monkeysession.setattr('nomad.config.services.admin_user_id', admin_user_id)


class KeycloakMock:
    def __init__(self):
        self.id_counter = 3
        self.users = dict(**users)

    def tokenauth(self, access_token: str):
        if access_token in self.users:
            return User(**self.users[access_token])
        else:
            raise infrastructure.KeycloakError('user does not exist')

    def add_user(self, user, *args, **kwargs):
        self.id_counter += 1
        user.user_id = fake_user_uuid(self.id_counter)
        user.username = (user.first_name[0] + user.last_name).lower()
        self.users[user.user_id] = dict(
            email=user.email,
            username=user.username,
            first_name=user.first_name,
            last_name=user.last_name,
            user_id=user.user_id,
        )

    def get_user(self, user_id=None, username=None, email=None):
        if user_id is not None:
            return User(**self.users[user_id])
        elif username is not None:
            for user_id, user_values in self.users.items():
                if user_values['username'] == username:
                    return User(**user_values)
            raise KeyError('Only test user usernames are recognized')
        elif email is not None:
            for user_id, user_values in self.users.items():
                if user_values['email'] == email:
                    return User(**user_values)
            raise KeyError('Only test user emails are recognized')
        else:
            assert False, 'no token based get_user during tests'

    def search_user(self, query):
        return [
            User(**user)
            for user in self.users.values()
            if query in ' '.join([str(value) for value in user.values()])
        ]

    def basicauth(self, username: str, password: str) -> str:
        for user in self.users.values():
            if user['username'] == username or user['email'] == username:
                return user['user_id']

        raise infrastructure.KeycloakError()


config.keycloak.realm_name = 'fairdi_nomad_test'
config.keycloak.password = 'password'

_keycloak = infrastructure.keycloak
_user_management = infrastructure.user_management


# use a session fixture in addition to the function fixture, to ensure mocked keycloak
# before other class, module, etc. scoped function are run
@pytest.fixture(scope='session', autouse=True)
def mocked_keycloak_session(monkeysession):
    monkeysession.setattr('nomad.infrastructure.keycloak', KeycloakMock())
    monkeysession.setattr('nomad.infrastructure.user_management', KeycloakMock())


@pytest.fixture(scope='function', autouse=True)
def mocked_keycloak(monkeypatch):
    monkeypatch.setattr('nomad.infrastructure.keycloak', KeycloakMock())
    monkeypatch.setattr('nomad.infrastructure.user_management', KeycloakMock())


@pytest.fixture(scope='function')
def keycloak(monkeypatch):
    monkeypatch.setattr('nomad.infrastructure.keycloak', _keycloak)
    monkeypatch.setattr('nomad.infrastructure.user_management', _user_management)


@pytest.fixture(scope='function')
def with_oasis_user_management(monkeypatch):
    from nomad.infrastructure import OasisUserManagement

    monkeypatch.setattr('nomad.infrastructure.user_management', OasisUserManagement())
    yield
    monkeypatch.setattr('nomad.infrastructure.user_management', _user_management)
