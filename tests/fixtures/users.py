import pytest

from nomad import infrastructure
from nomad.config import config
from nomad.datamodel import User
from tests.utils import test_user_uuid

admin_user_id = test_user_uuid(0)

test_users = {
    test_user_uuid(0): dict(username='admin', email='admin', user_id=test_user_uuid(0)),
    test_user_uuid(1): dict(
        username='scooper',
        email='sheldon.cooper@nomad-coe.eu',
        first_name='Sheldon',
        last_name='Cooper',
        user_id=test_user_uuid(1),
        is_oasis_admin=True,
    ),
    test_user_uuid(2): dict(
        username='lhofstadter',
        email='leonard.hofstadter@nomad-fairdi.tests.de',
        first_name='Leonard',
        last_name='Hofstadter',
        user_id=test_user_uuid(2),
    ),
}


@pytest.fixture(scope='session')
def test_user_molds():
    label_num = {
        'admin_user': 0,
        'test_user': 1,
        'other_test_user': 2,
    }
    return {label: test_users[test_user_uuid(num)] for label, num in label_num.items()}


@pytest.fixture(scope='session', autouse=True)
def configure_admin_user_id(monkeysession):
    monkeysession.setattr('nomad.config.services.admin_user_id', admin_user_id)


class KeycloakMock:
    def __init__(self):
        self.id_counter = 2
        self.users = dict(**test_users)

    def tokenauth(self, access_token: str):
        if access_token in self.users:
            return User(**self.users[access_token])
        else:
            raise infrastructure.KeycloakError('user does not exist')

    def add_user(self, user, *args, **kwargs):
        self.id_counter += 1
        user.user_id = test_user_uuid(self.id_counter)
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
            User(**test_user)
            for test_user in self.users.values()
            if query in ' '.join([str(value) for value in test_user.values()])
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
def proc_infra(worker, elastic_function, mongo_function, raw_files_function):
    """Combines all fixtures necessary for processing (elastic, worker, files, mongo)"""
    return dict(elastic=elastic_function)


@pytest.fixture(scope='function')
def with_oasis_user_management(monkeypatch):
    from nomad.infrastructure import OasisUserManagement

    monkeypatch.setattr('nomad.infrastructure.user_management', OasisUserManagement())
    yield
    monkeypatch.setattr('nomad.infrastructure.user_management', _user_management)


@pytest.fixture(scope='module')
def test_user():
    return User(**test_users[test_user_uuid(1)])


@pytest.fixture(scope='module')
def other_test_user():
    return User(**test_users[test_user_uuid(2)])


@pytest.fixture(scope='module')
def admin_user():
    return User(**test_users[test_user_uuid(0)])


@pytest.fixture(scope='module')
def test_users_dict(test_user, other_test_user, admin_user):
    return {
        'test_user': test_user,
        'other_test_user': other_test_user,
        'admin_user': admin_user,
    }
