import pytest
from urllib.parse import urlencode


def perform_get_token_test(client, http_method, status_code, username, password):
    if http_method == 'post':
        response = client.post(
            'auth/token',
            data=dict(username=username, password=password))
    else:
        response = client.get('auth/token?%s' % urlencode(
            dict(username=username, password=password)))

    assert response.status_code == status_code


@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_get_token(client, test_user, http_method):
    perform_get_token_test(client, http_method, 200, test_user.username, 'password')


@pytest.mark.parametrize('http_method', ['post', 'get'])
def test_get_token_bad_credentials(client, http_method):
    perform_get_token_test(client, http_method, 401, 'bad', 'credentials')
