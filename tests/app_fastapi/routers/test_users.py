
def test_me(client, test_user_auth):
    response = client.get('users/me', headers=test_user_auth)
    assert response.status_code == 200


def test_me_auth_required(client):
    response = client.get('users/me')
    assert response.status_code == 401


def test_me_auth_bad_token(client):
    response = client.get('users/me', headers={'Authentication': 'Bearer NOTATOKEN'})
    assert response.status_code == 401
