import pytest
from fastapi.testclient import TestClient

from nomad.app_fastapi.main import app


@pytest.fixture(scope='session')
def client():
    return TestClient(app, base_url='http://testserver/api/v1/')
