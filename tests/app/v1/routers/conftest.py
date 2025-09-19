import pytest


@pytest.fixture(scope='function')
def fastapi_cache(mongo_function, monkeypatch):
    from fastapi_cache import FastAPICache

    from nomad.app.main import MongoBackend

    FastAPICache.init(backend=MongoBackend())
    yield
    FastAPICache.reset()
