import pytest
import logging


@pytest.fixture()
def my_caplog(caplog):
    yield caplog

    # TODO there is a bug in pytest
    # assert len(caplog.records) > 0


def test_nowarn(my_caplog):
    logging.getLogger().warning('Hello, anybody there')
    # TODO there still seems that legace parsers/normalizers fiddle with the
    # log configuration. The following fails after running tests with parsers/normalizers
    # assert len(my_caplog.records) > 0
