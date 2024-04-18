import logging
import json
import pytest
import os.path

from nomad import config, utils
from nomad.utils import structlogging
from nomad.logtransfer import transfer_logs


@pytest.fixture(scope='function')
def log_handler(monkeypatch, raw_files_function):
    from structlog.stdlib import LoggerFactory

    test_root_logger = logging.getLogger()
    monkeypatch.setattr(
        'nomad.utils.structlogging.default_factory',
        LoggerFactory(test_root_logger),
    )
    structlogging.add_logtransfer_handler(test_root_logger)

    handler = structlogging.get_logtransfer_handler(test_root_logger)
    assert handler is not None
    assert isinstance(handler.formatter, structlogging.LogtransferFormatter)

    return handler


def test_logtransfer_handler(log_handler):
    test_logger = utils.get_logger('nomad.tests')
    test_logger.info('test event', data='test data')

    assert os.path.exists(log_handler.baseFilename)
    with open(log_handler.baseFilename, 'r') as f:
        logs = f.readlines()

    assert len(logs) == 1
    record = json.loads(logs[0])
    assert record['event'] == 'test event'
    assert record['nomad.tests.data'] == 'test data'
    assert all(
        [
            key in record
            for key in [
                '@timestamp',
                'level',
                'event',
                'nomad.version',
                'nomad.commit',
                'nomad.service',
            ]
        ]
    )


@pytest.mark.parametrize('size', [1, 2])
@pytest.mark.timeout(3)
def test_transfer_logs(log_handler, monkeypatch, api_v1, central_logstash_mock, size):
    monkeypatch.setattr(
        'nomad.config.config.oasis.central_nomad_deployment_url', config.client.url
    )
    monkeypatch.setattr('nomad.config.config.logtransfer.transfer_threshold', 0)
    monkeypatch.setattr('nomad.config.config.logtransfer.enabled', True)
    monkeypatch.setattr('nomad.config.config.logtransfer.file_rollover_wait_time', 0)

    test_logger = utils.get_logger('nomad.tests')
    for _ in range(size):
        test_logger.info('test event', data='test data')

    assert os.path.exists(log_handler.baseFilename)
    with open(log_handler.baseFilename, 'r') as f:
        logs = f.readlines()
    assert len(logs) == size

    monkeypatch.setattr(
        'nomad.config.config.logtransfer.transfer_capacity', len(logs[0]) + 1
    )

    with central_logstash_mock as server:
        transfer_logs()
        server.handle_request()

    assert len(central_logstash_mock.received_content) == 1

    received_json = json.loads(central_logstash_mock.received_content[0])

    assert 'ip_address' in received_json.keys()

    received_json.pop('ip_address')
    assert json.dumps(json.loads(logs[0]), sort_keys=True) == json.dumps(
        received_json, sort_keys=True
    )
