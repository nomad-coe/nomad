import requests
import json
import pytest
from copy import deepcopy

from tests.logtransfer.test_logtransfer import standard_test_log_record
from nomad.config import config
from nomad.logtransfer import gzip_bytes


@pytest.mark.parametrize(
    'is_gzip',
    [pytest.param(True, id='with gzip'), pytest.param(False, id='without gzip')],
)
@pytest.mark.timeout(3)
def test_single_logrecord(is_gzip, api_v1, central_logstash_mock):
    record = standard_test_log_record() + b'\n'

    if is_gzip:
        send_record_post = gzip_bytes(deepcopy(record))
        headers = {'Content-Encoding': 'gzip'}
    else:
        send_record_post = deepcopy(record)
        headers = None

    target = f'{config.client.url}/v1/federation/logs/'

    with central_logstash_mock as server:
        ret = requests.post(target, data=send_record_post, headers=headers)
        # print(f'send request {send_record_post}')
        server.handle_request()

    assert ret.status_code == 200
    assert int(json.loads(ret.content)['filesize']) > 0
    assert len(central_logstash_mock.received_content) == 1

    log_received_json = json.loads(deepcopy(central_logstash_mock.received_content[0]))
    # print(f'log_received_json={log_received_json}')
    log_sent_json = json.loads(record)

    assert 'ip_address' not in log_sent_json
    assert 'ip_address' in log_received_json.keys()
    assert log_received_json['ip_address'] == 'testclient'

    log_received_json.pop('ip_address')
    assert json.dumps(log_sent_json, sort_keys=True) == json.dumps(
        log_received_json, sort_keys=True
    )


@pytest.mark.parametrize(
    'send_header_ip,expected_ip',
    [
        pytest.param('127.0.0.1', '127.0.0.1', id='single ip'),
        pytest.param('127.0.0.1, 127.99.99.99', '127.0.0.1', id='concatenate ip'),
    ],
)
@pytest.mark.timeout(3)
def test_nginx_ipaddress_header(
    send_header_ip, expected_ip, api_v1, central_logstash_mock
):
    send_record = standard_test_log_record() + b'\n'
    target = f'{config.client.url}/v1/federation/logs/'

    with central_logstash_mock as server:
        header = {'X-Forwarded-For': send_header_ip}
        ret = requests.post(target, data=send_record, headers=header)
        server.handle_request()

    assert ret.status_code == 200
    assert len(central_logstash_mock.received_content) == 1
    actual_ip = json.loads(central_logstash_mock.received_content[0]).pop('ip_address')
    assert actual_ip == expected_ip


@pytest.mark.timeout(3)
def test_two_records_with_redundant_newline(api_v1, central_logstash_mock):
    target = f'{config.client.url}/v1/federation/logs/'

    log1 = standard_test_log_record(msg='testmsg1') + b'\n'
    log2 = standard_test_log_record(msg='testmsg2') + b'\n'

    # insert an unnecessary newline between the messages
    # this should simply be ignored and not error
    log_both = log1 + b'\n' + log2

    with central_logstash_mock as server:
        ret = requests.post(target, data=log_both)
        server.handle_request()

    assert ret.status_code == 200
    assert len(central_logstash_mock.received_content) == 2

    msg1_recv_json = json.loads(central_logstash_mock.received_content[0])
    msg2_recv_json = json.loads(central_logstash_mock.received_content[1])

    assert msg1_recv_json.pop('ip_address') == 'testclient'
    assert msg2_recv_json.pop('ip_address') == 'testclient'

    assert json.dumps(msg1_recv_json, sort_keys=False) == log1.decode().rstrip()
    assert json.dumps(msg2_recv_json, sort_keys=False) == log2.decode().rstrip()


@pytest.mark.parametrize(
    'msg,size',
    [
        pytest.param('testmsg', 7, id='no json format'),
        pytest.param('\n \n \t', 5, id='only whitespace'),
    ],
)
@pytest.mark.timeout(3)
def test_invalid_logrecord(msg, size, api_v1, central_logstash_mock):
    target = f'{config.client.url}/v1/federation/logs/'
    invalid = msg.encode()

    with central_logstash_mock as server:
        ret = requests.post(target, data=invalid)
        server.handle_request()

    assert ret.status_code == 200
    assert len(central_logstash_mock.received_content) == 0
    assert int(json.loads(ret.content)['filesize']) == size
