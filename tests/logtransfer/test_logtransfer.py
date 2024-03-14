import os.path
import time
import json
import pytest

from nomad import logtransfer
from logging import LogRecord
import logging
from nomad.utils.structlogging import LogstashFormatter
from nomad import config

logtransfer.enable_log_print = False


def standard_test_log_record(msg='testmsg') -> bytes:
    record = LogRecord(
        name='test',
        msg=msg,
        args=(),
        exc_info=None,
        lineno=1,
        level=logging.INFO,
        pathname='testpath',
    )
    return LogstashFormatter().format(record)


def test_logstash_format():
    logstash_log = standard_test_log_record('testmsg')
    logstash_formatted = json.loads(logstash_log)  # asserts valid json format
    assert logstash_formatted['event'] == 'testmsg'


def test_get_all_rotated_files():
    logtransfer.clear_logfiles()
    logger, handler = logtransfer._initialize_logger_and_handler()

    n_rollover = 20  # more than 10 to check also the correct sorting
    handler.backupCount = n_rollover

    expected = [
        os.path.abspath(
            os.path.join(config.fs.tmp, config.logtransfer.log_filename + f'.{i}')
        )
        for i in range(n_rollover, 0, -1)
    ]

    for i in range(n_rollover):
        logger.info(f'msg{i}')
        handler.doRollover()

    actual = logtransfer.get_all_rotated_logfiles()

    for i, a in enumerate(actual):
        assert a == expected[i]


def test_rotating_logfiles():
    logtransfer.clear_logfiles()

    # write some message to logger
    logger, handler = logtransfer._initialize_logger_and_handler()
    logger.info('msg')
    handler.doRollover()

    backup_files = logtransfer.get_all_rotated_logfiles()
    assert len(backup_files) == 1
    assert (
        os.path.basename(backup_files[0])
        == os.path.basename(logtransfer.get_log_filepath()) + '.1'
    )

    logtransfer.clear_logfiles()


@pytest.mark.timeout(3)
def test_logstash_submit_single_log(logtransfer_rollover_time):
    # the newline is important, because the server reads only lines
    logstash_line = standard_test_log_record() + b'\n'

    ret = logtransfer_rollover_time.send(logstash_line)

    assert ret > 0

    while logtransfer.is_empty_logfile():
        # wait until server has processed and written the logs to the file
        time.sleep(0.001)

    with open(logtransfer.get_log_filepath(), 'rb') as f:
        file_content = f.read()
        assert file_content.rstrip() == logstash_line.rstrip()


@pytest.mark.timeout(3)
def test_logstash_submit_multiple_logs(
    logtransfer_rollover_time, central_logstash_mock
):
    n_messages = 5
    log_event = [f'testlog{i}' for i in range(5)]

    for i in range(n_messages):
        logstash_line = standard_test_log_record(msg=log_event[i]) + b'\n'
        logtransfer_rollover_time.send(logstash_line)

    is_successful = False
    while not is_successful:
        # The logfile may also contain a subset of the messages that were sent
        # at this point leading to an AssertionError. This loop frequently
        # checks if at some all asserts are True and then terminates
        if not logtransfer.is_empty_logfile():
            with open(logtransfer.get_log_filepath(), 'rb') as f:
                file_content = f.read()

            try:
                for me in log_event:
                    assert f'"event": "{me}"'.encode() in file_content
            except AssertionError:
                break
            else:
                # is executed if no AssertionError error was raised
                is_successful = True

        if not is_successful:
            time.sleep(0.001)  # slow down checks


@pytest.mark.timeout(3)
def test_logstash_rollover_time(logtransfer_rollover_time):
    testlog = ['testmsg0', 'testmsg1']

    log1 = standard_test_log_record(msg=testlog[0]) + b'\n'
    ret1 = logtransfer_rollover_time.send(log1)
    assert ret1 > 0

    # this sleep is needed to trigger the rollover after some time interval
    # (check the fixture for the actual setting)
    time.sleep(0.15)

    log2 = standard_test_log_record(msg=testlog[1]) + b'\n'
    ret2 = logtransfer_rollover_time.send(log2)
    assert ret2 > 0

    is_successful_assert = False
    while not is_successful_assert:
        # finish the test as soon as all asserts are passed (this can vary
        # depending on when the server has processed the log records)
        try:
            rolled_log_filepath = logtransfer.get_log_filepath() + '.1'
            assert not logtransfer.is_empty_logfile()
            assert os.path.exists(rolled_log_filepath)

            with open(logtransfer.get_log_filepath(), 'rb') as f:
                # check active logfile
                file_content = f.read()
                assert testlog[0].encode() not in file_content
                assert testlog[1].encode() in file_content

            with open(rolled_log_filepath, 'rb') as f:
                # check backup file
                file_content = f.read()
                assert testlog[0].encode() in file_content
                assert testlog[1].encode() not in file_content
        except AssertionError:
            time.sleep(0.001)
        else:
            # executed if no AssertionError was raised
            is_successful_assert = True


def assert_equal_content_logfiles(expected_logs):
    # run until assertions pass
    while True:
        try:
            with open(logtransfer.get_log_filepath(), 'rb') as f:
                active_logfile = f.read()

            rotated_files = logtransfer.get_all_rotated_logfiles()

            # there should be at least one rollover logfile (due to exceeding filesize)
            assert len(rotated_files) > 0

            backup_content = []
            for filepath in rotated_files:
                with open(filepath, 'rb') as f:
                    backup_content.append(f.read())

            actual_logs = b''.join(backup_content) + active_logfile

            assert actual_logs == expected_logs
        except (FileNotFoundError, AssertionError):
            time.sleep(0.2)
        else:
            break  # break loop when no error was raised


@pytest.mark.timeout(10)
def test_logstash_rollover_space(logtransfer_rollover_space):
    n_bytes, idx = 0, 0
    expected_logs = []
    while True:
        log = standard_test_log_record(msg=f'testmsg{idx}') + b'\n'
        logtransfer_rollover_space.send(log)
        expected_logs.append(log)

        n_bytes += len(log)
        idx += 1

        # loop until there are three rolled-over logfiles
        if len(logtransfer.get_all_rotated_logfiles()) > 2:
            break
        else:
            # slow down logging a bit, so that server can perform rollovers
            time.sleep(0.01)

    expected_logs = b''.join(expected_logs)
    assert_equal_content_logfiles(expected_logs)


@pytest.mark.timeout(3)
def test_logstash_rollover_space_large_message(logtransfer_rollover_space):
    expected_logs = []

    n_bytes, idx = 0, 0
    while n_bytes < 10000:
        log = standard_test_log_record(msg=f'testmsg{idx}')
        expected_logs.append(log)
        n_bytes += len(log)

    # send one large messages
    expected_logs = b'\n'.join(expected_logs) + b'\n'
    ret = logtransfer_rollover_space.send(expected_logs)

    assert ret > 0
    assert_equal_content_logfiles(expected_logs)


@pytest.mark.timeout(3)
def test_logtransfer_to_federation_backend(
    api_v1, logtransfer_rollover_time, central_logstash_mock
):
    logstash_line = standard_test_log_record() + b'\n'

    with central_logstash_mock as server:
        # server is closed after context manager exits
        ret = logtransfer_rollover_time.send(logstash_line)
        server.handle_request()

    # only one message was sent
    assert len(central_logstash_mock.received_content) == 1

    # transform to json (also validates that both are valid)
    sent_json = json.loads(logstash_line.rstrip())
    received_json = json.loads(central_logstash_mock.received_content[0])

    # validate that there is an augmented key 'ip_address' on central
    assert 'ip_address' not in sent_json.keys()
    assert 'ip_address' in received_json.keys()

    # remove key and convert back to string to compare that both logs have
    # identical key/values
    received_json.pop('ip_address')
    assert json.dumps(sent_json, sort_keys=True) == json.dumps(
        received_json, sort_keys=True
    )
