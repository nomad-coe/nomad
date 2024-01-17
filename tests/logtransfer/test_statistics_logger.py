import json
import time

import pytest
import requests

import nomad.statistics_logger as stats


@pytest.mark.timeout(3)
def test_statistics_service(
    logstash_enabled, central_logstash_mock, collect_statistics, api_v1, elastic_infra
):
    # Note: the "elastic" fixture is required here so that info/ request at the
    # FastAPI backend actually collects statistics.

    statistics_process = stats.start_statistics_logger(fork_mode='process')
    # statistics_thread = stats.start_statistics_logger(fork_mode='thread')

    with central_logstash_mock as server:
        # timeout must be smaller than the number set in fixture collect_statistics
        timeout = 0.5
        server.set_request_timeout(timeout)

        # server closes when context manager is left
        server.handle_request()

    while statistics_process.is_alive():
        statistics_process.kill()
        time.sleep(0.01)

    # assert that all current statistics end up in logtransfer file
    collected_statistics = json.loads(requests.get('info').content)['statistics']

    received_content = []
    for item in central_logstash_mock.received_content:
        received_content += [json.loads(item.decode())]

    for expected_key, expected_value in collected_statistics.items():
        for log_entry in received_content:
            try:
                actual_key = log_entry['nomad.statistics_logger.key']
                actual_value = log_entry['nomad.statistics_logger.value']
                if actual_key == expected_key and actual_value == expected_value:
                    # success: found key-value - break inner loop
                    break

            except KeyError:
                # there may be other logs that do not contain the statistics keys
                pass
        else:
            # executes if inner-loop was *not* terminated with 'break'
            assert False, f'Could not find {expected_key=}, {expected_value=} pair.'
