#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import logging

import requests
import threading
import time
from multiprocessing import Process

from nomad.config import config
from nomad import utils


# The logger (typically with LogstashHandler enabled in the root logger) to
# submit the logs.
statistics_logstash_log = utils.get_logger(__name__)


class StatisticsLogger:
    """
    Service that frequently collects statistics and logs them in logstash format.

    Currently, the statistics are collected from the local FastAPI enpoint /info.
    The statistics logs are eventually transferred to
    * logstash directly (on central NOMAD)
    * logtransfer (on a NOMAD Oasis, see logtransfer.py).

    In the future it is also possible to collect statistics from different sources
    or perform aggregations.
    """

    def __init__(self):
        self.collect_interval = config.logtransfer.statistics_interval
        self.raise_exceptions = config.logtransfer.raise_unexpected_exceptions
        self.stop_thread_event = threading.Event()

    def _collect_statistics(self):
        """Collect statistics from the FastAPI /info endpoint."""
        info_response = requests.get(
            f'http://localhost:{config.services.api_port}/api/v1/info'
        )
        statistics = info_response.json()['statistics']
        return statistics

    def run_statistics_service(self):
        """Main loop of the statistics service."""

        while True:
            try:
                statistics = self._collect_statistics()

                for key, value in statistics.items():
                    statistics_logstash_log.info('statistics', key=key, value=value)

            except Exception as e:
                if self.raise_exceptions:
                    raise e

            # The event can be set to terminate the loop. This is mainly
            # required for testing.
            import sys

            sys.stdout.flush()
            time.sleep(self.collect_interval)


def start_statistics_logger(fork_mode='process'):
    """
    Initializes and starts a new process of the statistics logger.
    The returned object is the running process or thread.
    """

    statistics_logger = StatisticsLogger()

    if fork_mode == 'process':
        process = Process(target=statistics_logger.run_statistics_service)
        process.start()
        return process
    elif fork_mode == 'thread':
        statistics_thread = threading.Thread(
            target=statistics_logger.run_statistics_service
        )

        statistics_thread.setDaemon(True)
        statistics_thread.start()
        return statistics_thread
    else:
        raise ValueError(f'{fork_mode=} is not valid (select "process" or "thread").')
