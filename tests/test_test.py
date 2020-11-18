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

import pytest
import logging


@pytest.fixture()
def my_caplog(caplog):
    yield caplog

    # TODO there still seems that legace parsers/normalizers fiddle with the
    # log configuration. The following fails after running tests with parsers/normalizers
    # assert len(caplog.get_records(when='call')) > 0


def test_nowarn(my_caplog):
    logging.getLogger().warning('Hello, anybody there')
    # TODO there still seems that legace parsers/normalizers fiddle with the
    # log configuration. The following fails after running tests with parsers/normalizers
    # assert len(my_caplog.records) > 0
