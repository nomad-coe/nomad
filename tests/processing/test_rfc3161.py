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

import datetime

import httpx
import pytest
import rfc3161ng

from nomad.datamodel.datamodel import RFC3161Timestamp
from nomad.processing.data import get_rfc3161_token


@pytest.mark.parametrize('server,cert,result', [
    pytest.param('http://zeitstempel.dfn.de', None, True, id='zeitstempel.dfn.de'),
    pytest.param('http://timestamp.sectigo.com', None, True, id='timestamp.sectigo.com'),
    pytest.param('https://freetsa.org/tsr', 'https://freetsa.org/files/tsa.crt', True, id='freetsa.org/tsr'),
    pytest.param(
        'http://timestamp.digicert.com/',
        'https://knowledge.digicert.com/content/dam/digicertknowledgebase/attachments/time-stamp/TSACertificate.cer',
        True,
        id='timestamp.digicert.com-correct-cert'),
    pytest.param(
        'http://timestamp.digicert.com/',
        'https://freetsa.org/files/tsa.crt',
        False,
        id='timestamp.digicert.com-wrong-cert'),
])
def test_rfc3161ng_timestamp(server, cert, result, monkeysession):
    # this is due to requests being used by rfc3161ng
    # requests methods are modified in conftest.py which prohibits calling external servers
    monkeysession.setattr('requests.get', httpx.get)
    monkeysession.setattr('requests.post', httpx.post)

    token = get_rfc3161_token('test_hash', server=server, cert=cert)
    assert (token is not None) is result
    if result:
        rfc3161ng_time = rfc3161ng.get_timestamp(token)
        assert rfc3161ng_time < datetime.timedelta(seconds=5) + datetime.datetime.now()
        metadata = RFC3161Timestamp()
        metadata.token = token
        new_metadata = RFC3161Timestamp.m_from_dict(metadata.m_to_dict())
        assert new_metadata.token == token
        assert rfc3161ng.get_timestamp(new_metadata.token) == rfc3161ng_time
