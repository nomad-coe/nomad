# Copyright 2019  Alvin Noe Ladines, Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Interface to archive api. A dict of query parameters and a query schema similar to
the archive json format can be provided to filter archive.

q = ArchiveQuery({'atoms':'Si'})
metainfo = q.query()

for c in metainfo.calcs:
    print(c.section_run.section_single_configuration_calculation[0]({'energy_total':None}))
"""

import requests
import json

from nomad.app.api.common import query_api_url
from nomad.archive_library.metainfo import ArchiveMetainfo


class ArchiveQuery:
    def __init__(self, q_params, q_schema=None):
        self._q_params = q_params
        self._q_schema = q_schema
        self._archive_path = 'archive'
        self._query_path = 'query'
        self._q_params['scroll'] = q_params.get('scroll', True)
        self._q_params['per_page'] = q_params.get('per_page', 10)
        self._scroll_id = None
        self._data = []

    def _api_query(self):
        if self._scroll_id is not None:
            self._q_params['scroll_id'] = self._scroll_id
        url = query_api_url(
            self._archive_path, self._query_path, query_string=self._q_params)
        data = {'results': self._q_schema}

        response = requests.post(
            url, content_type='application/json', data=json.dumps(data))
        if response.status_code != 200:
            raise Exception('Query returned %s' % response.status_code)

        data = response.json
        if not isinstance(data, dict):
            data = data()
        results = data.get('results', None)
        scroll = data.get('scroll', None)
        if scroll:
            self._scroll_id = data.get('scroll_id', None)

        return results

    def _get_data(self):
        results = self._api_query()
        while results:
            self._data += results
            results = self._api_query()
            if self._scroll_id is None:
                break

    def query(self):
        self._get_data()
        if self._data:
            metainfo = ArchiveMetainfo(archive_data=self._data)
            return metainfo
