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
import os.path

from nomad import config
from nomad.archive_library.metainfo import ArchiveMetainfo


class ArchiveQuery:
    def __init__(self, *args, **kwargs):
        self._archive_path = 'archive'
        self._query_path = 'query'
        self._archive_data = []
        self._scroll_id = None
        self._page = None
        self._query_params = {}
        if args:
            self._query_params = args[0]
        if kwargs:
            self._query_params.update(kwargs)
        self._archive_schema = self._query_params.pop('archive_data', None)
        self._authentication = self._query_params.pop('authentication', None)
        self._max_n_pages = self._query_params.pop('max_n_pages', 3)

    def _get_value(self, name, in_dict):
        if not isinstance(in_dict, dict):
            return
        for key, val in in_dict.items():
            if key == name:
                res = val
            else:
                res = self._get_value(name, val)
        return res

    def _set_value(self, name, value, in_dict):
        if not isinstance(in_dict, dict):
            return
        for key, val in in_dict.items():
            if key == name:
                in_dict[name] = value
                return
            else:
                self._set_value(name, value, val)
        in_dict[name] = value

    def _api_query(self):
        url = os.path.join(config.api_url(False), self._archive_path, self._query_path)
        data = self._query_params
        if not isinstance(self._archive_schema, list):
            data['results'] = [self._archive_schema]
        else:
            data['results'] = self._archive_schema

        if self._page is not None:
            # increment the page number
            self._set_value('page', self._page + 1, data)
        if self._scroll_id is not None:
            self._set_value('scroll_id', self._scroll_id, data)

        response = requests.post(url, headers=self._authentication, json=data)
        if response.status_code != 200:
            raise Exception('Query returned %s' % response.status_code)

        data = response.json
        if not isinstance(data, dict):
            data = data()

        results = data.get('results', None)
        scroll = data.get('scroll', None)
        if scroll:
            self._scroll_id = scroll.get('scroll_id', None)
        pagination = data.get('pagination', None)
        if pagination:
            self._page = pagination.get('page', None)

        return results

    def _get_archive_data(self):
        results = self._api_query()
        n_page = 0
        while results:
            self._archive_data += results
            results = self._api_query()
            n_page += 1
            if n_page >= self._max_n_pages:
                break

    def query(self):
        self._get_archive_data()
        if self._archive_data:
            metainfo = ArchiveMetainfo(archive_data=self._archive_data, archive_schema=self._archive_schema)
            return metainfo
