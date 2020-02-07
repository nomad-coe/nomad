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

from nomad import config as nomad_config
from nomad.archive_library.metainfo import ArchiveMetainfo
from nomad.cli.client.client import KeycloakAuthenticator


class ArchiveQuery:
    def __init__(self, *args, **kwargs):
        self._archive_path = 'archive'
        self._query_path = 'query'
        self.archive_data = []
        self._scroll_id = None
        self._page = None
        self._query_params = {}
        if args:
            self._query_params = args[0]
        if kwargs:
            self._query_params.update(kwargs)
        self._archive_schema = self._query_params.pop('archive_data', None)
        if not isinstance(self._archive_schema, list):
            self._archive_schema = [self._archive_schema]
        self._max_n_pages = self._query_params.pop('max_n_pages', 100000)
        self._authentication = self._query_params.pop('authentication', None)
        self._url = self._query_params.pop('rul', None)
        self._user = self._query_params.pop('user', None)
        self._password = self._query_params.pop('password', None)
        if self._url:
            nomad_config.client.url = self._url
        if self._user:
            nomad_config.client.user = self._user
        if self._password:
            nomad_config.client.password = self._password

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

    def _get_authentication(self):
        if self._authentication is None:
            self._authentication = KeycloakAuthenticator(
                user=nomad_config.client.user,
                password=nomad_config.client.password,
                server_url=nomad_config.keycloak.server_external_url,
                realm_name=nomad_config.keycloak.realm_name,
                client_id=nomad_config.keycloak.public_client_id)
        if isinstance(self._authentication, KeycloakAuthenticator):
            return self._authentication.apply()
        else:
            return self._authentication

    def _api_query(self):
        url = os.path.join(nomad_config.client.url, self._archive_path, self._query_path)
        data = self._query_params
        data['results'] = self._archive_schema

        if self._page is not None:
            # increment the page number
            self._set_value('page', self._page + 1, data)
        if self._scroll_id is not None:
            self._set_value('scroll_id', self._scroll_id, data)

        response = requests.post(url, headers=self._get_authentication(), json=data)
        if response.status_code != 200:
            raise Exception('Query returned %s' % response.status_code)

        data = response.json
        if not isinstance(data, dict):
            data = data()

        results = data.get('results', [])
        scroll = data.get('Scroll', None)
        if scroll:
            self._scroll_id = scroll.get('scroll_id', None)
        pagination = data.get('Pagination', None)
        if pagination:
            self._page = pagination.get('page', None)

        return results

    def _get_archive_data(self):
        n_page = 0
        while True:
            results = self._api_query()
            self.archive_data += results
            n_page += 1
            if n_page >= self._max_n_pages:
                break
            if len(results) == 0:
                break

    def query(self):
        self._get_archive_data()
        if self.archive_data:
            self.metainfo = ArchiveMetainfo(archive_data=self.archive_data, archive_schema='*')
