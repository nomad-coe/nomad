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

'''
Contains the Python client side library to access the NOMAD archive.

# TODO
In module ``ArchiveMetainfo``, the data is provided either from raw
json data or as a filename of an existing msgpack database. The metainfo
can then queried by providing a schema.

.. code-block: python
    am = ArchiveMetainfo(archive_data)
    for calc in am.calcs:
        c.section_run.section_single_configuration_calculation[0]({'energy_total':None})

The ArchiveQuery enables a query interface to the archive data. A dict of query parameters
and a query schema similar to the archive json format can be provided to filter archive.

.. code-block: python
    q = ArchiveQuery({'atoms':'Si'})
    metainfo = q.query()
    for c in metainfo.calcs:
        print(c.section_run.section_single_configuration_calculation[0]({'energy_total':'*'}))
'''

from typing import Dict, Union, Any, List
from collections import Sequence
import requests
from urllib.parse import urlparse

from nomad import config
from nomad.cli.client.client import KeycloakAuthenticator
from nomad.datamodel import EntryArchive

# TODO this import is necessary to load all metainfo defintions that the parsers are using
from nomad import parsing  # pylint: disable=unused-import


class ArchiveQuery(Sequence):
    def __init__(
            self,
            query: dict = None, query_schema: dict = None,
            url: str = None, username: str = None, password: str = None,
            scroll: bool = False,
            authentication: Union[Dict[str, str], KeycloakAuthenticator] = None, **kwargs):

        self.scroll = scroll
        self._scroll_id = None
        self._page = 1

        self.query: Dict[str, Any] = {
            'query': {}
        }
        if query is not None:
            self.query['query'].update(query)
        if query_schema is not None:
            self.query['query_schema'] = query_schema
        self.query['query'].update(kwargs)

        self.password = password
        self.username = username
        self.url = config.client.url if url is None else url
        self._authentication = authentication

        self._total = -1
        self._results: List[dict] = []

    @property
    def authentication(self):
        if self._authentication is None and self.username is not None and self.password is not None:
            host = urlparse(self.url).netloc.split(':')[0]
            self._authentication = KeycloakAuthenticator(
                host=host,
                user=self.username,
                password=self.password,
                server_url=config.keycloak.server_url,
                realm_name=config.keycloak.realm_name,
                client_id=config.keycloak.client_id)

        if isinstance(self._authentication, KeycloakAuthenticator):
            return self._authentication.apply()

        else:
            return self._authentication

    def call_api(self):
        url = '%s/%s/%s' % (self.url, 'archive', 'query')

        if self.scroll:
            scroll_config = self.query.setdefault('scroll', {'scroll': True})
            if self._scroll_id is not None:
                scroll_config['scroll_id'] = self._scroll_id

        else:
            self.query.setdefault('pagination', {})['page'] = self._page

        response = requests.post(url, headers=self.authentication, json=self.query)
        if response.status_code != 200:
            raise response.raise_for_status()

        data = response.json
        if not isinstance(data, dict):
            data = data()

        if self.scroll:
            scroll = data['scroll']
            self._scroll_id = scroll['scroll_id']
            self._total = scroll['total']

        else:
            pagination = data['pagination']
            self._total = pagination['total']
            self._page = pagination['page'] + 1

        results = data.get('results', [])

        for result in results:
            archive = EntryArchive.m_from_dict(result['archive'])

            self._results.append(archive)

    def __getitem__(self, key):
        if key >= self.__len__():
            raise IndexError()

        while len(self._results) < key:
            self.call_api()

        return self._results[key]

    def __len__(self):  # pylint: disable=invalid-length-returned
        if self._total == -1:
            self.call_api()

        return self._total


def query_archive(*args, **kwargs):
    return ArchiveQuery(*args, **kwargs)


if __name__ == '__main__':
    run = query_archive()[1]
    run.section_system[1].atom_labels
