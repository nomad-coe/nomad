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
import collections.abc
import requests
from urllib.parse import urlparse
from bravado import requests_client as bravado_requests_client
import time
from keycloak import KeycloakOpenID
from io import StringIO

from nomad import config
from nomad import metainfo as mi
from nomad.datamodel import EntryArchive

# TODO this import is necessary to load all metainfo defintions that the parsers are using
from nomad import parsing  # pylint: disable=unused-import


class KeycloakAuthenticator(bravado_requests_client.Authenticator):
    def __init__(self, host, user, password, **kwargs):
        super().__init__(host=host)
        self.user = user
        self.password = password
        self.token = None
        self.__oidc = KeycloakOpenID(**kwargs)

    def apply(self, request=None):
        if self.token is None:
            self.token = self.__oidc.token(username=self.user, password=self.password)
            self.token['time'] = time.time()
        elif self.token['expires_in'] < int(time.time()) - self.token['time'] + 10:
            try:
                self.token = self.__oidc.refresh_token(self.token['refresh_token'])
                self.token['time'] = time.time()
            except Exception:
                self.token = self.__oidc.token(username=self.user, password=self.password)
                self.token['time'] = time.time()

        if request:
            request.headers.setdefault('Authorization', 'Bearer %s' % self.token['access_token'])
            return request
        else:
            return dict(Authorization='Bearer %s' % self.token['access_token'])


class ApiStatistics(mi.MSection):

    nentries = mi.Quantity(
        type=int,
        description='Number queries entries')

    last_response_nentries = mi.Quantity(
        type=int,
        description='Number of entries loaded in the last api call')

    last_response_data_size = mi.Quantity(
        type=int, unit=mi.units.bytes,
        description='Bytes loaded in the last api call')

    loaded_data_size = mi.Quantity(
        type=int, unit=mi.units.bytes,
        description='Bytes loaded from this query')

    loaded_nentries = mi.Quantity(
        type=int,
        description='Number of downloaded entries')

    napi_calls = mi.Quantity(
        type=int,
        description='Number of made api calls')

    def __repr__(self):
        out = StringIO()
        for quantity in self.m_def.all_quantities.values():
            out.write('%s: %s\n' % (quantity.description, self.m_get(quantity)))

        return out.getvalue()


class ArchiveQuery(collections.abc.Sequence):
    def __init__(
            self,
            query: dict = None, required: dict = None,
            url: str = None, username: str = None, password: str = None,
            scroll: bool = False, per_page: int = 10, max: int = None,
            authentication: Union[Dict[str, str], KeycloakAuthenticator] = None):

        self.scroll = scroll
        self._scroll_id = None
        self.page = 1
        self.per_page = per_page
        self.max = max

        self.query: Dict[str, Any] = {
            'query': {}
        }
        if query is not None:
            self.query['query'].update(query)
        if required is not None:
            self.query['query_schema'] = required

        self.password = password
        self.username = username
        self.url = config.client.url if url is None else url
        self._authentication = authentication

        self._total = -1
        self._capped_total = -1
        self._results: List[dict] = []
        self._statistics = ApiStatistics()

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
            self.query.setdefault('pagination', {}).update(
                page=self.page, per_page=self.per_page)

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
            self.page = pagination['page'] + 1

        if self.max is not None:
            self._capped_total = min(self.max, self._total)
        else:
            self._capped_total = self._total

        results = data.get('results', [])

        for result in results:
            archive = EntryArchive.m_from_dict(result['archive'])

            self._results.append(archive)

        try:
            data_size = len(response.content)
            self._statistics.last_response_data_size = data_size
            self._statistics.loaded_data_size += data_size
            self._statistics.nentries = self._total
            self._statistics.last_response_nentries = len(results)
            self._statistics.loaded_nentries = len(self._results)
            self._statistics.napi_calls += 1
        except Exception:
            # fails in test due to mocked requests library
            pass

    def __repr__(self):
        if self._total == -1:
            self.call_api()

        return str(self._statistics)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return [self[i] for i in range(*key.indices(len(self)))]

        if key >= self.__len__():
            raise IndexError()

        while len(self._results) < key + 1:
            self.call_api()

        return self._results[key]

    def __len__(self):  # pylint: disable=invalid-length-returned
        if self._capped_total == -1:
            self.call_api()

        return self._capped_total

    @property
    def total(self):
        if self._total == -1:
            self.call_api()

        return self._total

    @property
    def statistics(self):
        if self._total == -1:
            self.call_api()

        return self._statistics


def query_archive(*args, **kwargs):
    return ArchiveQuery(*args, **kwargs)


if __name__ == '__main__':
    run = query_archive()[1]
    run.section_system[1].atom_labels
