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
Contains interfaces to the archive metainfo and query.

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
"""

import numpy as np
import requests
import os.path
from urllib.parse import urlparse
from typing import Dict, List, Any

from nomad.metainfo import MSection, Quantity, SubSection
from nomad.metainfo.metainfo import MObjectMeta
from nomad import config as nomad_config
from nomad.cli.client.client import KeycloakAuthenticator


class ArchiveMetainfo:
    """
    Converts archive data in json format to the new nomad metainfo model
    Arguments:
        archive_data: the archive data in json format
    """
    def __init__(self, archive_data: List[Dict[str, Any]]):
        self._archive_data = archive_data
        self.metainfo = None
        self._metacls = None
        self._calcs: Dict[str, MSection] = {}
        self._calc_ids: List[str] = []
        self._base_metacls = None
        self._base_data = None
        self._prefix = 'calc'
        self._init_calcs()

    def _init_calcs(self):
        for calc in self._archive_data:
            calc_id = list(calc.keys())[0]
            data = calc[calc_id]
            self._calc_ids.append(calc_id)
            self._calcs[calc_id] = self._build_meta_cls(data, calc_id).m_from_dict(data)

    def __getitem__(self, key):
        if isinstance(key, str):
            calc = self._calcs.get(key, None)
            if calc:
                calc.calc_id = key
            return calc
        elif isinstance(key, int):
            calc_id = self._calc_ids[key]
            calc = self._calcs[calc_id]
            calc.calc_id = calc_id
            return calc
        else:
            calc_ids = self._calc_ids[key]
            calcs = []
            for calc_id in calc_ids:
                calc = self._calcs[calc_id]
                calc.calc_id = calc_id
                calcs.append(calc)
            return calcs

    def __len__(self):
        return len(self._calcs)

    def __iter__(self):
        self._n = -1
        return self

    def __next__(self):
        self._n += 1
        if self._n >= len(self):
            raise StopIteration
        calc = list(self._calcs.values())[self._n]
        calc.calc_id = list(self._calcs.keys())[self._n]
        return calc

    @property
    def calcs(self):
        """
        Calculations in metainfo form which can be actively queried by using the get
        functionality and providing a schema
        """
        if not self._calcs:
            self._init_calcs()
        for calc_id, calc in self._calcs.items():
            calc.calc_id = calc_id
            yield calc

    @property
    def base_data(self):
        if self._base_data is None:
            calc_id = self._calc_ids[0]
            self._base_data = self._archive_data[calc_id]
        return self._base_data

    @property
    def base_metacls(self):
        """
        The base metaclass to apply a calculation
        """
        if self._base_metacls is None:
            name = self._prefix
            self._base_metacls = self._build_meta_cls(self.base_data, name)
        return self._base_metacls

    def _get_dtype(self, data):
        if isinstance(data, np.ndarray):
            if len(data) == 0:
                dtype = int
            else:
                dtype = self._get_dtype(data[0])
        else:
            dtype = type(data)
        return dtype

    def _to_meta_obj(self, content):
        if isinstance(content, Quantity):
            return content
        if isinstance(content, MObjectMeta):
            return SubSection(sub_section=content, repeats=content.repeats)
        else:
            if isinstance(content, list):
                content = np.array(content)
            dtype = self._get_dtype(content)
            if isinstance(content, np.ndarray):
                dtype = np.dtype(dtype)
                shape = np.shape(content)
                return Quantity(type=dtype, shape=shape)
            else:
                return Quantity(type=dtype)

    def _create_section(self, name, contents):
        section = type(name.title(), (MSection,), contents)
        return section

    def _build_meta_cls(self, data, name, return_section=True):
        if isinstance(data, dict):
            contents = {}
            for key, val in data.items():
                content = self._build_meta_cls(val, key, True)
                content = self._to_meta_obj(content)
                contents[key] = content
            if return_section:
                section = self._create_section(name, contents)
                section.repeats = False
                return section
            else:
                return contents

        elif isinstance(data, list):
            if not data:
                return self._to_meta_obj(data)
            if not isinstance(data[0], dict):
                return self._to_meta_obj(data)
            contents = {}
            for i in range(len(data)):
                content = self._build_meta_cls(data[i], name, False)
                contents.update(content)
            section = self._create_section(name, contents)
            section.repeats = True
            return section

        else:
            return self._to_meta_obj(data)


class ArchiveQuery:
    def __init__(self, *args, **kwargs):
        self.archive_data = []
        self._scroll_id = None
        self._page = None
        self._query_params = {}
        if args:
            self._query_params = args[0]
        if kwargs:
            self._query_params.update(kwargs)
        self._max_n_pages = self._query_params.pop('max_n_pages', 100000)
        self._authentication = self._query_params.pop('authentication', None)
        self._url = self._query_params.pop('url', None)
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
            host = urlparse(nomad_config.client.url).netloc.split(':')[0]
            self._authentication = KeycloakAuthenticator(
                host=host,
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
        url = os.path.join(nomad_config.client.url, 'archive', 'query')

        if self._scroll_id is not None:
            self._query_params['scroll']['scroll_id'] = self._scroll_id
        elif self._page is not None:
            self._query_params['pagination']['page'] = self._page + 1

        response = requests.post(url, headers=self._get_authentication(), json=self._query_params)
        if response.status_code != 200:
            raise response.raise_for_status()

        data = response.json
        if not isinstance(data, dict):
            data = data()

        results = data.get('results', [])
        self._scroll_id = data.get('scroll', {}).get('scroll_id', None)
        self._page = data.get('pagination', {}).get('page', None)

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
            self.metainfo = ArchiveMetainfo(archive_data=self.archive_data)


def query(*args, **kwargs):
    archive_query_obj = ArchiveQuery(*args, **kwargs)
    archive_query_obj.query()
    return archive_query_obj.metainfo
