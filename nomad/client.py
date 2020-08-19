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
.. _install-client:

Install the NOMAD client library
________________________________

The NOMAD client library is a Python module (part of the nomad Python package) that
allows to access the NOMAD archive to retrieve and analyse (large amounts) of NOMAD's
archive data. It allows to use queries to filter for desired entries, bulk download
the required parts of the respective archives, and navigate the results using NOMAD's
metainfo Python API.

To install the NOMAD Python package, you can use ``pip install`` to install our
source distribution

.. parsed-literal::

    pip install nomad-lab


First example
_____________

.. literalinclude:: ../examples/client.py
    :language: python

This script should yield a result like this:

.. code::

    Number queries entries: 7628
    Number of entries loaded in the last api call: 10
    Bytes loaded in the last api call: 118048
    Bytes loaded from this query: 118048
    Number of downloaded entries: 10
    Number of made api calls: 1

    Cd2O2: energy -11467.827149010665 hartree
    Sr2O2: energy -6551.45699684026 hartree
    Sr2O2: energy -6551.461104765451 hartree
    Be2O2: energy -178.6990610734937 hartree
    Ca2O2: energy -1510.3938165430286 hartree
    Ca2O2: energy -1510.3937761449583 hartree
    Ba2O2: energy -16684.667362890417 hartree
    Mg2O2: energy -548.9736595672932 hartree
    Mg2O2: energy -548.9724185656775 hartree
    Ca2O2: energy -1510.3908614326358 hartree

Let's discuss the different elements here. First, we have a set of imports. The NOMAD source
codes comes with various sub-modules. The `client` module contains everything related
to what is described here; the `metainfo` is the Python interface to NOMAD's common
archive data format and its data type definitions; the `config` module simply contains
configuration values (like the URL to the NOMAD API).

Next, we create an :class:`ArchiveQuery` instance. This object will be responsible for talking
to NOMAD's API for us in a transparent and lazy manner. This means, it will not download
all data right away, but do so when we are actually iterating through the results.

The archive query takes several parameters:

- The ``query`` is a dictionary of search criteria. The query is used to filter all of NOMAD's
  entry down to a set of desired entries. You can use NOMAD's GUI to create queries and
  copy their Python equivalent with the ``<>``-code button on the result list.
- The ``required`` part, allows to specify what parts of the archive should be downloaded.
  Leave it out to download the whole archives. Based on NOMAD's Metainfo (the 'schema' of
  all archives), you can determine what sections to include and which to leave out. Here,
  we are interested in the first run (usually entries only have one run) and the first
  calculation result.
- With the optional ``per_page`` you can determine, how many results are downloaded at
  a time. For bulk downloading many results, we recommend ~100. If you are just interested
  in the first results a lower number might increase performance.
- With the optional ``max``, we limit the maximum amount of entries that are downloaded,
  just to avoid accidentely iterating through a result set of unknown and potentially large
  size.

When you print the archive query object, you will get some basic statistics about the
query and downloaded data.

The archive query object can be treated as a Python list-like. You use indices and ranges
to select results. Here we iterate through a slice and print the calculated energies
from the first calculation of the entries. Each result is a Python object with attributes
governed by the NOMAD Metainfo. Quantities yield numbers, string, or numpy arrays, while
sub-sections return lists of further objects. Here we navigate the sections ``section_run`` and
sub-section ``section_system`` to access the quantity ``energy_total``. This quantity is a
number with an attached unit (Joule), which can be converted to something else (e.g. Hartree).

The create query object keeps all results in memory. Keep this in mind, when you are
accessing a large amount of query results. You should use :func:`ArchiveQuery.clear`
to remove unnecessary results.

The NOMAD Metainfo
__________________

You can imagine the NOMAD Metainfo as a complex schema for hiearchically organized scientific
data. In this sense, the NOMAD Metainfo is a set of data type definitions. These definitions
then govern how the archive for an data entry in NOMAD might look like. You can browse the
hierarchy of definitions in our `Metainfo browser <../metainfo>`_.

Be aware, that the definitions entail everything that an entry could possibly contain, but
not all entries contain all sections and all quantities. What an entry contains depends
on the information that the respective uploaded data contained, what could be extracted,
and of course what was calculated in the first place. To see what the archive of an concrete
entry looks like, you can use the `search interface <../search>`_, select an entry from the
list fo search results, and click on the *Archive* tab.

To *see inside* an archive object in Python, you can use :func:`nomad.metainfo.MSection.m_to_dict`
which is provided by all archive objects. This will convert a (part of an) archive into a
regular, JSON-serializable Python dictionary.

For more details on the metainfo Python interface, consult the `metainfo documentation <metainfo.html>`_.

The ArchiveQuery class
______________________

.. autoclass:: ArchiveQuery

Working with private data
_________________________

Public NOMAD data can be accessed without any authentication; everyone can use our API
without the need for an account or login. However, if you want to work with your own
data that is not yet published, or embargoed data was shared with you, you need to
authenticate before accessing this data. Otherwise, you will simply not find it with
your queries. To authenticate simply provide your NOMAD username and password to the
:class:`ArchiveQuery` constructor.

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


class QueryError(Exception):
    pass


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
        type=int, default=0,
        description='Number queries entries')

    last_response_nentries = mi.Quantity(
        type=int, default=0,
        description='Number of entries loaded in the last api call')

    last_response_data_size = mi.Quantity(
        type=int, unit=mi.units.bytes, default=0,
        description='Bytes loaded in the last api call')

    loaded_data_size = mi.Quantity(
        type=int, unit=mi.units.bytes, default=0,
        description='Bytes loaded from this query')

    loaded_nentries = mi.Quantity(
        type=int, default=0,
        description='Number of downloaded entries')

    napi_calls = mi.Quantity(
        type=int, default=0,
        description='Number of made api calls')

    def __repr__(self):
        out = StringIO()
        for quantity in self.m_def.all_quantities.values():
            out.write('%s: %s\n' % (quantity.description, self.m_get(quantity)))

        return out.getvalue()


class ArchiveQuery(collections.abc.Sequence):
    '''
    Object of this class represent a query on the NOMAD Archive. It is solely configured
    through its constructor. After creation, it implements the
    Python ``Sequence`` interface and therefore acts as a sequence of query results.

    Not all results are downloaded at once, expect that this class will continuesly pull
    results from the API, while you access or iterate to the far side of the result list.

    Attributes:
        query: A dictionary of search parameters. Consult the search API to get a
            comprehensive list of parameters.
        required: A potentially nested dictionary of sections to retrieve.
        url: Optional, override the default NOMAD API url.
        username: Optional, allows authenticated access.
        password: Optional, allows authenticated access.
        per_page: Determine how many results are downloaded per page (or scroll window).
            Default is 10.
        max: Optionally determine the maximum amount of downloaded archives. The iteration
            will stop even if more results are available. Default is unlimited.
        raise_errors: There situations where archives for certain entries are unavailable.
            If set to True, this cases will raise an Exception. Otherwise, the entries
            with missing archives are simply skipped (default).
        authentication: Optionally provide detailed authentication information. Usually,
            providing ``username`` and ``password`` should suffice.
    '''
    def __init__(
            self,
            query: dict = None, required: dict = None,
            url: str = None, username: str = None, password: str = None,
            per_page: int = 10, max: int = None,
            raise_errors: bool = False,
            authentication: Union[Dict[str, str], KeycloakAuthenticator] = None):

        self._after = None
        self.page = 1
        self.per_page = per_page
        self.max = max

        self.query: Dict[str, Any] = {
            'query': {},
            'raise_errors': raise_errors
        }
        if query is not None:
            self.query['query'].update(query)
        if required is not None:
            self.query['query_schema'] = required
            # We try to add all required properties to the query to ensure that only
            # results with those properties are returned.
            section_run_key = next(key for key in required if key.split('[')[0] == 'section_run')
            if section_run_key is not None:
                # add all quantities in required to the query part
                quantities = {'section_run'}
                stack = []
                section_run = required[section_run_key]
                if isinstance(section_run, dict):
                    stack.append(section_run)
                while len(stack) > 0:
                    required_dict = stack.pop()
                    for key, value in required_dict.items():
                        if isinstance(value, dict):
                            stack.append(value)
                        quantities.add(key.split('[')[0])
                self.query['query'].setdefault('dft.quantities', []).extend(quantities)
                self.query['query']['domain'] = 'dft'

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
        '''
        The authentication information that is used, if username or password were
        provided.
        '''
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
        '''
        Calls the API to retrieve the next set of results. Is automatically called, if
        not yet downloaded entries are accessed.
        '''
        url = '%s/%s/%s' % (self.url, 'archive', 'query')

        aggregation = self.query.setdefault('aggregation', {'per_page': self.per_page})
        if self._after is not None:
            aggregation['after'] = self._after

        response = requests.post(url, headers=self.authentication, json=self.query)
        if response.status_code != 200:
            if response.status_code == 400:
                message = response.json().get('message')
                errors = response.json().get('errors')
                if message:
                    raise QueryError('%s: %s' % (message, errors))

                raise QueryError('The query is invalid for unknown reasons.')

            raise response.raise_for_status()

        data = response.json
        if not isinstance(data, dict):
            data = data()

        aggregation = data['aggregation']
        self._after = aggregation.get('after')
        self._total = aggregation['total']

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

        if self._after is None:
            # there are no more search results, we need to avoid further calls
            self._capped_total = len(self._results)
            self._total = len(self._results)

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
        ''' The total ammount of search results. '''
        if self._total == -1:
            self.call_api()

        return self._total

    @property
    def statistics(self):
        ''' A metainfo object with a basic set of query statistics. '''
        if self._total == -1:
            self.call_api()

        return self._statistics

    def clear(self, index: int = None):
        '''
        Remove caches results. The results are replaced with None in this object. If you
        keep references to the results elsewhere, the garbage collection might not catch
        those.

        Arguments:
            index: Remove all results upto and including the giving index. Default is to
                remove all results.
        '''
        for i, _ in enumerate(self._results[:index]):
            print(i)
            self._results[i] = None


def query_archive(*args, **kwargs):
    return ArchiveQuery(*args, **kwargs)


if __name__ == '__main__':
    run = query_archive()[1]
    run.section_system[1].atom_labels
