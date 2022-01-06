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

'''
Retrieve and analyse (large amounts) of NOMAD's
archive data. This allows to use queries to filter for desired entries, bulk download
the required parts of the respective archives, and navigate the results using NOMAD's
metainfo Python API.
'''

from typing import Dict, Any, List
import collections.abc
import requests
from io import StringIO
import math
import multiprocessing
import json

from nomad import config
from nomad import metainfo as mi
from nomad.datamodel import EntryArchive

# TODO this import is necessary to load all metainfo defintions that the parsers are using
from nomad import parsing  # pylint: disable=unused-import

from .api import Auth


class QueryError(Exception):
    pass


class ApiStatistics(mi.MSection):
    total = mi.Quantity(
        type=int, default=0,
        description='Total number of entries that fulfil the query')

    nentries = mi.Quantity(
        type=int, default=0,
        description='Number queried entries')

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


class ProcState:
    '''
    A basic pickable data-class that holds the state of one parallel running
    processes that loads archive API data.
    '''
    def __init__(self, archive_query: 'ArchiveQuery'):
        self.url = archive_query.url
        self.query_and_list = archive_query.query_and_list
        self.request: Dict[str, Any] = dict(
            owner=archive_query.owner,
            required=archive_query.required)
        self.per_page = archive_query.per_page
        self.authentication = archive_query.authentication

        self.upload_ids: List[str] = []
        self.nentries = 0
        self.total = None
        self.after = None
        self.results = None
        self.error: Exception = None
        self.data_size = 0


def _run_proc(proc_state: ProcState) -> ProcState:
    '''
    The main function for a process that retrieves data from the archive API based
    on its state. Will create a new state. Otherwise it is completely stateless.
    '''
    try:
        url = f'{proc_state.url}/v1/entries/archive/query'

        # create the query
        query_and_list = list(proc_state.query_and_list)
        query_and_list.append({'upload_id:any': proc_state.upload_ids})
        proc_state.request['query'] = {'and': query_and_list}
        proc_state.request.setdefault('pagination', {'page_size': proc_state.per_page})
        if proc_state.after is not None:
            proc_state.request['pagination']['page_after_value'] = proc_state.after

        # run the query
        response = requests.post(url, auth=proc_state.authentication, json=proc_state.request)
        response_json = response.json()
        if response.status_code != 200:
            if response.status_code == 422:
                raise QueryError(f'The query is invalid (422): {response_json["detail"]}')
            if response.status_code == 400:
                message = response_json.get('message')
                errors = response_json.get('errors')
                if message:
                    raise QueryError('%s: %s' % (message, errors))
                raise QueryError('The query is invalid for unknown reasons (400).')
            raise QueryError(
                'The query is invalid for unknown reasons (%d).' % response.status_code)

        # update the state
        proc_state.data_size += len(response.content)
        proc_state.results = response_json['data']
        proc_state.after = response_json['pagination'].get('next_page_after_value', None)
        proc_state.total = response_json['pagination']['total']
    except Exception as e:
        proc_state.error = e

    return proc_state


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
            will stop if max is surpassed even if more results are available. Default is 10.000.
            None value will set it to unlimited.
        authentication: Optionally provide detailed authentication information. Usually,
            providing ``username`` and ``password`` should suffice.
        parallel: Number of processes to use to retrieve data in parallel. Only data
            from different uploads can be retrieved in parallel. Default is 1. The
            argument ``per_page`` will refer to archived retrieved in one process per
            call.
    '''
    def __init__(
            self, owner: str = 'visible',
            query: dict = None, required: dict = None,
            url: str = None, username: str = None, password: str = None,
            parallel: int = 1, per_page: int = 10, max: int = 10000,
            authentication: Auth = None, auth: Auth = None):

        self.owner = owner
        self.page = 1
        self.parallel = parallel
        self.per_page = per_page
        self.max = max

        self.query_and_list: List[dict] = []
        if query is not None:
            self.query_and_list.append(query)

        self.required = required if required is not None else dict(run='*')

        # We try to add all required properties to the query to ensure that only
        # results with those properties are returned.
        quantities = set()

        def collect(required, parent_section: mi.Section, parent_path: str = None):
            if not isinstance(required, dict):
                return

            for key, value in required.items():
                def_name = key.split('[')[0]
                definition = parent_section.all_properties.get(def_name)

                if definition is None:
                    raise KeyError(f'{def_name} is not a property of {parent_section}')

                if parent_path:
                    qualified_def_name = f'{parent_path}.{def_name}'
                else:
                    qualified_def_name = def_name
                quantities.add(qualified_def_name)

                if isinstance(definition, mi.SubSection):
                    collect(
                        value,
                        parent_section=definition.section_def,
                        parent_path=qualified_def_name)
                elif isinstance(definition, mi.Quantity) and isinstance(definition.type, mi.Reference):
                    next_parent_section = definition.type.target_section_def.m_resolved()
                    parent_path = next_parent_section.path
                    if parent_path in ['__ambiguous__', '__no_archive_path__']:
                        continue
                    collect(
                        value,
                        parent_section=next_parent_section,
                        parent_path=parent_path)

        collect(required, parent_section=EntryArchive.m_def)
        self.query_and_list.append({'quantities': list(quantities)})

        self.password = password
        self.username = username
        self.url = config.client.url if url is None else url
        self._authentication = authentication
        if not self._authentication:
            self._authentication = auth

        self._total = -1
        self._results: List[dict] = []
        self._statistics = ApiStatistics()
        self._proc_states: List[ProcState] = None

    @property
    def query(self):
        return {'and': self.query_and_list}

    @property
    def authentication(self):
        '''
        The authentication information that is used, if username or password were
        provided.
        '''
        if self._authentication is None and self.username is not None and self.password is not None:
            self._authentication = Auth(self.username, self.password)
        else:
            return self._authentication

    def _create_initial_proc_state(self):
        '''
        Does preliminary queries to the repo API to determine the distribution of queried
        entries over uploads and creates initial state for the processes that collect
        data from the archive API in parallel.
        '''
        uploads: Dict[str, Any] = dict()
        nentries = 0

        # acquire all uploads and how many entries they contain
        url = f'{self.url}/v1/entries/query'
        after: str = None

        while True:
            uploads_request = {
                'owner': self.owner,
                'query': self.query,
                'pagination': {
                    'page_size': 0
                },
                'aggregations': {
                    'uploads': {
                        'terms': {
                            'quantity': 'upload_id',
                            'pagination': {
                                'page_size': 100,
                                'page_after_value': after
                            }
                        }
                    }
                }
            }
            response = requests.post(url, json=uploads_request, auth=self.authentication)
            response_json = response.json()

            if response.status_code != 200:
                if response.status_code == 400:
                    raise Exception(response_json['description'])
                if response.status_code == 422:
                    raise Exception(response_json['detail'])
                raise Exception(
                    'Error requesting NOMAD API: HTTP %d' % response.status_code)

            total = response_json['pagination']['total']
            agg_data = response_json['aggregations']['uploads']['terms']
            after = agg_data['pagination'].get('next_page_after_value', None)
            values = {bucket['value']: bucket for bucket in agg_data['data']}

            if len(values) == 0:
                break

            uploads.update(values)
            for upload in values.values():
                nentries += upload['count']

            if self.max is not None and nentries >= self.max:
                break

        # distribute uploads to processes
        if self.parallel is None:
            self.parallel = 1

        # TODO This implements a simplified distribution, where an upload is fully
        # handled by an individual process. This works because of the high likely hood
        # that popular analysis queries (e.g. AFLOW) have results spread over
        # many uploads. In other use-cases, e.g. analysing data from an individual user,
        # this might not work well, because all entries might be contained in one
        # upload.
        self._proc_states = []
        nentries_per_proc = math.ceil(nentries / self.parallel)
        proc_state = ProcState(self)

        for upload_id, upload_data in uploads.items():
            if proc_state.nentries >= nentries_per_proc:
                self._proc_states.append(proc_state)
                proc_state = ProcState(self)

            proc_state.upload_ids.append(upload_id)
            proc_state.nentries += upload_data['count']

        self._proc_states.append(proc_state)
        self._total = nentries
        self._statistics.total = total
        self._statistics.nentries = nentries

    def call_api(self):
        '''
        Calls the API to retrieve the next set of results. Is automatically called, if
        not yet downloaded entries are accessed.
        '''
        if self._proc_states is None:
            self._create_initial_proc_state()

        # run the necessary processes
        nproc_states = len(self._proc_states)
        if nproc_states == 1:
            self._proc_states[0] = _run_proc(self._proc_states[0])
        elif nproc_states > 1:
            with multiprocessing.Pool(nproc_states) as pool:
                self._proc_states = pool.map(_run_proc, self._proc_states)
        else:
            assert False, 'archive query was not stopped before running out of things to query'

        # grab the results from the processes
        new_states: List[ProcState] = []
        self._statistics.last_response_nentries = 0
        self._statistics.last_response_data_size = 0
        for proc_state in self._proc_states:
            if proc_state.error:
                raise proc_state.error

            self._statistics.last_response_data_size += proc_state.data_size
            self._statistics.loaded_data_size += proc_state.data_size
            self._statistics.last_response_nentries += len(proc_state.results)

            self._results.extend([
                EntryArchive.m_from_dict(result['archive'])
                for result in proc_state.results])

            proc_state.results = None
            if proc_state.after is not None:
                new_states.append(proc_state)

        self._proc_states = new_states
        self._statistics.loaded_nentries = len(self._results)
        self._statistics.napi_calls += 1

        if self.max is not None and len(self._results) >= self.max:
            # artificially end the query
            self._proc_states = []

        if len(self._proc_states) == 0:
            self._total = len(self._results)

    def __repr__(self):
        if self._total == -1:
            self.call_api()

        return f'Query: {json.dumps(self.query, indent=2)}\n{self._statistics}'

    def __getitem__(self, key):
        if isinstance(key, slice):
            return [self[i] for i in range(*key.indices(len(self)))]

        if key >= self.__len__():
            raise IndexError()

        while len(self._results) < key + 1:
            self.call_api()

        return self._results[key]

    def __len__(self):  # pylint: disable=invalid-length-returned
        if self._total == -1:
            self.call_api()

        return self._total

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
    run.system[1].atoms.labels
