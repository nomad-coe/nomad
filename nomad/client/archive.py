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

import asyncio
from asyncio import Semaphore
from typing import Any, Dict, List, Tuple
import threading

import httpx
from httpx import Timeout
from keycloak import KeycloakOpenID

from nomad import config, metainfo as mi
from nomad.datamodel import EntryArchive, ClientContext


class RunThread(threading.Thread):
    def __init__(self, func, args, kwargs):
        self.func = func
        self.args = args
        self.kwargs = kwargs
        self.result = None
        super().__init__()

    def run(self):
        self.result = asyncio.run(self.func(*self.args, **self.kwargs))


def run_async(func, *args, **kwargs):
    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None
    if loop and loop.is_running():
        # In jupyter there is already a loop running
        thread = RunThread(func, args, kwargs)
        thread.start()
        thread.join()
        return thread.result
    else:
        # Create our own loop
        return asyncio.run(func(*args, **kwargs))


def _collect(required, parent_section: mi.Section, parent_path: str = None) -> set:
    '''
    Flatten required quantities for uncoupled query
    '''

    quantities: set = set()

    if not isinstance(required, dict):
        return quantities

    for key, value in required.items():
        # some keys may index, get the exact name
        definition_name = key.split('[')[0]

        # validate definition name
        definition = parent_section.all_properties.get(definition_name)
        if definition is None:
            raise KeyError(f'{definition_name} is not a property of {parent_section}.')

        if parent_path:
            definition_name = f'{parent_path}.{definition_name}'

        quantities.add(definition_name)

        if isinstance(definition, mi.SubSection):
            quantities.update(_collect(value, definition.section_def, definition_name))
        elif isinstance(definition, mi.Quantity) and isinstance(definition.type, mi.Reference):
            next_parent_section = definition.type.target_section_def.m_resolved()
            parent_path = next_parent_section.path
            if parent_path in ['__ambiguous__', '__no_archive_path__']:
                continue
            quantities.update(_collect(value, next_parent_section, parent_path))

    return quantities


class ArchiveQuery:
    '''
    The async implementation works well with a large number of uploads, each has a small number of
    entries.

    Authentication is created by using valid `username` and `password`.
    If any is invalid, authenticated access is not available.

    Params:
        owner (str): ownership scope
        query (dict): query
        required (dict): required properties
        url (str): server url, if not specified, the default official NOMAD one is used
        after (str): specify the starting upload id to query, if users have knowledge of uploads,
            they may wish to start from a specific upload, default: ''
        results_max (int): maximum results to query, default: 1000
        page_size (int): size of page in each query, default: 100
        username (str): username for authenticated access, default: ''
        password (str): password for authenticated access, default: ''
        retry (int): number of retry when fetching uploads, default: 4
        sleep_time (float): sleep time for retry, default: 1.
    '''

    def __init__(
            self, owner: str = 'visible', query: dict = None, required: dict = None,
            url: str = None, after: str = None, results_max: int = 1000, page_size: int = 10,
            username: str = None, password: str = None, retry: int = 4, sleep_time: float = 4,
            from_api: bool = False):
        self._owner: str = owner
        self._required = required if required else dict(run='*')
        self._query_list: List[dict] = []
        if query:
            self._query_list.append(query)
        self._query_list.append({'quantities': list(_collect(self._required, EntryArchive.m_def))})
        self._url: str = url if url else config.client.url + '/v1'
        self._after: str = after
        self._results_max: int = results_max if results_max > 0 else 1000
        self._page_size: int = page_size if page_size > 0 else 10
        self._retry: int = retry if retry >= 0 else 4
        self._sleep_time: float = sleep_time if sleep_time > 0. else 1.

        from nomad.client import Auth
        self._auth = Auth(user=username, password=password, from_api=from_api)

        self._oidc = KeycloakOpenID(
            server_url=config.keycloak.public_server_url, realm_name=config.keycloak.realm_name,
            client_id=config.keycloak.client_id)

        # local data storage
        self._uploads: List[Tuple[str, int]] = []
        self._current_after: str = self._after
        self._current_results: int = 0

        # check if url has the form of http(s)://<hostname>/api/v1
        # http://nomad-lab.eu/prod/v1/api/v1
        if self._url.endswith('/'):
            self._url = self._url[:-1]
        if not self._url.endswith('/api/v1'):
            self._url += '/v1' if self._url.endswith('/api') else '/api/v1'

    @property
    def _query(self) -> dict:
        return {'and': self._query_list}

    @property
    def _fetch_url(self) -> str:
        return f'{self._url}/entries/query'

    @property
    def _download_url(self) -> str:
        return f'{self._url}/entries/archive/query'

    @property
    def _fetch_request(self) -> dict:
        '''
        Generate fetch request.
        '''

        request = {
            'owner': self._owner,
            'query': self._query,
            'pagination': {
                'page_size': 0
            },
            'aggregations': {
                'uploads': {
                    'terms': {
                        'quantity': 'upload_id',
                        'pagination': {
                            'page_size': self._page_size,
                            'page_after_value': self._current_after
                        }
                    }
                }
            }
        }

        # print(f'Current request: {request}')

        return request

    def _download_request(self, upload_id: str, upload_count: int) -> dict:
        '''
        Generate download request.
        '''

        request: Dict[str, Any] = dict(owner=self._owner, required=self._required)
        request['query'] = {'and': []}
        for t_list in self._query_list:
            request['query']['and'].append(t_list)
        request['query']['and'].append({'upload_id': upload_id})
        request.setdefault('pagination', {'page_size': upload_count})

        # print(f'Current request: {request}')

        return request

    def clear(self):
        '''
        Clear all fetched and downloaded data. Users can then call .fetch() and .download() again.
        '''

        self._uploads: List[Tuple[str, Any]] = []  # (id, count)
        self._current_after = self._after
        self._current_results = 0

    async def _fetch_async(self, number: int = 0) -> int:
        '''
        There is no need to perform fetching asynchronously as the required number of uploads
        depends on previous queries.

        It is just a wrapper to avoid the `requests` library.

        Params:
            number (int): approx. number of **entries** to be fetched

        Returns:
            The number of entries fetched
        '''

        # if maximum number of entries have been previously fetched
        # not going to fetch more entries
        if self._current_results >= self._results_max:
            return 0

        # get all entries at once
        if number == 0:
            number = self._results_max

        num_retry: int = 0
        num_entry: int = 0

        async with httpx.AsyncClient(timeout=Timeout(timeout=300)) as session:
            while True:
                response = await session.post(
                    self._fetch_url, json=self._fetch_request, headers=self._auth.headers())

                if response.status_code >= 400:
                    if response.status_code < 500:
                        response_json = response.json()
                        reason = response_json.get("description") or response_json.get(
                            "detail") or "unknown reason"
                        raise ValueError(f'Server returns {response.status_code}: {reason}')
                    if response.status_code in [500, 502, 504]:
                        if num_retry > self._retry:
                            print('Maximum retry reached.')
                            break
                        else:
                            print(f'Retrying in {self._sleep_time} seconds...')
                            await asyncio.sleep(self._sleep_time)
                            num_retry += 1
                            continue

                response_json = response.json()

                data = response_json['aggregations']['uploads']['terms']
                header = [(bucket['value'], int(bucket['count'])) for bucket in data['data']]
                self._current_after = data['pagination'].get('next_page_after_value', None)

                current_size: int = sum([count for _, count in header])

                # no more entries
                if current_size == 0:
                    break

                if self._current_results + current_size >= self._results_max:
                    # current query has sufficient entries to exceed the limit
                    for upload_id, count in header:
                        self._current_results += count
                        num_entry += count
                        # required number of entries have been acquired
                        if self._current_results >= self._results_max:
                            entry_difference = self._current_results - self._results_max
                            self._current_results -= entry_difference
                            num_entry -= entry_difference
                            self._uploads.append((upload_id, count - entry_difference))
                            self._current_after = upload_id
                            break
                        else:
                            self._uploads.append((upload_id, count))
                    break
                else:
                    # current query should be added
                    num_entry += current_size
                    self._current_results += current_size
                    self._uploads.extend(header)

                    # if exceeds the required number, exit
                    # `self._current_after` is automatically set
                    if num_entry >= number:
                        break

        print(f'{num_entry} entries are qualified and added to the download list.')

        return num_entry

    async def _download_async(self, number: int = 0) -> List[EntryArchive]:
        '''
        Download required entries asynchronously.

        Params:
            number (int): number of **entries** to download

        Returns:
            A list of EntryArchive
        '''

        num_entry: int = 0
        num_upload: int = 0
        for _, count in self._uploads:
            num_entry += count
            num_upload += 1
            if num_entry >= number:
                break

        semaphore = Semaphore(30)

        async with httpx.AsyncClient(timeout=Timeout(timeout=300)) as session:
            tasks = [asyncio.create_task(
                self._acquire(
                    upload, session, semaphore)) for upload in self._uploads[:num_upload]]
            results = await asyncio.gather(*tasks)

        # flatten 2D list
        return [result for sub_results in results for result in sub_results]

    async def _acquire(self, upload: Tuple[str, int], session, semaphore) -> List[EntryArchive]:
        '''
        Perform the download task.

        Params:
            upload (Tuple[str, int]): upload

            session (httpx.AsyncClient): httpx client

            semaphore (asyncio.Semaphore): semaphore

        Returns:
            A list of EntryArchive
        '''

        request = self._download_request(upload[0], upload[1])

        async with semaphore:
            response = await session.post(self._download_url, json=request, headers=self._auth.headers())
            if response.status_code >= 400:
                print(
                    f'Request with upload id {upload[0]} returns {response.status_code},'
                    f' will retry in the next download call...')
                return []

            # successfully downloaded data
            # extract and remove the corresponding id in the list
            response_json = response.json()

            self._uploads.remove(upload)

            context = ClientContext(self._url, upload_id=upload[0], auth=self._auth)
            result = [EntryArchive.m_from_dict(
                result['archive'], m_context=context) for result in response_json['data']]

            if not result:
                print(f'No result returned for id {upload[0]}, is the query proper?')

            return result

    def fetch(self, number: int = 0) -> int:
        '''
        Fetch uploads from remote.

        Params:
            number (int): number of **entries** to fetch

        Returns:
            The number of entries fetched
        '''

        print('Fetching remote uploads...')

        return run_async(self._fetch_async, number)

    def download(self, number: int = 0) -> List[EntryArchive]:
        '''
        Download fetched entries from remote.
        Automatically call .fetch() if not fetched.

        Params:
            number (int): number of **entries** to download at a single time

        Returns:
            A list of downloaded EntryArchive
        '''

        pending_size: int = sum([count for _, count in self._uploads])

        # download all at once
        if number == 0:
            number = pending_size

        if number == 0:
            # empty list, fetch as many as possible
            number = self.fetch()
        elif pending_size < number:
            # if not sufficient fetched entries, fetch first
            self.fetch(number - pending_size)

        print('Downloading required data...')

        return run_async(self._download_async, number)

    async def async_fetch(self, number: int = 0) -> int:
        '''
        Asynchronous interface for use in a running event loop.
        '''

        print('Fetching remote uploads...')

        return await self._fetch_async(number)

    async def async_download(self, number: int = 0) -> List[EntryArchive]:
        '''
        Asynchronous interface for use in a running event loop.
        '''

        pending_size: int = sum([count for _, count in self._uploads])

        # download all at once
        if number == 0:
            number = pending_size

        if number == 0:
            # empty list, fetch as many as possible
            number = await self.async_fetch()
        elif pending_size < number:
            # if not sufficient fetched entries, fetch first
            await self.async_fetch(number - pending_size)

        print('Downloading required data...')

        return await self._download_async(number)

    def upload_list(self) -> List[Tuple[str, int]]:
        return self._uploads
