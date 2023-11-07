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
from __future__ import annotations

import asyncio
from asyncio import Semaphore
from typing import Any, Dict, List
import threading

from click import progressbar
from httpx import Timeout, AsyncClient
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


def _collect(required, parent_section: mi.Section = None, parent_path: str = None) -> set:
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
        definition = None
        if parent_section:
            definition = parent_section.all_properties.get(definition_name)

        if parent_path:
            definition_name = f'{parent_path}.{definition_name}'

        if not definition:
            # We have to stop here because we cannot further statically analyze the
            # required. We have to assume the remainder is based on a polymorphic subsection,
            # and we cannot know the exact subsection type.
            return quantities

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

    Setting a high value for `semaphore` may cause the server to return 500, 502, 504 errors.

    Params:
        owner (str): ownership scope
        query (dict): query
        required (dict): required properties
        url (str): server url, if not specified, the default official NOMAD one is used
        after (str): specify the starting upload id to query, if users have knowledge of uploads,
            they may wish to start from a specific upload, default: ''
        results_max (int): maximum results to query, default: 1000
        page_size (int): size of page in each query, cannot exceed the limit 10000, default: 100
        username (str): username for authenticated access, default: ''
        password (str): password for authenticated access, default: ''
        retry (int): number of retry when fetching uploads, default: 4
        sleep_time (float): sleep time for retry, default: 4.
        semaphore (int): number of concurrent downloads, this depends on server settings, default: 4
    '''

    def __init__(
            self, owner: str = 'visible', query: dict = None, required: dict = None,
            url: str = None, after: str = None, results_max: int = 1000, page_size: int = 100,
            username: str = None, password: str = None, retry: int = 4, sleep_time: float = 4.,
            from_api: bool = False, semaphore: int = 4):
        self._owner: str = owner
        self._required = required if required else dict(run='*')
        self._query_list: list[dict] = [{'quantities': list(_collect(self._required, EntryArchive.m_def))}]
        if query:
            self._query_list.append(query)
        self._url: str = url if url else config.client.url + '/v1'
        self._after: str = after
        self._results_max: int = results_max if results_max > 0 else 1000
        self._page_size: int = min(page_size, 9999) if page_size > 0 else 100
        if self._page_size > self._results_max:
            self._page_size = self._results_max
        self._retry: int = retry if retry >= 0 else 4
        self._sleep_time: float = sleep_time if sleep_time > 0. else 4.
        self._semaphore = semaphore

        from nomad.client import Auth
        self._auth = Auth(user=username, password=password, from_api=from_api)

        self._oidc = KeycloakOpenID(
            server_url=config.keycloak.public_server_url, realm_name=config.keycloak.realm_name,
            client_id=config.keycloak.client_id)

        # local data storage
        self._entries: list[tuple[str, str]] = []
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

        request: dict = {
            'owner': self._owner,
            'query': self._query,
            'pagination': {'page_size': self._page_size},
            "required": {"include": ["entry_id", "upload_id"]}
        }

        if self._current_after:
            request['pagination']['page_after_value'] = self._current_after

        # print(f'Current request: {request}')

        return request

    def _download_request(self, entry_id: str) -> dict:
        '''
        Generate download request.
        '''

        request: Dict[str, Any] = dict(owner=self._owner, required=self._required)
        request['query'] = {'and': []}
        for t_list in self._query_list:
            request['query']['and'].append(t_list)
        request['query']['and'].append({'entry_id': entry_id})
        request.setdefault('pagination', {'page_size': 1})

        # print(f'Current request: {request}')

        return request

    def clear(self):
        '''
        Clear all fetched and downloaded data. Users can then call .fetch() and .download() again.
        '''

        self._entries = []
        self._current_after = self._after
        self._current_results = 0

    async def _fetch_async(self, number: int) -> int:
        '''
        There is no need to perform fetching asynchronously as the required number of uploads
        depends on previous queries.

        It is just a wrapper to avoid the `requests` library.

        Params:
            number (int): approx. number of **entries** to be fetched

        Returns:
            The number of entries fetched
        '''

        # if the maximum number of entries has been previously fetched
        # not going to fetch more entries
        if self._current_results >= self._results_max:
            return 0

        # get all entries at once
        if number == 0:
            number = self._results_max

        num_retry: int = 0
        num_entry: int = 0

        async with AsyncClient(timeout=Timeout(timeout=300)) as session:
            while True:
                response = await session.post(
                    self._fetch_url, json=self._fetch_request, headers=self._auth.headers())

                if response.status_code >= 400:
                    if response.status_code < 500:
                        response_json = response.json()
                        reason = response_json.get("description") or response_json.get(
                            "detail") or "unknown reason"
                        raise ValueError(f'Server returns {response.status_code}: {reason}')
                    if response.status_code in (500, 502, 504):
                        if num_retry > self._retry:
                            print('Maximum retry reached.')
                            break
                        else:
                            print(f'Retrying in {self._sleep_time} seconds...')
                            await asyncio.sleep(self._sleep_time)
                            num_retry += 1
                            continue

                response_json = response.json()

                self._current_after = response_json['pagination'].get('next_page_after_value', None)

                data = [(entry['entry_id'], entry['upload_id']) for entry in response_json['data']]
                current_size: int = len(data)

                # no more entries
                if current_size == 0:
                    break

                if self._current_results + current_size > self._results_max:
                    # current query has sufficient entries to exceed the limit
                    data = data[:self._results_max - self._current_results]
                    self._current_results += len(data)
                    self._entries.extend(data)
                    break
                else:
                    # current query should be added
                    num_entry += current_size
                    self._current_results += current_size
                    self._entries.extend(data)

                    # if exceeds the required number, exit
                    # `self._current_after` is automatically set
                    if num_entry >= number:
                        break

                if self._current_after is None:
                    break

        print(f'{num_entry} entries are qualified and added to the download list.')

        return num_entry

    async def _download_async(self, number: int) -> List[EntryArchive]:
        '''
        Download required entries asynchronously.

        Params:
            number (int): number of **entries** to download

        Returns:
            A list of EntryArchive
        '''
        semaphore = Semaphore(self._semaphore)

        with progressbar(length=number, label=f'Downloading {number} entries...') as bar:
            async with AsyncClient(timeout=Timeout(timeout=300)) as session:
                tasks = [asyncio.create_task(
                    self._acquire(
                        ids, session, semaphore, bar)) for ids in self._entries[:number]]
                results = await asyncio.gather(*tasks)

        return [result for result in results if result]

    async def _acquire(
            self, ids: tuple[str, str],
            session: AsyncClient,
            semaphore: Semaphore,
            bar
    ) -> EntryArchive | None:
        '''
        Perform the download task.

        Params:
            upload (Tuple[str, int]): upload

            session (httpx.AsyncClient): httpx client

            semaphore (asyncio.Semaphore): semaphore

        Returns:
            A list of EntryArchive
        '''

        entry_id, upload_id = ids

        request = self._download_request(entry_id)

        async with semaphore:
            response = await session.post(self._download_url, json=request, headers=self._auth.headers())
            bar.update(1)
            self._entries.remove(ids)

            if response.status_code >= 400:
                print(
                    f'Request with entry id {entry_id} returns {response.status_code},'
                    f' will retry in the next download call...')
                self._entries.append(ids)
                return None

            # successfully downloaded data
            context = ClientContext(self._url, upload_id=upload_id, auth=self._auth)
            result = EntryArchive.m_from_dict(response.json()['data'][0]['archive'], m_context=context)

            if not result:
                print(f'No result returned for id {entry_id}, is the query proper?')

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

        pending_size: int = len(self._entries)

        # download all at once
        if number == 0:
            number = pending_size

        if number == 0:
            # empty list, fetch as many as possible
            number = self.fetch()
        elif pending_size < number:
            # if not sufficient fetched entries, fetch first
            self.fetch(number - pending_size)

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

        pending_size: int = len(self._entries)

        # download all at once
        if number == 0:
            number = pending_size

        if number == 0:
            # empty list, fetch as many as possible
            number = await self.async_fetch()
        elif pending_size < number:
            # if not sufficient fetched entries, fetch first
            await self.async_fetch(number - pending_size)

        return await self._download_async(number)

    def entry_list(self) -> list[tuple[str, str]]:
        return self._entries
