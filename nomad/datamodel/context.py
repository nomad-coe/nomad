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

from typing import Dict
from urllib.parse import urlsplit, urlunsplit
import re
import os.path

import requests

from nomad import utils, config
from nomad.metainfo import Context as MetainfoContext, MSection, Quantity, MetainfoReferenceError
from nomad.datamodel import EntryArchive
from nomad.parsing.parser import ArchiveParser


class Context(MetainfoContext):
    '''
    The nomad implementation of a metainfo context.
    '''

    def __init__(self, installation_url: str = None):
        # take installation_url and ensure it has no trailing slash
        if installation_url is None:
            self.installation_url = config.api_url(api='api/v1')
        elif installation_url.endswith('/'):
            self.installation_url = installation_url[:-1]
        else:
            self.installation_url = installation_url

        # now check if the installation_url points to 'api/v1'
        if not self.installation_url.endswith('/api/v1'):
            self.installation_url += '/v1' if self.installation_url.endswith('/api') else '/api/v1'

        self.archives: Dict[str, MSection] = {}
        self.urls: Dict[MSection, str] = {}

    @property
    def upload_id(self):
        return None

    @staticmethod
    def _get_ids(root: MSection) -> tuple:
        if isinstance(root, EntryArchive) and root.metadata:
            return root.metadata.upload_id, root.metadata.entry_id

        return None, None

    @staticmethod
    def _normalize_fragment(fragment: str) -> str:
        if fragment and fragment != '' and not fragment.startswith('/'):
            fragment = f'/{fragment}'
        return fragment

    def create_reference(self, section: MSection, quantity_def: Quantity, value: MSection) -> str:
        fragment = value.m_path()

        source_root: MSection = section.m_root()
        target_root: MSection = value.m_root()

        if source_root == target_root:
            return f'#{fragment}'

        if target_root in self.urls:
            return f'{self.urls[target_root]}#{fragment}'

        upload_id, entry_id = self._get_ids(target_root)
        assert entry_id is not None, 'Only archives with entry_id can be referenced'
        assert upload_id is not None, 'Only archives with upload_id can be referenced'
        assert source_root.m_context == self, 'Can only create references from archives with this context'
        installation_url = getattr(target_root.m_context, 'installation_url', self.installation_url)

        if installation_url != self.installation_url:
            return f'{installation_url}/uploads/{upload_id}/archive/{entry_id}#{fragment}'

        source_upload_id, source_entry_id = self._get_ids(source_root)

        if entry_id == source_entry_id:
            return f'#{fragment}'

        if upload_id == source_upload_id:
            return f'../upload/archive/{entry_id}#{fragment}'

        return f'../uploads/{upload_id}/archive/{entry_id}#{fragment}'

    def normalize_reference(self, source: MSection, url: str) -> str:
        '''
        Replace mainfile references with entry-based references.
        '''
        url_parts = urlsplit(url)
        fragment = self._normalize_fragment(url_parts.fragment)
        path = url_parts.path
        match = re.search(r'/archive/mainfile/(.*)$', path)
        if not match:
            return urlunsplit(url_parts[0:4] + (fragment,))

        mainfile = match.group(1)
        upload_id = self.upload_id
        if upload_id is None:
            root_section: MSection = source.m_root()
            upload_id = root_section.metadata.upload_id
        assert upload_id is not None, 'Only archives with upload_id can be referenced'
        entry_id = utils.generate_entry_id(upload_id, mainfile)
        path = path.replace(f'/archive/mainfile/{mainfile}', f'/archive/{entry_id}')
        return urlunsplit((url_parts.scheme, url_parts.netloc, path, url_parts.query, fragment,))

    def load_archive(self, entry_id: str, upload_id: str, installation_url: str) -> EntryArchive:
        ''' Loads the archive for the given identification. '''
        raise NotImplementedError()

    def load_raw_file(self, path: str, upload_id: str, installation_url: str) -> MSection:
        ''' Loads a raw file based on the given upload and path. Interpret as metainfo data. '''
        raise NotImplementedError()

    def _parse_url(self, url: str) -> tuple:
        installation_url = self.installation_url
        upload_id = self.upload_id

        url_match = re.match(r'^../uploads?/(\w*)/?(archive|raw)/([^?]+)$', url)
        if url_match:
            if url_match.group(1) is not '':
                upload_id = url_match.group(1)
            kind = url_match.group(2)
            path = url_match.group(3)
            return installation_url, upload_id, kind, path

        url_match = re.search(r'(?<!\.)/uploads/(\w+)/(archive|raw)/([^?]+)$', url)
        if url_match:
            installation_url = url.replace(url_match.group(0), '')
            upload_id = url_match.group(1)
            kind = url_match.group(2)
            path = url_match.group(3)
            return installation_url, upload_id, kind, path

        raise MetainfoReferenceError(f'the url {url} is not a valid metainfo reference')

    def load_url(self, url: str) -> MSection:
        installation_url, upload_id, kind, path = self._parse_url(url)
        if kind == 'archive':
            if path.startswith('mainfile/'):
                entry_id = utils.generate_entry_id(upload_id, path.replace('mainfile/', ''))
            elif '/' in path:
                raise MetainfoReferenceError(f'the url {url} is not a valid metainfo reference')
            else:
                entry_id = path

            return self.load_archive(entry_id, upload_id, installation_url)

        if kind == 'raw':
            return self.load_raw_file(path, upload_id, installation_url)

        raise MetainfoReferenceError(f'the url {url} is not a valid metainfo reference')

    def resolve_archive_url(self, url: str) -> MSection:
        if url not in self.archives:
            archive = self.load_url(url)
            self.archives[url] = archive
            self.urls[archive] = url

        return self.archives[url]


class ServerContext(Context):
    def __init__(self, upload=None):
        super().__init__()
        self.upload = upload

    @property
    def upload_files(self):
        if self.upload:
            return self.upload.upload_files

    @property
    def upload_id(self):
        if self.upload:
            return self.upload.upload_id

    def _get_upload_files(self, upload_id: str, installation_url: str):
        if installation_url and self.installation_url != installation_url:
            raise NotImplementedError()

        if self.upload_files and self.upload_files.upload_id == upload_id:
            return self.upload_files

        # delayed import, context is part of datamodel which should be available in
        # base install, files however requires [infrastructure].
        # TODO move server context to some infrastructure package!
        from nomad import files
        return files.UploadFiles.get(upload_id)

    def load_archive(self, entry_id: str, upload_id: str, installation_url: str) -> EntryArchive:
        upload_files = self._get_upload_files(upload_id, installation_url)

        try:
            archive_dict = upload_files.read_archive(entry_id)[entry_id].to_dict()
        except KeyError:
            raise MetainfoReferenceError(f'archive does not exist {entry_id}')

        return EntryArchive.m_from_dict(archive_dict, m_context=self)

    def load_raw_file(self, path: str, upload_id: str, installation_url: str) -> EntryArchive:
        upload_files = self._get_upload_files(upload_id, installation_url)

        try:
            archive = EntryArchive(m_context=self)
            from nomad.parsing.parser import ArchiveParser
            with upload_files.raw_file(path, 'rt') as f:
                ArchiveParser().parse_file(path, f, archive)
            return archive
        except Exception:
            raise MetainfoReferenceError(f'Could not load {path}.')

    def raw_file(self, *args, **kwargs):
        return self.upload_files.raw_file(*args, **kwargs)


class ClientContext(Context):
    '''
    Since it is a client side context, use config.client.url by default.
    Otherwise, if invoked from ArchiveQuery, use the url provided by user.

    Arguments:
        - installation_url: The installation_url that should be used for intra installation
            references.
        - local_dir: For intra "upload" references, files will be looked up here.
    '''
    def __init__(self, installation_url: str = None, local_dir: str = None):
        super().__init__(config.client.url + '/v1' if installation_url is None else installation_url)
        self.local_dir = local_dir

    def load_archive(self, entry_id: str, upload_id: str, installation_url: str) -> EntryArchive:
        # TODO currently upload_id might be None
        if installation_url is None:
            installation_url = config.api_url(api='api/v1')

        response = requests.get(f'{installation_url}/uploads/{upload_id}/archive/{entry_id}')

        if response.status_code != 200:
            raise MetainfoReferenceError(f'cannot retrieve archive {entry_id} from {installation_url}')

        return EntryArchive.m_from_dict(response.json()['data']['archive'], m_context=self)

    def load_raw_file(self, path: str, upload_id: str, installation_url: str) -> MSection:
        # TODO currently upload_id might be None
        if upload_id is None:
            # try to find a local file, useful when the context is used for local parsing
            file_path = os.path.join(self.local_dir, path)
            if os.path.exists(file_path):
                with open(file_path, 'rt') as f:
                    archive = EntryArchive(m_context=self)
                    ArchiveParser().parse_file(file_path, f, archive)
                    return archive

            raise MetainfoReferenceError(f'cannot retrieve raw file without upload id')

        if installation_url is None:
            installation_url = config.api_url(api='api/v1')

        url = f'{installation_url}/uploads/{upload_id}/raw/{path}'
        response = requests.get(url)

        if response.status_code != 200:
            raise MetainfoReferenceError(f'cannot retrieve raw file {path} from {url}')

        archive_data = response.json()

        if 'm_def' not in archive_data:
            return EntryArchive.m_from_dict(archive_data, m_context=self)
        else:
            return MSection.from_dict(archive_data, m_context=self)
