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

from typing import Tuple, Dict
from urllib.parse import urljoin, urlsplit, urlunsplit
import re
import json

from nomad import utils, config
from nomad.metainfo import Context as MetainfoContext, MSection, Quantity, MetainfoReferenceError
from nomad.datamodel import EntryArchive


class Context(MetainfoContext):
    '''
    The nomad implementation of a metainfo context.
    '''

    def __init__(self, installation_url: str = None):
        self.installation_url = installation_url
        if self.installation_url is None:
            self.installation_url = config.api_url(api='api/v1')

        self.archives: Dict[str, MSection] = {}
        self.urls: Dict[MSection, str] = {}

    @property
    def upload_id(self):
        return None

    def _get_ids(self, root: MSection) -> Tuple[str, str]:
        if not isinstance(root, EntryArchive):
            return None, None

        if root.metadata:
            return root.metadata.upload_id, root.metadata.entry_id

        return None, None

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
            return urljoin(installation_url, f'/api/v1/uploads/{upload_id}/archive/{entry_id}#{fragment}')

        source_upload_id, source_entry_id = self._get_ids(source_root)

        if entry_id == source_entry_id:
            return f'#{fragment}'

        if upload_id == source_upload_id:
            return f'../upload/archive/{entry_id}#{fragment}'

        return f'../uploads/{upload_id}/archive/{entry_id}#{fragment}'

    def _normalize_fragment(self, fragment):
        if fragment and fragment != '' and not fragment.startswith('/'):
            fragment = f'/{fragment}'
        return fragment

    def normalize_reference(self, source: MSection, url: str) -> str:
        ''' Replace mainfile references with entry-based references. '''
        url_parts = urlsplit(url)
        path = url_parts.path
        match = re.search(r'/archive/mainfile/(.*)$', path)
        if not match:
            return urlunsplit(url_parts[0:4] + (self._normalize_fragment(url_parts.fragment),))

        mainfile = match.group(1)
        upload_id = self.upload_id
        if upload_id is None:
            root_section: MSection = source.m_root()
            upload_id = root_section.metadata.upload_id
        assert upload_id is not None, 'Only archives with upload_id can be referenced'
        entry_id = utils.generate_entry_id(upload_id, mainfile)
        path = path.replace(f'/archive/mainfile/{mainfile}', f'/archive/{entry_id}')
        return urlunsplit((
            url_parts.scheme,
            url_parts.netloc,
            path,
            url_parts.query,
            self._normalize_fragment(url_parts.fragment),))

    def load_archive(self, entry_id: str, upload_id: str, installation_url: str) -> EntryArchive:
        ''' Loads the archive for the given identification. '''
        raise NotImplementedError()

    def load_raw_file(self, path: str, upload_id: str, installation_url: str) -> MSection:
        ''' Loads a raw file based on the given upload and path. Interprets is as metainfo data. '''
        raise NotImplementedError()

    def _parse_url(self, url: str) -> Tuple[str, str, str, str]:
        installation_url = self.installation_url
        upload_id = self.upload_id

        url_match = re.match(r'^../upload/(archive|raw)/([^\?]+)$', url)
        if url_match:
            kind = url_match.group(1)
            path = url_match.group(2)
            return installation_url, upload_id, kind, path

        url_match = re.match(r'^../uploads/([A-Za-z0-9_]+)/(archive|raw)/([^\?]+)$', url)
        if url_match:
            upload_id = url_match.group(1)
            kind = url_match.group(2)
            path = url_match.group(3)
            return installation_url, upload_id, kind, path

        url_match = re.search(r'/api/v1/uploads/([A-Za-z0-9_]+)/(archive|raw)/([^\?]+)$', url)
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
    def __init__(self, *args, upload=None, **kwargs):
        super().__init__(*args, **kwargs)
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
        else:
            # delayed import, context is part of datamodel which should be available in
            # base install, files however requires [insfrastructure].
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
            with upload_files.raw_file(path, 'rt') as f:
                archive_data = json.load(f)
        except Exception:
            raise MetainfoReferenceError(f'Could not load {path}.')

        if 'm_def' not in archive_data:
            return EntryArchive.m_from_dict(archive_data, m_context=self)
        else:
            return MSection.from_dict(archive_data, m_context=self)


class ClientContext(Context):
    def load_archive(self, entry_id: str, upload_id: str, installation_url: str) -> EntryArchive:
        raise NotImplementedError()

    def load_raw_file(self, path: str, upload_id: str, installation_url: str) -> MSection:
        raise NotImplementedError()
