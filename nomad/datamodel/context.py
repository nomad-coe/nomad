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
from nomad.datamodel.util import parse_path
from nomad.datamodel.datamodel import EntryMetadata
from nomad.metainfo import Context as MetainfoContext, MSection, Quantity, MetainfoReferenceError
from nomad.datamodel import EntryArchive


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
    def _get_ids(root: MSection, required: bool = True) -> tuple:
        upload_id, entry_id = None, None
        if isinstance(root, EntryArchive) and root.metadata:
            upload_id, entry_id = root.metadata.upload_id, root.metadata.entry_id
        if required:
            assert entry_id is not None, 'Only archives with entry_id can be referenced'
            assert upload_id is not None, 'Only archives with upload_id can be referenced'
        return upload_id, entry_id

    @staticmethod
    def _normalize_fragment(fragment: str) -> str:
        if fragment and fragment != '' and not fragment.startswith('/'):
            fragment = f'/{fragment}'
        return fragment

    def create_reference(
            self, section: MSection, quantity_def: Quantity, value: MSection,
            global_reference: bool = False
    ) -> str:
        fragment = value.m_path()
        target_root: MSection = value.m_root()

        if global_reference:
            upload_id, entry_id = self._get_ids(target_root, required=True)

            return f'../uploads/{upload_id}/archive/{entry_id}#{fragment}'

        source_root: MSection = section.m_root()

        if source_root == target_root:
            return f'#{fragment}'

        if target_root in self.urls:
            return f'{self.urls[target_root]}#{fragment}'

        # when the package is loaded from mongo, there is no metadata
        # and since it is versioned, it is not ideal to create a reference based on processed data
        if getattr(target_root, 'metadata', None) is None:
            return None  # type: ignore

        upload_id, entry_id = self._get_ids(target_root, required=True)
        assert source_root.m_context == self, 'Can only create references from archives with this context'
        installation_url = getattr(target_root.m_context, 'installation_url', self.installation_url)

        if installation_url != self.installation_url:
            return f'{installation_url}/uploads/{upload_id}/archive/{entry_id}#{fragment}'

        source_upload_id, source_entry_id = self._get_ids(source_root, required=False)

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

    def load_raw_file(self, path: str, upload_id: str, installation_url: str, url: str = None) -> MSection:
        ''' Loads a raw file based on the given upload and path. Interpret as metainfo data. '''
        raise NotImplementedError()

    def raw_path_exists(self, path: str) -> bool:
        ''' Use to check if a raw path already exists. '''
        raise NotImplementedError()

    def raw_path(self) -> str:
        ''' The path to the uploads raw files directory. '''
        return os.path.curdir

    def process_updated_raw_file(self, path, allow_modify=False):
        '''
        Use when parsing or normalizing a file causes another file to be added or modified.
        Call this method from the parse/normalize method for the first file, after adding/updating
        the other file is complete. The provided path should denote the added/modified file.
        If we have a ServerContext and the added/modified file matches a parser, we will
        trigger processing of this file (synchronously or asynchronously, depending on if
        we're processing locally or not). For non-ServerContexts, the method does nothing.

        Note, that *modifying* existing files is discouraged, as this needs to be done with
        care to avoid infinite loops of files modifying each other etc. We would thus recommend
        to only use this method for *adding* files. If you still want to modify existing
        files, you must set the `allow_modify` flag, otherwise the call will raise an exception
        if the entry already exists. Also note that this method should not be used to modify
        the same file (i.e. the file that you're currently parsing/normalizing), only when
        adding/modifying other files.
        '''
        pass

    def _parse_url(self, url: str) -> tuple:
        url_results = parse_path(f'{url}#/placeholder', self.upload_id)
        if url_results is None:
            raise MetainfoReferenceError(f'the url {url} is not a valid metainfo reference')

        installation_url, upload_id, entry_id_or_mainfile, kind, _fragment = url_results  # _fragment === /placeholder

        if installation_url is None:
            installation_url = self.installation_url

        return installation_url, upload_id, kind, entry_id_or_mainfile

    def load_url(self, url: str) -> MSection:
        installation_url, upload_id, kind, entry_id_or_mainfile = self._parse_url(url)

        if kind not in ('archive', 'raw'):
            raise MetainfoReferenceError(f'the url {url} is not a valid metainfo reference')

        if kind == 'archive':
            return self.load_archive(entry_id_or_mainfile, upload_id, installation_url)

        return self.load_raw_file(entry_id_or_mainfile, upload_id, installation_url, url)

    def resolve_archive_url(self, url: str) -> MSection:
        if url not in self.archives:
            self.cache_archive(url, self.load_url(url))

        return self.archives[url]

    def cache_archive(self, url: str, archive):
        self.archives[url] = archive
        self.urls[archive] = url


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
        upload_files = files.UploadFiles.get(upload_id)
        assert upload_files and upload_files.upload_id == upload_id
        return upload_files

    def load_archive(self, entry_id: str, upload_id: str, installation_url: str) -> EntryArchive:
        upload_files = self._get_upload_files(upload_id, installation_url)

        try:
            archive_dict = upload_files.read_archive(entry_id)[entry_id].to_dict()
        except KeyError:
            if upload_id != self.upload_id:
                raise MetainfoReferenceError(f'Referencing another Upload is not allowed.')
            from nomad.processing import Entry
            if entry := Entry.objects(entry_id=entry_id).first():
                return self.load_raw_file(entry.mainfile, upload_id, installation_url)
            raise MetainfoReferenceError(f'Could not load {entry_id}.')

        context = self
        if upload_id != self.upload_id:
            from nomad.processing import Upload
            context = ServerContext(Upload(upload_id=upload_id))

        return EntryArchive.m_from_dict(archive_dict, m_context=context)

    def load_raw_file(self, path: str, upload_id: str, installation_url: str, url: str = None) -> EntryArchive:
        upload_files = self._get_upload_files(upload_id, installation_url)

        try:
            # Make sure the archive has proper entry_id, even though we are just
            # loading a raw file. This is important to serialize references into this
            # archive!
            archive = EntryArchive(
                m_context=self,
                metadata=EntryMetadata(
                    upload_id=upload_id,
                    mainfile=path,
                    entry_id=utils.generate_entry_id(upload_id, path)))
            from nomad.parsing.parser import ArchiveParser
            parser = ArchiveParser()
            with upload_files.raw_file(path, 'rt') as f:
                parser.parse_file(path, f, archive)
            if url:
                self.cache_archive(url, archive)
            parser.validate_defintions(archive)
            return archive
        except Exception:
            raise MetainfoReferenceError(f'Could not load {path}.')

    def raw_path(self):
        return self.upload_files._raw_dir.os_path

    def raw_file(self, *args, **kwargs):
        return self.upload_files.raw_file(*args, **kwargs)

    def raw_path_exists(self, path) -> bool:
        return self.upload_files.raw_path_exists(path)

    def process_updated_raw_file(self, path, allow_modify=False):
        self.upload.process_updated_raw_file(path, allow_modify)

    def retrieve_package_by_section_definition_id(self, definition_reference: str, definition_id: str) -> dict:
        if '://' not in definition_reference:
            # not a valid url, may be just a plain python name or reference name
            # use information on the current server
            from nomad.app.v1.routers.metainfo import get_package_by_section_definition_id
            return get_package_by_section_definition_id(definition_id)

        try:
            url_parts = urlsplit(definition_reference)
        except ValueError:
            raise MetainfoReferenceError(f'cannot retrieve section {definition_id} from {definition_reference}')

        # appears to be a valid url
        # build the corresponding request to retrieve the definition
        # important: here we assume the original reference url has the following form:
        #     https://example.nomad.site/some/prefix?possible=query#<definition_contains_id>
        # The definition_id is extracted from the url and used to build the request.
        # The target endpoint is assumed to be
        #     https://example.nomad.site/some/prefix/metainfo/<definition_id>
        response = requests.get(
            urlunsplit((url_parts.scheme, url_parts.netloc, url_parts.path, f'metainfo/{definition_id}', '',)))

        if response.status_code >= 400:
            raise MetainfoReferenceError(f'cannot retrieve section {definition_id} from {definition_reference}')

        return response.json()['data']


def _validate_url(url):
    return config.api_url(api='api/v1') if url is None else url


class ClientContext(Context):
    '''
    Since it is a client side context, use config.client.url by default.
    Otherwise, if invoked from ArchiveQuery, use the url provided by user.

    Arguments:
        - installation_url: The installation_url that should be used for intra installation
            references.
        - local_dir: For intra "upload" references, files will be looked up here.
    '''

    def __init__(
        self, installation_url: str = None, local_dir: str = None, upload_id: str = None,
        username: str = None, password: str = None, auth=None
    ):
        super().__init__(config.client.url + '/v1' if installation_url is None else installation_url)
        self._installation_url = installation_url
        self.local_dir = local_dir
        if auth:
            self._auth = auth
        else:
            from nomad.client import Auth
            self._auth = Auth(user=username, password=password)
        self._upload_id = upload_id

    @property
    def upload_id(self):
        return self._upload_id

    def load_archive(self, entry_id: str, upload_id: str, installation_url: str) -> EntryArchive:
        # TODO currently upload_id might be None
        if upload_id is None:
            url = f'{_validate_url(installation_url)}/entries/{entry_id}/archive'
        else:
            url = f'{_validate_url(installation_url)}/uploads/{upload_id}/archive/{entry_id}'

        response = requests.get(url, auth=self._auth)

        if response.status_code != 200:
            raise MetainfoReferenceError(f'cannot retrieve archive {entry_id} from {installation_url}')

        context = self
        if upload_id != self.upload_id:
            context = ClientContext(installation_url=self._installation_url, local_dir=self.local_dir, upload_id=upload_id, auth=self._auth)

        return EntryArchive.m_from_dict(response.json()['data']['archive'], m_context=context)

    def load_raw_file(self, path: str, upload_id: str, installation_url: str, url: str = None) -> MSection:
        # TODO currently upload_id might be None
        if upload_id is None:
            # try to find a local file, useful when the context is used for local parsing
            file_path = os.path.join(self.local_dir, path) if self.local_dir else path
            if os.path.exists(file_path):
                from nomad.parsing.parser import ArchiveParser
                with open(file_path, 'rt') as f:
                    archive = EntryArchive(m_context=self)
                    ArchiveParser().parse_file(file_path, f, archive)
                    return archive

            raise MetainfoReferenceError(f'cannot retrieve raw file without upload id')

        entry_id = utils.generate_entry_id(upload_id, path)
        return self.load_archive(entry_id, upload_id, installation_url)

    def raw_file(self, path, *args, **kwargs):
        file_path = os.path.join(self.local_dir, path)
        return open(file_path, *args, **kwargs)

    def create_reference(
        self, section: MSection, quantity_def: Quantity, value: MSection,
        global_reference: bool = False
    ) -> str:
        try:
            return super().create_reference(section, quantity_def, value, global_reference)
        except AssertionError:
            return f'<unavailable url>/#{value.m_path()}'

    def retrieve_package_by_section_definition_id(self, definition_reference: str, definition_id: str) -> dict:
        if definition_reference.startswith('http'):
            try:
                url_parts = urlsplit(definition_reference)
                # it appears to be a valid remote url
                # we assume the netloc is the installation_url
                url = urlunsplit((url_parts.scheme, url_parts.netloc, f'api/v1/metainfo/{definition_id}', '', '',))
            except ValueError:
                # falls back to default installation_url
                url = f'{self.installation_url}/metainfo/{definition_id}'
        else:
            # falls back to default installation_url
            url = f'{self.installation_url}/metainfo/{definition_id}'

        response = requests.get(url)

        if response.status_code >= 400:
            raise MetainfoReferenceError(f'cannot retrieve section {definition_id} from {definition_reference}')

        return response.json()['data']
