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
import pytest
import json

from nomad import utils
from nomad.datamodel import Context
from nomad.datamodel.datamodel import EntryArchive, EntryMetadata


@pytest.fixture(scope='module')
def context():
    class TestContext(Context):
        def load_archive(self, entry_id: str, upload_id: str, installation_url: str) -> EntryArchive:
            assert installation_url is None
            return EntryArchive(metadata=EntryMetadata(entry_id=entry_id, upload_id=upload_id))

    return TestContext()


@pytest.mark.parametrize('url, result', [
    pytest.param('#/root', (None, None, None, '/root'), id='fragment'),
    pytest.param('#root', (None, None, None, '/root'), id='fragment-slash'),
    pytest.param('../upload/archive/entry_id', (None, None, 'entry_id', None), id='only-path'),
    pytest.param('../upload/archive/entry_id#root', (None, None, 'entry_id', '/root'), id='entry'),
    pytest.param('../uploads/upload_id/archive/entry_id#root', (None, 'upload_id', 'entry_id', '/root'), id='upload'),
    pytest.param(
        'https://oasis.de/nomad/api/v1/uploads/upload_id/archive/entry_id#root',
        ('https://oasis.de/nomad', 'upload_id', 'entry_id', '/root'), id='oasis')
])
def test_parse_url(context, url, result):
    context = Context()
    assert context.parse_url(url) == result


@pytest.mark.parametrize('url, result', [
    pytest.param('#/root', '#/root', id='fragment'),
    pytest.param('#root', '#/root', id='fragment-slash'),
    pytest.param('../upload/archive/entry_id#root', '../upload/archive/entry_id#/root', id='entry'),
    pytest.param(
        '../upload/archive/mainfile/path#root',
        f'../upload/archive/{utils.generate_entry_id("test_id", "path")}#/root',
        id='mainfile')
])
def test_normalize_reference(context, url, result):
    root_section = EntryArchive(metadata=EntryMetadata(upload_id='test_id'))
    assert context.normalize_reference(root_section, url) == result


@pytest.mark.parametrize('source, target_archive, target_path, result', [
    pytest.param(
        '''{ "run": [{ "system": [{}] }]}''',
        None,
        '/run/0/system/0', '#/run/0/system/0', id='intra-archive'
    ),
    pytest.param(
        '''{ "metadata": { "upload_id": "source", "entry_id": "source" }}''',
        '''{ "metadata": { "upload_id": "source", "entry_id": "target" }, "run": [{ "system": [{}] }]}''',
        '/run/0/system/0', '../upload/archive/target#/run/0/system/0', id='intra-upload'
    ),
    pytest.param(
        '''{ "metadata": { "upload_id": "source", "entry_id": "source" }}''',
        '''{ "metadata": { "upload_id": "target", "entry_id": "target" }, "run": [{ "system": [{}] }]}''',
        '/run/0/system/0', '../uploads/target/archive/target#/run/0/system/0', id='intra-oasis'
    )
])
def test_create_reference(context, source, target_archive, target_path, result):
    source = EntryArchive.m_from_dict(json.loads(source))
    source.m_context = context

    if target_archive is None:
        target_archive = source
    else:
        target_archive = EntryArchive.m_from_dict(json.loads(target_archive))

    target = target_archive.m_resolve(target_path)

    assert context.create_reference(source, target_archive, target) == result


@pytest.mark.parametrize('url', [
    pytest.param('../upload/archive/entry', id='intra-upload'),
    pytest.param('../uploads/upload/archive/entry', id='intra-oasis'),
])
def test_resolve_archive(context, url):
    target = context.resolve_archive(url)
    assert target is not None
    _, upload_id, entry_id, _ = context.parse_url(url)
    assert target.metadata.entry_id == entry_id
    assert target.metadata.upload_id == upload_id
