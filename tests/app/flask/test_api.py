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

from typing import Any, Iterable
import pytest
import time
import json
import zipfile
import io
import inspect
import datetime
import os.path
from urllib.parse import urlencode
import base64
import itertools
from hashlib import md5

from nomad.app.flask.common import rfc3339DateTime
from nomad.app.flask.api.auth import generate_upload_token
from nomad import files, config, utils, infrastructure
from nomad.search.v0 import SearchRequest
from nomad.search.v1 import search
from nomad.metainfo import search_extension
from nomad.files import UploadFiles, PublicUploadFiles
from nomad.processing import Upload, Calc, SUCCESS
from nomad.datamodel import EntryMetadata, User, Dataset

from tests.test_files import example_file, example_file_mainfile, example_file_contents
from tests.test_files import create_staging_upload, create_public_upload, assert_upload_files
from tests.search import assert_search_upload
from tests.processing import test_data as test_processing
from tests.processing.test_data import oasis_publishable_upload  # pylint: disable=unused-import
from tests.conftest import clear_elastic, clear_raw_files
from tests.utils import ExampleData

from ..conftest import create_auth_headers
from .test_app import BlueprintClient

logger = utils.get_logger(__name__)

example_statistics = ['atoms', 'authors', 'dft.system']


@pytest.fixture(scope='function')
def api(client):
    return BlueprintClient(client, '/api')


@pytest.fixture(scope='function')
def test_user_signature_token(api, test_user_auth):
    rv = api.get('/auth/', headers=test_user_auth)
    assert rv.status_code == 200
    return json.loads(rv.data)['signature_token']


def get_upload_entries_metadata(upload: dict) -> Iterable[EntryMetadata]:
    ''' Create a iterable of :class:`EntryMetadata` from a API upload json record. '''
    return [
        EntryMetadata(domain='dft', calc_id=entry['calc_id'], mainfile=entry['mainfile'])
        for entry in upload['calcs']['results']]


def assert_zip_file(rv, files: int = -1, basename: bool = None):
    assert rv.status_code == 200
    assert len(rv.data) > 0
    with zipfile.ZipFile(io.BytesIO(rv.data)) as zip_file:
        assert zip_file.testzip() is None
        zip_files = zip_file.namelist()
        if files >= 0:
            assert len(zip_files) == files
        if basename is not None:
            if basename:
                assert all(
                    os.path.basename(name) == name
                    for name in zip_files if name != 'manifest.json')
            else:
                assert all(
                    os.path.basename(name) != name
                    for name in zip_files for name in zip_files if name != 'manifest.json')


class TestInfo:
    def test_info(self, api, elastic):
        rv = api.get('/info/')
        assert rv.status_code == 200

        data = json.loads(rv.data)
        assert 'codes' in data
        assert 'parsers' in data
        assert 'statistics' in data
        assert len(data['parsers']) >= len(data['codes'])
        assert len(data['domains']) >= 1
        assert rv.status_code == 200

        rv = api.get('/info/')
        assert rv.status_code == 200


class TestKeycloak:
    def test_auth_wo_credentials(self, api, keycloak, no_warn):
        rv = api.get('/auth/')
        assert rv.status_code == 401

    @pytest.fixture(scope='function')
    def auth_headers(self, api, keycloak):
        basic_auth = base64.standard_b64encode(b'sheldon.cooper@nomad-coe.eu:password')
        rv = api.get('/auth/', headers=dict(Authorization='Basic %s' % basic_auth.decode('utf-8')))
        assert rv.status_code == 200
        auth = json.loads(rv.data)
        assert 'access_token' in auth
        assert auth['access_token'] is not None
        return dict(Authorization='Bearer %s' % auth['access_token'])

    def test_auth_with_password(self, api, auth_headers):
        pass

    def test_auth_with_access_token(self, api, auth_headers):
        rv = api.get('/auth/', headers=auth_headers)
        assert rv.status_code == 200

    def assert_sheldon(self, user):
        assert user.email is not None
        assert user.name == 'Sheldon Cooper'
        assert user.first_name == 'Sheldon'
        assert user.last_name == 'Cooper'
        assert user.created is not None
        assert user.affiliation is not None
        assert user.affiliation_address is not None

    def test_get_user(self, keycloak):
        user = infrastructure.keycloak.get_user(username='scooper')
        self.assert_sheldon(user)

    def test_search_user(self, keycloak):
        users = infrastructure.keycloak.search_user(query='Sheldon')
        assert len(users) == 1
        self.assert_sheldon(users[0])


class TestAuth:
    def test_auth_wo_credentials(self, api, no_warn):
        rv = api.get('/auth/')
        assert rv.status_code == 401

    def test_auth_with_token(self, api, test_user_auth):
        rv = api.get('/auth/', headers=test_user_auth)
        assert rv.status_code == 200
        self.assert_auth(api, json.loads(rv.data))

    def assert_auth(self, api, auth):
        assert 'user' not in auth
        assert 'access_token' in auth
        assert 'upload_token' in auth
        assert 'signature_token' in auth

    def test_signature_token(self, test_user_signature_token, no_warn):
        assert test_user_signature_token is not None

    def test_users(self, api):
        rv = api.get('/auth/users?query=Sheldon')
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert len(data['users'])
        user = data['users'][0]
        keys = user.keys()
        required_keys = ['name', 'email', 'user_id']
        assert all(key in keys for key in required_keys)
        for key in keys:
            assert user.get(key) is not None

    def test_invite(self, api, test_user_auth, no_warn):
        rv = api.put(
            '/auth/users', headers=test_user_auth, content_type='application/json',
            data=json.dumps({
                'first_name': 'John',
                'last_name': 'Doe',
                'affiliation': 'Affiliation',
                'email': 'john.doe@affiliation.edu'
            }))
        assert rv.status_code == 200
        data = json.loads(rv.data)
        keys = data.keys()
        required_keys = ['name', 'email', 'user_id']
        assert all(key in keys for key in required_keys)


class TestUploads:

    def assert_uploads(self, upload_json_str, count=0, **kwargs):
        data = json.loads(upload_json_str)
        assert 'pagination' in data
        assert 'page' in data['pagination']

        data = data['results']
        assert isinstance(data, list)
        assert len(data) == count

        if count > 0:
            self.assert_upload(json.dumps(data[0]), **kwargs)

    def assert_upload(self, upload_json_str, id=None, **kwargs):
        data = json.loads(upload_json_str)
        assert 'upload_id' in data
        if id is not None:
            assert id == data['upload_id']
        assert 'create_time' in data

        for key, value in kwargs.items():
            assert data.get(key, None) == value

        return data

    def assert_processing(self, api, test_user_auth, upload_id):
        upload_endpoint = '/uploads/%s' % upload_id

        # poll until completed
        upload = self.block_until_completed(api, upload_id, test_user_auth)

        assert len(upload['tasks']) == 4
        assert upload['tasks_status'] == SUCCESS
        assert upload['current_task'] == 'cleanup'
        assert not upload['process_running']

        calcs = upload['calcs']['results']
        for calc in calcs:
            assert calc['tasks_status'] == SUCCESS
            assert calc['current_task'] == 'archiving'
            assert len(calc['tasks']) == 3

            assert 'atoms' in calc['metadata']
            assert api.get('/archive/logs/%s/%s' % (calc['upload_id'], calc['calc_id']), headers=test_user_auth).status_code == 200

        if upload['calcs']['pagination']['total'] > 1:
            rv = api.get('%s?page=2&per_page=1&order_by=tasks_status' % upload_endpoint, headers=test_user_auth)
            assert rv.status_code == 200
            upload = self.assert_upload(rv.data)
            assert len(upload['calcs']['results']) == 1

        entries = get_upload_entries_metadata(upload)
        assert_upload_files(upload_id, entries, files.StagingUploadFiles)
        assert_search_upload(entries, index='v0', additional_keys=['atoms', 'dft.system'])

    def assert_published(self, api, test_user_auth, upload_id, proc_infra, metadata={}):
        rv = api.get('/uploads/%s' % upload_id, headers=test_user_auth)
        upload = self.assert_upload(rv.data)

        rv = api.post(
            '/uploads/%s' % upload_id,
            headers=test_user_auth,
            data=json.dumps(dict(operation='publish', metadata=metadata)),
            content_type='application/json')
        assert rv.status_code == 200
        upload = self.assert_upload(rv.data)
        assert upload['current_process'] == 'publish_upload'
        assert upload['process_running']

        additional_keys = ['with_embargo']
        if 'external_id' in metadata:
            additional_keys.append('external_id')

        self.block_until_completed(api, upload_id, test_user_auth)

        upload_proc = Upload.objects(upload_id=upload_id).first()
        assert upload_proc is not None
        assert upload_proc.published is True
        assert upload_proc.embargo_length == min(36, metadata.get('embargo_length', 36))

        with upload_proc.entries_metadata() as entries:
            for entry in entries:
                for key, transform in {
                        'comment': lambda e: e.comment,
                        'with_embargo': lambda e: e.with_embargo,
                        'references': lambda e: e.references,
                        'coauthors': lambda e: [u.user_id for u in e.coauthors],
                        '_uploader': lambda e: e.uploader.user_id,
                        '_pid': lambda e: str(e.pid),
                        'external_id': lambda e: e.external_id}.items():
                    if key in metadata:
                        assert transform(entry) == metadata[key], key

        assert_upload_files(upload_id, entries, files.PublicUploadFiles, published=True)
        assert_search_upload(entries, additional_keys=additional_keys, published=True)

    def block_until_completed(self, api, upload_id: str, test_user_auth):
        while True:
            time.sleep(0.1)
            rv = api.get('/uploads/%s' % upload_id, headers=test_user_auth)
            if rv.status_code == 200:
                upload = self.assert_upload(rv.data)
                if not upload['process_running'] and not upload['tasks_running']:
                    return upload
            elif rv.status_code == 404:
                return None
            else:
                raise Exception(
                    'unexpected status code while blocking for upload processing: %s' %
                    str(rv.status_code))

    def assert_upload_does_not_exist(self, api, upload_id: str, test_user_auth):
        self.block_until_completed(api, upload_id, test_user_auth)

        rv = api.get('/uploads/%s' % upload_id, headers=test_user_auth)
        assert rv.status_code == 404

        assert Upload.objects(upload_id=upload_id).first() is None
        assert Calc.objects(upload_id=upload_id).count() is 0

        mongo_db = infrastructure.mongo_client[config.mongo.db_name]
        mongo_collection = mongo_db['archive']
        assert mongo_collection.find({}).count() == 0

        upload_files = UploadFiles.get(upload_id)
        assert upload_files is None or isinstance(upload_files, PublicUploadFiles)

    def test_get_command(self, api, test_user_auth, no_warn):
        rv = api.get('/uploads/command', headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert 'upload_command' in data
        assert '/api/uploads' in data['upload_command']
        assert 'upload_url' in data

    def test_get_empty(self, api, test_user_auth, no_warn):
        rv = api.get('/uploads/', headers=test_user_auth)

        assert rv.status_code == 200
        self.assert_uploads(rv.data, count=0)

    def test_get_not_existing(self, api, test_user_auth, no_warn):
        rv = api.get('/uploads/123456789012123456789012', headers=test_user_auth)
        assert rv.status_code == 404

    def test_put_upload_token(self, api, non_empty_example_upload, test_user):
        url = '/uploads/?token=%s&local_path=%s&name=test_upload' % (
            generate_upload_token(test_user), non_empty_example_upload)
        rv = api.put(url)
        assert rv.status_code == 200, rv.data
        assert 'Thanks for uploading' in rv.data.decode('utf-8')

    @pytest.mark.parametrize('mode', ['multipart', 'stream', 'local_path'])
    @pytest.mark.parametrize('name', [None, 'test_name'])
    def test_put(self, api, test_user_auth, proc_infra, example_upload, mode, name):
        file = example_upload
        if name:
            url = '/uploads/?name=%s' % name
        else:
            url = '/uploads/'

        if mode == 'multipart':
            rv = api.put(
                url, data=dict(file=(open(file, 'rb'), 'the_name')), headers=test_user_auth)
            if not name:
                name = 'the_name'
        elif mode == 'stream':
            with open(file, 'rb') as f:
                rv = api.put(url, data=f.read(), headers=test_user_auth)
        elif mode == 'local_path':
            url += '&' if name else '?'
            url += 'local_path=%s' % file
            rv = api.put(url, headers=test_user_auth)
        else:
            assert False

        assert rv.status_code == 200
        if mode == 'local_path':
            upload = self.assert_upload(rv.data, upload_path=file, name=name)
        else:
            upload = self.assert_upload(rv.data, name=name)
        assert upload['tasks_running']

        self.assert_processing(api, test_user_auth, upload['upload_id'])

    @pytest.mark.timeout(config.tests.default_timeout)
    def test_upload_limit(self, api, mongo, test_user, test_user_auth, proc_infra):
        for _ in range(0, config.services.upload_limit):
            Upload.create(user=test_user)
        file = example_file
        rv = api.put('/uploads/?local_path=%s' % file, headers=test_user_auth)
        assert rv.status_code == 400
        assert Upload.user_uploads(test_user).count() == config.services.upload_limit

    def test_delete_not_existing(self, api, test_user_auth, no_warn):
        rv = api.delete('/uploads/123456789012123456789012', headers=test_user_auth)
        assert rv.status_code == 404

    @pytest.fixture(scope='function')
    def slow_processing(self, monkeypatch):
        old_cleanup = Upload.cleanup

        def slow_cleanup(self):
            time.sleep(0.5)
            old_cleanup(self)

        monkeypatch.setattr('nomad.processing.data.Upload.cleanup', slow_cleanup)
        yield True
        monkeypatch.setattr('nomad.processing.data.Upload.cleanup', old_cleanup)

    def test_delete_published(self, api, test_user_auth, proc_infra, no_warn):
        rv = api.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(api, test_user_auth, upload['upload_id'])
        self.assert_published(api, test_user_auth, upload['upload_id'], proc_infra)
        rv = api.delete('/uploads/%s' % upload['upload_id'], headers=test_user_auth)
        assert rv.status_code == 400

    def test_delete(self, api, test_user_auth, proc_infra, no_warn):
        rv = api.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(api, test_user_auth, upload['upload_id'])
        rv = api.delete('/uploads/%s' % upload['upload_id'], headers=test_user_auth)
        assert rv.status_code == 200
        self.assert_upload_does_not_exist(api, upload['upload_id'], test_user_auth)

    def test_post_empty(self, api, test_user_auth, empty_upload, proc_infra, no_warn):
        rv = api.put('/uploads/?local_path=%s' % empty_upload, headers=test_user_auth)
        assert rv.status_code == 200
        upload = self.assert_upload(rv.data)
        self.assert_processing(api, test_user_auth, upload['upload_id'])
        rv = api.post(
            '/uploads/%s' % upload['upload_id'], headers=test_user_auth,
            data=json.dumps(dict(operation='publish')),
            content_type='application/json')
        assert rv.status_code == 400

    def test_post(self, api, test_user_auth, non_empty_example_upload, proc_infra, no_warn):
        rv = api.put('/uploads/?local_path=%s' % non_empty_example_upload, headers=test_user_auth)
        assert rv.status_code == 200
        upload = self.assert_upload(rv.data)
        self.assert_processing(api, test_user_auth, upload['upload_id'])
        self.assert_published(api, test_user_auth, upload['upload_id'], proc_infra)

        # still visible
        assert api.get('/uploads/%s' % upload['upload_id'], headers=test_user_auth).status_code == 200
        # still listed with all=True
        rv = api.get('/uploads/?state=all', headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)['results']
        assert len(data) > 0
        assert any(item['upload_id'] == upload['upload_id'] for item in data)
        # not listed with all=False
        rv = api.get('/uploads/', headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)['results']
        assert not any(item['upload_id'] == upload['upload_id'] for item in data)

    def test_post_metadata(
            self, api, proc_infra, admin_user_auth, test_user_auth, test_user,
            other_test_user, no_warn, example_user_metadata):
        rv = api.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(api, test_user_auth, upload['upload_id'])
        metadata = dict(**example_user_metadata)
        metadata['_upload_time'] = datetime.datetime.utcnow().isoformat()
        self.assert_published(api, admin_user_auth, upload['upload_id'], proc_infra, metadata)

    def test_post_metadata_forbidden(self, api, proc_infra, test_user_auth, no_warn):
        rv = api.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(api, test_user_auth, upload['upload_id'])
        rv = api.post(
            '/uploads/%s' % upload['upload_id'],
            headers=test_user_auth,
            data=json.dumps(dict(operation='publish', metadata=dict(_pid='256'))),
            content_type='application/json')
        assert rv.status_code == 401

    def test_post_metadata_and_republish(
            self, api, proc_infra, admin_user_auth, test_user_auth, test_user,
            other_test_user, no_warn, example_user_metadata):
        rv = api.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(api, test_user_auth, upload['upload_id'])
        metadata = dict(**example_user_metadata)
        metadata['_upload_time'] = datetime.datetime.utcnow().isoformat()
        self.assert_published(api, admin_user_auth, upload['upload_id'], proc_infra, metadata)
        self.assert_published(api, admin_user_auth, upload['upload_id'], proc_infra, {})

    def test_post_re_process(self, api, published, test_user_auth, monkeypatch):
        monkeypatch.setattr('nomad.config.meta.version', 're_process_test_version')
        monkeypatch.setattr('nomad.config.meta.commit', 're_process_test_commit')

        upload_id = published.upload_id
        rv = api.post(
            '/uploads/%s' % upload_id,
            headers=test_user_auth,
            data=json.dumps(dict(operation='re-process')),
            content_type='application/json')

        assert rv.status_code == 200, rv.data
        assert self.block_until_completed(api, upload_id, test_user_auth) is not None

    # TODO validate metadata (or all input models in API for that matter)
    # def test_post_bad_metadata(self, api, proc_infra, test_user_auth):
    #     rv = api.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
    #     upload = self.assert_upload(rv.data)
    #     self.assert_processing(api, test_user_auth, upload['upload_id'])
    #     rv = api.post(
    #         '/uploads/%s' % upload['upload_id'],
    #         headers=test_user_auth,
    #         data=json.dumps(dict(operation='publish', metadata=dict(doesnotexist='hi'))),
    #         content_type='application/json')
    #     assert rv.status_code == 400

    @pytest.mark.parametrize('upload_file, ending', [
        ('examples_potcar.zip', ''),
        ('examples_potcar_gz.tgz', '.gz'),
        ('examples_potcar_xz.tgz', '.xz')])
    def test_potcar(self, api, proc_infra, test_user_auth, upload_file, ending):
        # only the owner, shared with people are supposed to download the original potcar file
        example_file = 'tests/data/proc/%s' % upload_file
        rv = api.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)

        upload = self.assert_upload(rv.data)
        upload_id = upload['upload_id']
        self.assert_processing(api, test_user_auth, upload_id)
        self.assert_published(api, test_user_auth, upload_id, proc_infra)
        rv = api.get('/raw/%s/examples_potcar/POTCAR%s' % (upload_id, ending))
        assert rv.status_code == 401
        rv = api.get('/raw/%s/examples_potcar/POTCAR%s' % (upload_id, ending), headers=test_user_auth)
        assert rv.status_code == 200
        rv = api.get('/raw/%s/examples_potcar/POTCAR%s.stripped' % (upload_id, ending))
        assert rv.status_code == 200

    def test_put_publish_directly(self, api, test_user_auth, non_empty_example_upload, proc_infra, no_warn):
        rv = api.put('/uploads/?%s' % urlencode(dict(
            local_path=non_empty_example_upload,
            publish_directly=True)), headers=test_user_auth)
        assert rv.status_code == 200, rv.data
        upload = self.assert_upload(rv.data)
        upload_id = upload['upload_id']

        # poll until completed
        upload = self.block_until_completed(api, upload_id, test_user_auth)

        assert len(upload['tasks']) == 4
        assert upload['tasks_status'] == SUCCESS
        assert upload['current_task'] == 'cleanup'
        assert not upload['process_running']

        upload_proc = Upload.objects(upload_id=upload_id).first()
        assert upload_proc.published

        entries = get_upload_entries_metadata(upload)
        assert_upload_files(upload_id, entries, files.PublicUploadFiles)
        assert_search_upload(entries, index='v0', additional_keys=['atoms', 'dft.system'])

    def test_post_from_oasis_admin(self, api, non_empty_uploaded, other_test_user_auth, test_user, no_warn):
        url = '/uploads/?%s' % urlencode(dict(
            local_path=non_empty_uploaded[1], oasis_upload_id='oasis_upload_id',
            oasis_uploader_id=test_user.user_id, oasis_deployment_id='an_id'))
        assert api.put(url, headers=other_test_user_auth).status_code == 401

    def test_post_from_oasis_duplicate(self, api, non_empty_uploaded, test_user, test_user_auth, no_warn):
        Upload.create(upload_id='oasis_upload_id', user=test_user).save()
        url = '/uploads/?%s' % urlencode(dict(
            local_path=non_empty_uploaded[1], oasis_upload_id='oasis_upload_id',
            oasis_uploader_id=test_user.user_id, oasis_deployment_id='an_id'))
        assert api.put(url, headers=test_user_auth).status_code == 400

    def test_post_from_oasis_missing_parameters(self, api, non_empty_uploaded, test_user_auth, test_user, no_warn):
        url = '/uploads/?%s' % urlencode(dict(
            local_path=non_empty_uploaded[1], oasis_upload_id='oasis_upload_id',
            oasis_uploader_id=test_user.user_id))
        assert api.put(url, headers=test_user_auth).status_code == 400

        url = '/uploads/?%s' % urlencode(dict(
            local_path=non_empty_uploaded[1], oasis_upload_id='oasis_upload_id',
            oasis_deployment_id='an_id'))
        assert api.put(url, headers=test_user_auth).status_code == 400

        url = '/uploads/?%s' % urlencode(dict(
            local_path=non_empty_uploaded[1],
            oasis_uploader_id=test_user.user_id, oasis_deployment_id='an_id'))
        assert api.put(url, headers=test_user_auth).status_code == 400

    def test_post_from_oasis(self, api, test_user_auth, test_user, oasis_example_upload, proc_infra, no_warn):
        rv = api.put('/uploads/?%s' % urlencode(dict(
            local_path=oasis_example_upload,
            oasis_upload_id='oasis_upload_id',
            oasis_deployment_id='an_id',
            oasis_uploader_id=test_user.user_id)), headers=test_user_auth)
        assert rv.status_code == 200, rv.data
        upload = self.assert_upload(rv.data)
        upload_id = upload['upload_id']
        assert upload_id == 'oasis_upload_id'

        # poll until completed
        upload = self.block_until_completed(api, upload_id, test_user_auth)

        assert len(upload['tasks']) == 4
        assert upload['tasks_status'] == SUCCESS
        assert upload['current_task'] == 'cleanup'
        assert not upload['process_running']

        upload_proc = Upload.objects(upload_id=upload_id).first()
        assert upload_proc.published
        assert upload_proc.from_oasis

        entries = get_upload_entries_metadata(upload)
        assert_upload_files(upload_id, entries, files.PublicUploadFiles)
        assert_search_upload(entries, index='v0', additional_keys=['atoms', 'dft.system'])

    def test_post_publish_from_oasis(
            self, api, oasis_publishable_upload, other_test_user_auth, monkeypatch, no_warn):
        monkeypatch.setattr('nomad.config.keycloak.oasis', True)
        cn_upload_id, upload = oasis_publishable_upload
        upload_id = upload.upload_id

        rv = api.post(
            '/uploads/%s' % upload_id,
            headers=other_test_user_auth,
            data=json.dumps(dict(operation='publish-to-central-nomad')),
            content_type='application/json')

        assert rv.status_code == 200, rv.data
        upload = self.assert_upload(rv.data)
        assert upload['current_process'] == 'publish_from_oasis'
        assert upload['process_running']

        self.block_until_completed(api, upload_id, other_test_user_auth)

        cn_upload = Upload.objects(upload_id=cn_upload_id).first()
        cn_upload.block_until_complete()
        assert len(cn_upload.errors) == 0
        assert cn_upload.current_process == 'process_upload'


today = datetime.datetime.utcnow().date()
today_datetime = datetime.datetime(*today.timetuple()[:6])


class UploadFilesBasedTests:

    @staticmethod
    def fix_signature(func, wrapper):
        additional_args = list(inspect.signature(func).parameters.values())[4:]
        wrapper_sig = inspect.signature(wrapper)
        wrapper_args = list(wrapper_sig.parameters.values())[:3] + additional_args
        wrapper_sig = wrapper_sig.replace(parameters=tuple(wrapper_args))
        wrapper.__signature__ = wrapper_sig

    @staticmethod
    def check_authorization(func):
        @pytest.mark.parametrize('test_data', [
            [True, None, True],      # in staging for upload
            [True, None, False],     # in staging for different user
            [True, None, None],      # in staging for guest
            [True, None, 'admin'],   # in staging, for admin
            [False, True, True],     # in public, restricted for uploader
            [False, True, False],    # in public, restricted for different user
            [False, True, None],     # in public, restricted for guest
            [False, True, 'admin'],  # in public, restricted for admin
            [False, False, True],    # in public, public, for uploader
            [False, False, False],   # in public, public, for different user
            [False, False, None]     # in public, public, for guest
        ], indirect=True)
        def wrapper(self, api, test_data, *args, **kwargs):
            upload, authorized, auth_headers = test_data
            try:
                func(self, api, upload, auth_headers, *args, **kwargs)
            except AssertionError as assertion:
                assertion_str = str(assertion)
                if not authorized:
                    if '0 == 5' in assertion_str:
                        # the user is not authorized an gets an empty zip as expected
                        return
                    if '401' in assertion_str:
                        # the user is not authorized and gets a 401 as expected
                        return
                raise assertion

            if not authorized:
                assert False
        UploadFilesBasedTests.fix_signature(func, wrapper)
        return wrapper

    @staticmethod
    def ignore_authorization(func):
        @pytest.mark.parametrize('test_data', [
            [True, None, True],      # in staging
            [False, False, None],    # in public
        ], indirect=True)
        def wrapper(self, api, test_data, *args, **kwargs):
            upload, _, auth_headers = test_data
            func(self, api, upload, auth_headers, *args, **kwargs)
        UploadFilesBasedTests.fix_signature(func, wrapper)
        return wrapper

    @pytest.fixture(scope='function')
    def test_data(self, request, mongo, raw_files, no_warn, test_user, other_test_user, admin_user):
        # delete potential old test files
        for _ in [0, 1]:
            upload_files = UploadFiles.get('test_upload')
            if upload_files:
                upload_files.delete()

        in_staging, restricted, for_uploader = request.param

        if in_staging:
            authorized = for_uploader is True or for_uploader == 'admin'
        else:
            authorized = not restricted or for_uploader is True or for_uploader == 'admin'

        if for_uploader is True:
            auth_headers = create_auth_headers(test_user)
        elif for_uploader is False:
            auth_headers = create_auth_headers(other_test_user)
        elif for_uploader == 'admin':
            auth_headers = create_auth_headers(admin_user)
        else:
            auth_headers = None

        calc_specs = 'r' if restricted else 'p'
        Upload.create(user=test_user, upload_id='test_upload')
        if in_staging:
            _, _, upload_files = create_staging_upload('test_upload', calc_specs=calc_specs)
        else:
            _, _, upload_files = create_public_upload('test_upload', calc_specs=calc_specs)

        yield 'test_upload', authorized, auth_headers

        upload_files.delete()


class TestArchive(UploadFilesBasedTests):
    @UploadFilesBasedTests.check_authorization
    def test_get(self, api, upload, auth_headers):
        rv = api.get('/archive/%s/0' % upload, headers=auth_headers)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert data is not None
        assert 'section_metadata' in data
        assert 'section_run' in data

    @UploadFilesBasedTests.ignore_authorization
    def test_get_signed(self, api, upload, _, test_user_signature_token):
        rv = api.get('/archive/%s/0?signature_token=%s' % (upload, test_user_signature_token))
        assert rv.status_code == 200
        assert json.loads(rv.data) is not None

    @UploadFilesBasedTests.check_authorization
    def test_get_calc_proc_log(self, api, upload, auth_headers):
        rv = api.get('/archive/logs/%s/0' % upload, headers=auth_headers)
        assert rv.status_code == 200
        assert len(rv.data) > 0

    @UploadFilesBasedTests.ignore_authorization
    def test_get_calc_proc_log_signed(self, api, upload, _, test_user_signature_token):
        rv = api.get('/archive/logs/%s/0?signature_token=%s' % (upload, test_user_signature_token))
        assert rv.status_code == 200
        assert len(rv.data) > 0

    @UploadFilesBasedTests.ignore_authorization
    def test_get_non_existing_archive(self, api, upload, auth_headers):
        rv = api.get('/archive/%s' % 'doesnt/exist', headers=auth_headers)
        assert rv.status_code == 404

    @pytest.mark.parametrize('compress', [False, True])
    def test_archive_zip_dowload_upload_id(self, api, non_empty_processed, test_user_auth, compress):
        url = '/archive/download?upload_id=%s&compress=%s' % (non_empty_processed.upload_id, 'true' if compress else 'false')
        rv = api.get(url, headers=test_user_auth)

        assert rv.status_code == 200
        assert_zip_file(rv, files=2)

    @pytest.mark.parametrize('query_params', [
        {'atoms': 'Si'},
        {'authors': 'Sheldon Cooper'}
    ])
    def test_archive_zip_dowload(self, api, processeds, test_user_auth, query_params):

        url = '/archive/download?%s' % urlencode(query_params)
        rv = api.get(url, headers=test_user_auth)

        assert rv.status_code == 200
        assert_zip_file(rv, files=len(processeds) + 1)
        with zipfile.ZipFile(io.BytesIO(rv.data)) as zip_file:
            with zip_file.open('manifest.json', 'r') as f:
                manifest = json.load(f)
                assert len(manifest) == len(processeds)

    def test_archive_zip_dowload_empty(self, api, elastic):
        url = '/archive/download?upload_id=doesNotExist'
        rv = api.get(url)

        assert rv.status_code == 200
        assert_zip_file(rv, files=1)

    @pytest.mark.parametrize('required', [{'section_workflow': '*'}, {'section_run': '*'}])
    def test_archive_query_paginated(self, api, published_wo_user_metadata, required):
        data = {'required': required, 'pagination': {'per_page': 5}}
        uri = '/archive/query'
        rv = api.post(uri, content_type='application/json', data=json.dumps(data))

        assert rv.status_code == 200
        data = rv.get_json()

        assert data
        results = data.get('results', None)
        assert len(results) > 0
        for result in results:
            assert 'calc_id' in result
            assert 'parser_name' in result
            assert 'archive' in result

        # TODO assert archive contents

    def test_archive_not_exists(self, api, published_wo_user_metadata):
        entry_metadata = EntryMetadata(
            domain='dft', upload_id=published_wo_user_metadata.upload_id,
            calc_id='test_id', published=True, with_embargo=False)
        entry_metadata.a_elastic.index(refresh=True)

        data = {}
        uri = '/archive/query'
        rv = api.post(uri, content_type='application/json', data=json.dumps(data))

        rv = api.post(uri, content_type='application/json', data=json.dumps(dict(per_page=5, raise_errors=True)))
        assert rv.status_code == 404

        rv = api.post(uri, content_type='application/json', data=json.dumps(dict(per_page=5, raise_errors=False)))
        assert rv.status_code == 200

    def test_archive_query_aggregated(self, api, published_wo_user_metadata):
        uri = '/archive/query'
        schema = {
            'section_run': {
                'section_single_configuration_calculation': {
                    'energy_total': '*'}}}

        query = {'required': schema, 'aggregation': {'per_page': 1}}

        count = 0
        while True:
            rv = api.post(uri, content_type='application/json', data=json.dumps(query))
            assert rv.status_code == 200
            data = rv.get_json()
            results = data.get('results', None)
            count += len(results)
            after = data['aggregation']['after']
            if after is None:
                break

            query['aggregation']['after'] = after

        assert count > 0

    @pytest.mark.timeout(config.tests.default_timeout)
    @pytest.fixture(scope='function')
    def example_upload(self, proc_infra, test_user):
        path = 'tests/data/proc/example_vasp_with_binary.zip'
        results = []
        for uid in range(2):
            upload_id = 'vasp_%d' % uid
            processed = test_processing.run_processing((upload_id, path), test_user)

            processed.publish_upload()
            try:
                processed.block_until_complete(interval=.01)
            except Exception:
                pass

            results.append(processed)

        return results

    @pytest.mark.parametrize('query_expression, nresults', [
        pytest.param({}, 4, id='empty'),
        pytest.param({'dft.system': 'bulk'}, 4, id='match'),
        pytest.param({'$gte': {'n_atoms': 1}}, 4, id='comparison'),
        pytest.param({
            '$and': [
                {'dft.system': 'bulk'}, {'$not': [{'dft.compound_type': 'ternary'}]}
            ]
        }, 2, id="and-with-not"),
        pytest.param({
            '$or': [
                {'upload_id': ['vasp_0']}, {'$gte': {'n_atoms': 1}}
            ]
        }, 4, id="or-with-gte"),
        pytest.param({
            '$not': [{'dft.spacegroup': 221}, {'dft.spacegroup': 227}]
        }, 0, id="not"),
        pytest.param({
            '$and': [
                {'dft.code_name': 'VASP'},
                {'$gte': {'n_atoms': 3}},
                {'$lte': {'dft.workflow.section_geometry_optimization.final_energy_difference': 1e-24}}
            ]}, 0, id='client-example')
    ])
    def test_post_archive_query(self, api, example_upload, query_expression, nresults):
        data = {'pagination': {'per_page': 5}, 'query': query_expression}
        rv = api.post('/archive/query', content_type='application/json', data=json.dumps(data))

        assert rv.status_code == 200
        data = rv.get_json()

        assert data
        results = data.get('results', None)
        assert len(results) == nresults

    @pytest.mark.parametrize('query', [
        pytest.param({'$bad_op': {'n_atoms': 1}}, id='bad-op')
    ])
    def test_post_archive_bad_query(self, api, query):
        rv = api.post(
            '/archive/query', content_type='application/json',
            data=json.dumps(dict(query=query)))

        assert rv.status_code == 400


class TestMetainfo():
    @pytest.mark.parametrize('package', ['common', 'vasp'])
    def test_regular(self, api, package):
        rv = api.get('/metainfo/%s' % package)
        assert rv.status_code == 200
        assert len(rv.get_json()) > 0

    def test_full_name(self, api):
        rv = api.get('/metainfo/nomad.datamodel.metainfo.common')
        assert rv.status_code == 200

    def test_extension(self, api):
        rv = api.get('/metainfo/common.json')
        assert rv.status_code == 200

        rv = api.get('/metainfo/legacy/common.json')
        assert rv.status_code == 200

        rv = api.get('/metainfo/legacy/common.nomadmetainfo.json')
        assert rv.status_code == 200

    @pytest.mark.parametrize('package', ['common', 'vasp'])
    def test_legacy(self, api, package):
        rv = api.get('/metainfo/legacy/%s' % package)
        assert rv.status_code == 200
        assert len(rv.get_json().get('metaInfos')) > 0

    def test_does_not_exist(self, api):
        rv = api.get('/metainfo/doesnotexist')
        assert rv.status_code == 404
        rv = api.get('/metainfo/legacy/doesnotexist')
        assert rv.status_code == 404

    def test_all(self, api):
        rv = api.get('/metainfo/')
        rv.status_code == 200
        assert len(rv.get_json()) > 0


class TestRepo():
    @pytest.fixture(scope='class')
    def example_elastic_calcs(
            self, elastic_infra, raw_files_infra, normalized,
            test_user: User, other_test_user: User):
        clear_elastic(elastic_infra)

        example_dataset = Dataset(
            dataset_id='ds_id', name='ds_name', user_id=test_user.user_id, doi='ds_doi')
        example_dataset.a_mongo.create()

        entry_metadata = EntryMetadata(
            domain='dft', upload_id='example_upload_id', calc_id='0', upload_time=today_datetime)
        entry_metadata.files = ['test/mainfile.txt']
        entry_metadata.apply_domain_metadata(normalized)
        entry_metadata.encyclopedia = normalized.section_metadata.encyclopedia

        entry_metadata.m_update(datasets=[example_dataset.dataset_id])

        entry_metadata.m_update(
            calc_id='1', uploader=test_user.user_id, published=True, with_embargo=False)
        entry_metadata.a_elastic.index(refresh=True)

        entry_metadata.m_update(
            calc_id='2', uploader=other_test_user.user_id, published=True,
            with_embargo=False, pid='2', upload_time=today_datetime - datetime.timedelta(days=5),
            external_id='external_2')
        entry_metadata.m_update(
            atoms=['Fe'], comment='this is a specific word', formula='AAA')
        entry_metadata.dft.basis_set = 'zzz'
        entry_metadata.a_elastic.index(refresh=True)

        entry_metadata.m_update(
            calc_id='3', uploader=other_test_user.user_id, published=False,
            with_embargo=False, pid='3', external_id='external_3')
        entry_metadata.a_elastic.index(refresh=True)

        entry_metadata.m_update(
            calc_id='4', uploader=other_test_user.user_id, published=True,
            with_embargo=True, pid='4', external_id='external_4')
        entry_metadata.a_elastic.index(refresh=True)

        yield

        example_dataset.a_mongo.delete()
        clear_raw_files()

    def assert_search(self, rv: Any, number_of_calcs: int) -> dict:
        if rv.status_code != 200:
            print(rv.data)
        assert rv.status_code == 200

        data = json.loads(rv.data)

        results = data.get('results', None)
        assert results is not None
        assert isinstance(results, list)
        assert len(results) == number_of_calcs

        return data

    def test_own_calc(self, api, example_elastic_calcs, no_warn, test_user_auth):
        rv = api.get('/repo/0/1', headers=test_user_auth)
        assert rv.status_code == 200

    def test_get_code(self, api, example_elastic_calcs, no_warn, test_user_auth):
        rv = api.get('/repo/0/1', headers=test_user_auth)
        assert rv.status_code == 200
        data = rv.json
        assert data['code']['curl'] is not None
        assert data['code']['python'] is not None
        assert data['code']['clientlib'] is not None

    def test_get_defaults(self, api, example_elastic_calcs, no_warn, test_user_auth):
        rv = api.get('/repo/0/1', headers=test_user_auth)
        assert rv.status_code == 200
        data = rv.json
        assert data['license'] is not None

    def test_public_calc(self, api, example_elastic_calcs, no_warn, other_test_user_auth):
        rv = api.get('/repo/0/1', headers=other_test_user_auth)
        assert rv.status_code == 200

    def test_embargo_calc(self, api, example_elastic_calcs, no_warn, test_user_auth):
        rv = api.get('/repo/0/4', headers=test_user_auth)
        assert rv.status_code == 401

    def test_own_embargo_calc(self, api, example_elastic_calcs, no_warn, other_test_user_auth):
        rv = api.get('/repo/0/4', headers=other_test_user_auth)
        assert rv.status_code == 200

    def test_staging_calc(self, api, example_elastic_calcs, no_warn, test_user_auth):
        rv = api.get('/repo/0/3', headers=test_user_auth)
        assert rv.status_code == 401

    def test_own_staging_calc(self, api, example_elastic_calcs, no_warn, other_test_user_auth):
        rv = api.get('/repo/0/3', headers=other_test_user_auth)
        assert rv.status_code == 200

    def test_non_existing_calcs(self, api, example_elastic_calcs, test_user_auth):
        rv = api.get('/repo/0/10', headers=test_user_auth)
        assert rv.status_code == 404

    def test_search_datasets(self, api, example_elastic_calcs, no_warn, other_test_user_auth):
        rv = api.get('/repo/?owner=all&datasets_grouped=true', headers=other_test_user_auth)
        data = self.assert_search(rv, 4)

        datasets = data.get('datasets_grouped', None)
        assert datasets is not None
        values = datasets['values']
        assert values['ds_id']['total'] == 4
        assert values['ds_id']['examples'][0]['datasets'][0]['dataset_id'] == 'ds_id'
        assert 'after' in datasets
        assert 'datasets' in data['statistics']['total']['all']
        assert data['statistics']['total']['all']['datasets'] > 0

    def test_search_uploads(self, api, example_elastic_calcs, no_warn, other_test_user_auth):
        rv = api.get('/repo/?owner=all&uploads_grouped=true', headers=other_test_user_auth)
        data = self.assert_search(rv, 4)

        uploads = data.get('uploads_grouped', None)
        assert uploads is not None
        values = uploads['values']

        assert values['example_upload_id']['total'] == 3  # just three because of inconsistent upload time
        assert values['example_upload_id']['examples'][0]['upload_id'] == 'example_upload_id'
        assert 'after' in uploads
        assert 'uploads' in data['statistics']['total']['all']

    @pytest.mark.parametrize('calcs, owner, auth', [
        (3, 'all', 'none'),
        (2, 'public', 'none'),
        (3, 'all', 'test_user'),
        (2, 'public', 'test_user'),
        (4, 'all', 'other_test_user'),
        (1, 'user', 'test_user'),
        (3, 'user', 'other_test_user'),
        (0, 'staging', 'test_user'),
        (1, 'staging', 'other_test_user')
    ])
    def test_search_owner(self, api, example_elastic_calcs, no_warn, test_user_auth, other_test_user_auth, calcs, owner, auth):
        auth = dict(none=None, test_user=test_user_auth, other_test_user=other_test_user_auth).get(auth)
        rv = api.get('/repo/?owner=%s' % owner, headers=auth)
        data = self.assert_search(rv, calcs)
        if calcs > 0:
            results = data.get('results', None)
            result = utils.flat(results[0])
            for key in ['uploader.name', 'calc_id', 'formula', 'upload_id']:
                assert key in result

    @pytest.mark.parametrize('calcs, start, end', [
        (2, today - datetime.timedelta(days=6), today),
        (2, today - datetime.timedelta(days=5), today),
        (1, today - datetime.timedelta(days=4), today),
        (1, today, today),
        (1, today - datetime.timedelta(days=6), today - datetime.timedelta(days=5)),
        (0, today - datetime.timedelta(days=7), today - datetime.timedelta(days=6)),
        (2, None, None),
        (1, today, None),
        (2, None, today)
    ])
    def test_search_time(self, api, example_elastic_calcs, no_warn, calcs, start, end):
        query_string = '?owner=visible'
        if start is not None:
            query_string += '&from_time=%s' % rfc3339DateTime.format(start)
        if end is not None:
            query_string += '&until_time=%s' % rfc3339DateTime.format(end)

        rv = api.get('/repo/%s' % query_string)
        self.assert_search(rv, calcs)

    @pytest.mark.parametrize('calcs, quantity, value, user', [
        (2, 'dft.system', 'bulk', 'test_user'),
        (0, 'dft.system', 'atom', 'test_user'),
        (1, 'atoms', 'Br', 'test_user'),
        (1, 'atoms', 'Fe', 'test_user'),
        (0, 'atoms', ['Fe', 'Br', 'A', 'B'], 'test_user'),
        (0, 'only_atoms', ['Br', 'Si'], 'test_user'),
        (1, 'only_atoms', ['Fe'], 'test_user'),
        (1, 'only_atoms', ['Br', 'K', 'Si'], 'test_user'),
        (1, 'only_atoms', ['Br', 'Si', 'K'], 'test_user'),
        (1, 'comment', 'specific', 'test_user'),
        (1, 'authors', 'Leonard Hofstadter', 'test_user'),
        (2, 'files', 'test/mainfile.txt', 'test_user'),
        (2, 'paths', 'mainfile.txt', 'test_user'),
        (2, 'paths', 'test', 'test_user'),
        (2, 'dft.quantities', ['wyckoff_letters_primitive', 'hall_number'], 'test_user'),
        (0, 'dft.quantities', 'dos', 'test_user'),
        (2, 'external_id', 'external_2,external_3', 'other_test_user'),
        (1, 'external_id', 'external_2', 'test_user'),
        (1, 'external_id', 'external_2,external_3', 'test_user'),
        (0, 'external_id', 'external_x', 'test_user')
    ])
    def test_search_parameters(
            self, api, example_elastic_calcs, no_warn, test_user_auth,
            other_test_user_auth, calcs, quantity, value, user):
        user_auth = test_user_auth if user == 'test_user' else other_test_user_auth
        query_string = urlencode({
            'owner': 'public' if user == 'test_user' else 'all',
            quantity: value,
            'statistics': ['dft.system']}, doseq=True)

        rv = api.get('/repo/?%s' % query_string, headers=user_auth)
        logger.debug('run search quantities test', query_string=query_string)
        data = self.assert_search(rv, calcs)

        statistics = data.get('statistics', None)
        assert statistics is not None
        if quantity == 'system' and calcs != 0:
            # for simplicity we only assert on quantities for this case
            assert 'dft.system' in statistics
            assert len(statistics['dft.system']) == 1
            assert value in statistics['dft.system']

    def test_search_exclude(self, api, example_elastic_calcs, no_warn):
        rv = api.get('/repo/?exclude=atoms,only_atoms,dft.optimade,dft.quantities')
        assert rv.status_code == 200
        result = utils.flat(json.loads(rv.data)['results'][0])
        assert 'atoms' not in result
        assert 'only_atoms' not in result
        assert 'dft.basis_set' in result

    metrics_permutations = [[], search_extension.metrics_by_index[config.elastic.index_name]] + [[metric] for metric in search_extension.metrics_by_index[config.elastic.index_name]]

    def test_search_admin(self, api, example_elastic_calcs, no_warn, admin_user_auth):
        rv = api.get('/repo/?owner=admin', headers=admin_user_auth)
        self.assert_search(rv, 4)

    def test_search_admin_auth(self, api, example_elastic_calcs, no_warn, test_user_auth):
        rv = api.get('/repo/?owner=admin', headers=test_user_auth)
        assert rv.status_code == 401

        rv = api.get('/repo/?owner=admin')
        assert rv.status_code == 401

    @pytest.mark.parametrize('metrics, statistics', itertools.product(metrics_permutations, [True, False]))
    def test_search_total_metrics(self, api, example_elastic_calcs, no_warn, metrics, statistics):
        query_params = dict(metrics=metrics)
        if statistics:
            query_params['statistics'] = example_statistics
        url = '/repo/?%s' % urlencode(query_params, doseq=True)
        rv = api.get(url)
        assert rv.status_code == 200, str(rv.data)
        data = json.loads(rv.data)
        total_metrics = data.get('statistics', {}).get('total', {}).get('all', None)
        assert total_metrics is not None
        assert 'code_runs' in total_metrics
        for metric in metrics:
            assert metric in total_metrics

    @pytest.mark.parametrize('metrics', metrics_permutations)
    def test_search_aggregation_metrics(self, api, example_elastic_calcs, no_warn, metrics):
        rv = api.get('/repo/?%s' % urlencode({
            'metrics': metrics,
            'group_statistics': True,
            'dft.groups_grouped': True,
            'encyclopedia.material.materials_grouped': True,
            'datasets_grouped': True,
            'uploads_grouped': True}, doseq=True))

        assert rv.status_code == 200
        data = json.loads(rv.data)

        for name, quantity in data.get('statistics').items():
            for metrics_result in quantity.values():
                assert 'code_runs' in metrics_result
                if name != 'authors':
                    for metric in metrics:
                        assert metric in metrics_result
                else:
                    assert len(metrics_result) == 1  # code_runs is the only metric for authors

        for group in ['dft.groups_grouped', 'uploads_grouped', 'datasets_grouped', 'encyclopedia.material.materials_grouped']:
            assert group in data
            assert 'after' in data[group]
            assert 'values' in data[group]
            # assert len(data[group]['values']) == data['statistics']['total']['all'][group]

    @pytest.mark.parametrize('query, nbuckets', [
        (dict(interval='1M', metrics='dft.total_energies'), 1),
        (dict(interval='1d', metrics='dft.quantities'), 6),
        (dict(interval='1y', from_time='2019-03-20T12:43:54.566414'), 1),
        (dict(until_time='2010-03-20T12:43:54.566414'), 0),
        (dict(interval='1m', from_time='2020-02-20T12:43:54.566414', metrics='dft.calculations'), 7201)
    ])
    def test_search_date_histogram(self, api, example_elastic_calcs, no_warn, query, nbuckets):
        rv = api.get('/repo/?date_histogram=true&%s' % urlencode(query))
        assert rv.status_code == 200
        data = json.loads(rv.data)
        histogram = data.get('statistics').get('date_histogram')
        assert len(histogram) >= nbuckets

    def test_search_date_histogram_metrics(self, api, example_elastic_calcs, no_warn):
        rv = api.get('/repo/?date_histogram=true&metrics=unique_entries')
        assert rv.status_code == 200
        data = json.loads(rv.data)
        histogram = data.get('statistics').get('date_histogram')
        bucket = histogram[list(histogram.keys())[0]]
        assert 'code_runs' in bucket
        assert 'unique_entries' in bucket

    @pytest.mark.parametrize('n_results, page, per_page', [(2, 1, 5), (1, 1, 1), (0, 2, 3)])
    def test_search_pagination(self, api, example_elastic_calcs, no_warn, n_results, page, per_page):
        rv = api.get('/repo/?&page=%d&per_page=%d&statistics=dft.system' % (page, per_page))
        assert rv.status_code == 200
        data = json.loads(rv.data)
        results = data.get('results', None)
        assert data['pagination']['total'] == 2
        assert results is not None
        assert len(results) == n_results

    def test_post_search_query(self, api, example_elastic_calcs, no_warn):
        query_expression = {
            '$and': [
                {'dft.system': 'bulk'},
                {'$not': [{'$lt': {'upload_time': '2020-01-01'}}]}
            ]
        }
        data = {
            'pagination': {'page': 1, 'per_page': 5}, 'query': query_expression,
            'statistics_required': ['dft.system'],
            'uploads_grouped': {}}
        rv = api.post('/repo/', content_type='application/json', data=json.dumps(data))
        assert rv.status_code == 200
        data = json.loads(rv.data)
        results = data.get('results', None)
        assert len(results) > 0
        assert len(data['uploads_grouped']['values']) > 0

    @pytest.mark.parametrize('first, order_by, order', [
        ('1', 'formula', -1), ('2', 'formula', 1),
        ('2', 'dft.basis_set', -1), ('1', 'dft.basis_set', 1),
        (None, 'authors', -1)])
    def test_search_order(self, api, example_elastic_calcs, no_warn, first, order_by, order):
        rv = api.get('/repo/?&order_by=%s&order=%d' % (order_by, order))
        assert rv.status_code == 200
        data = json.loads(rv.data)
        results = data.get('results', None)
        assert data['pagination']['total'] == 2
        assert len(results) == 2
        if first is not None:
            assert results[0]['calc_id'] == first

    @pytest.mark.parametrize('n_results, size', [(2, None), (2, 5), (1, 1)])
    def test_search_scroll(self, api, example_elastic_calcs, no_warn, n_results, size):
        if size is not None:
            rv = api.get('/repo/?scroll=1m,&per_page=%d' % size)
        else:
            rv = api.get('/repo/?scroll=1m')

        assert rv.status_code == 200
        data = json.loads(rv.data)
        results = data.get('results', None)
        assert data.get('scroll', {}).get('size', -1) > 0
        assert results is not None
        assert len(results) == n_results
        scroll_id = data.get('scroll', {}).get('scroll_id', None)
        assert scroll_id is not None

        has_another_page = False
        while scroll_id is not None:
            rv = api.get('/repo/?scroll=1m&scroll_id=%s' % scroll_id)
            data = json.loads(rv.data)
            scroll_info = data.get('scroll', {})
            assert scroll_info['scroll']
            scroll_id = scroll_info.get('scroll_id', None)
            has_another_page |= len(data.get('results')) > 0

        if n_results < 2:
            assert has_another_page

    def test_search_user_authrequired(self, api, example_elastic_calcs, no_warn):
        rv = api.get('/repo/?owner=user')
        assert rv.status_code == 401

    @pytest.mark.parametrize('suggestions, quantity, value', [
        (1, 'dft.system', 'bulk'),
        (1, 'dft.system', 'ulk'),
        (1, 'dft.system', 'ul'),
        (0, 'dft.system', 'notbulk'),
        (1, 'dft.system', None)
    ])
    def test_suggestions_search(self, api, example_elastic_calcs, no_warn, test_user_auth, suggestions, quantity, value):
        url = '/repo/suggestions/%s' % quantity
        if value is not None:
            url = url + '?include=%s' % value
        rv = api.get(url, headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        values = data['suggestions']
        assert len(values) == suggestions

    @pytest.mark.parametrize('calcs, quantity, value', [
        (2, 'dft.system', 'bulk'),
        (0, 'dft.system', 'atom'),
        (1, 'atoms', 'Br'),
        (1, 'atoms', 'Fe'),
        (1, 'authors', 'Leonard Hofstadter'),
        (2, 'files', 'test/mainfile.txt'),
        (0, 'dft.quantities', 'dos'),
        (2, 'dft.compound_type', 'ternary'),
        (0, 'dft.labels_springer_compound_class', 'intermetallic')
    ])
    def test_quantity_search(self, api, example_elastic_calcs, no_warn, test_user_auth, calcs, quantity, value):
        rv = api.get('/repo/quantity/%s' % quantity, headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        values = data['quantity']['values']
        assert (value in values) == (calcs > 0)
        if value in values:
            assert values[value]['total'] == calcs
        else:
            assert 0 == calcs

    def test_quantity_search_after(self, api, example_elastic_calcs, no_warn, test_user_auth):
        rv = api.get('/repo/quantity/atoms?size=1')
        assert rv.status_code == 200
        data = json.loads(rv.data)

        quantity = data['quantity']
        assert 'after' in quantity
        after = quantity['after']
        assert len(quantity['values']) == 1
        value = list(quantity['values'].keys())[0]

        while True:
            rv = api.get('/repo/quantity/atoms?size=1&after=%s' % after)
            assert rv.status_code == 200
            data = json.loads(rv.data)

            quantity = data['quantity']

            if quantity.get('after') is None:
                assert len(quantity['values']) == 0
                break
            assert len(quantity['values']) == 1
            assert value != list(quantity['values'].keys())[0]
            assert after != quantity['after']
            after = quantity['after']

    def test_quantities_search(self, api, example_elastic_calcs, no_warn, test_user_auth):
        rv = api.get(
            '/repo/quantities?%s' % urlencode(
                dict(quantities=['dft.system', 'atoms'], size=1), doseq=True),
            headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert 'dft.system' in data['quantities']
        assert 'atoms' in data['quantities']

    @pytest.mark.parametrize('pid_or_handle, with_login, success', [
        ('2', True, True), ('2', False, True),
        ('3', True, True), ('3', False, False),
        ('4', True, True), ('4', False, False),
        ('21.11132/2', True, True)])
    def test_resolve_pid(
            self, api, example_elastic_calcs, other_test_user_auth, pid_or_handle, with_login,
            success, no_warn):
        rv = api.get(
            '/repo/pid/%s' % pid_or_handle,
            headers=other_test_user_auth if with_login else {})
        assert rv.status_code == 200 if success else 404

        try:
            pid = str(int(pid_or_handle))
        except ValueError:
            pid = str(utils.decode_handle_id(pid_or_handle.split('/')[1]))

        if success:
            assert json.loads(rv.data)['calc_id'] == '%s' % pid
            assert json.loads(rv.data)['upload_id'] == 'example_upload_id'

    @pytest.mark.timeout(config.tests.default_timeout)
    def test_raw_id(self, api, test_user, test_user_auth, proc_infra):
        example_upload_file = 'tests/data/proc/with_raw_id.zip'
        example_upload_id = os.path.basename(example_upload_file).replace('.zip', '')
        test_processing.run_processing((example_upload_id, example_upload_file), test_user)

        rv = api.get(
            '/repo/?%s' % urlencode(dict(owner='all', raw_id='C61A2F88-A0EA-4F0B-AA47-A715868B2E26')),
            headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert data['pagination']['total'] == 1
        assert data['results'][0]['raw_id'] == 'C61A2F88-A0EA-4F0B-AA47-A715868B2E26'

    def test_optimade(self, api, non_empty_processed, test_user_auth):
        rv = api.get(
            '/repo/?%s' % urlencode({'owner': 'all', 'dft.optimade': 'nelements >= 1'}),
            headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert data['pagination']['total'] > 0

        rv = api.get(
            '/repo/?%s' % urlencode({'owner': 'all', 'dft.optimade': 'nelements = 23'}),
            headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert data['pagination']['total'] == 0

        rv = api.get(
            '/repo/?%s' % urlencode({'owner': 'all', 'dft.optimade': 'this is not optimade'}),
            headers=test_user_auth)
        assert rv.status_code == 400

    def test_labels(self, api, non_empty_processed, test_user_auth):
        rv = api.get(
            '/repo/?%s' % urlencode(dict(owner='all', labels=['nonmetal', 'semiconductor']), doseq=True),
            headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert data['pagination']['total'] > 0

    def test_get_code_from_query(self, api, example_elastic_calcs, test_user_auth):
        rv = api.get('/repo/?code_name=VASP', headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert data['code']['curl'] is not None
        assert data['code']['python'] is not None
        assert data['code']['clientlib'] is not None


class TestEditRepo():

    @pytest.fixture(autouse=True)
    def class_api(self, session_client, elastic_infra, mongo_infra):
        clear_elastic(elastic_infra)
        mongo_infra.drop_database('test_db')

        self.api = BlueprintClient(session_client, '/api')
        yield

        mongo_infra.drop_database('test_db')
        clear_elastic(elastic_infra)

    @pytest.fixture(autouse=True)
    def example_datasets(self, test_user, other_test_user):
        self.example_dataset = Dataset(
            dataset_id='example_ds', name='example_ds', user_id=test_user.user_id)
        self.example_dataset.a_mongo.create()

        self.other_example_dataset = Dataset(
            dataset_id='other_example_ds', name='other_example_ds',
            user_id=other_test_user.user_id)
        self.other_example_dataset.a_mongo.create()

        yield

        self.example_dataset.a_mongo.delete()
        self.other_example_dataset.a_mongo.delete()

    @pytest.fixture(autouse=True)
    def remove_new_dataset(self):
        yield 'new_ds'
        Dataset.m_def.a_mongo.objects(name='new_ds').delete()

    @pytest.fixture(autouse=True)
    def example_data(self, class_api, test_user, other_test_user, raw_files):
        # TODO
        example_data = ExampleData()

        example_data.create_entry(
            upload_id='upload_1', uploader=test_user, published=True, with_embargo=False)
        example_data.create_entry(
            upload_id='upload_2', uploader=test_user, published=True, with_embargo=True)
        example_data.create_entry(
            upload_id='upload_2', uploader=test_user, published=False, with_embargo=False)
        example_data.create_entry(
            upload_id='upload_3', uploader=other_test_user, published=True, with_embargo=False)

        example_data.save()

        # for i, entry in enumerate(entries):
        #     upload = uploads[entry.pop('upload_id')]
        #     user = entry.pop('user')
        #     metadata = dict(uploader=user.user_id, **entry)
        #     upload.create_test_structure(i + 1, 2, 1, [], 0, metadata=metadata)

        # for upload in uploads.values():
        #     upload.create_upload_files()

        # search_refresh()

    @pytest.fixture(autouse=True)
    def auth(self, test_user_auth):
        self.test_user_auth = test_user_auth

    def perform_edit(self, query=None, verify=False, **kwargs):
        actions = {}
        for key, value in kwargs.items():
            if isinstance(value, list):
                actions[key] = [dict(value=i) for i in value]
            else:
                actions[key] = dict(value=value)

        data = dict(actions=actions)
        if query is not None:
            data.update(query=query)
        if verify:
            data.update(verify=verify)

        return self.api.post(
            '/repo/edit', headers=self.test_user_auth, content_type='application/json',
            data=json.dumps(data))

    def assert_edit(self, rv, quantity: str, success: bool, message: bool, status_code: int = 200):
        assert rv.status_code == status_code, rv.data
        data = json.loads(rv.data)
        actions = data.get('actions')
        assert actions is not None
        assert [quantity] == list(actions.keys())
        quantity_actions = actions[quantity]
        if not isinstance(quantity_actions, list):
            quantity_actions = [quantity_actions]
        has_failure = False
        has_message = False
        for action in quantity_actions:
            has_failure = has_failure or not action['success']
            has_message = has_message or ('message' in action)
        assert not has_failure == success
        assert has_message == message

    def mongo(self, *args, edited: bool = True, **kwargs):
        for calc_id in args:
            calc = Calc.objects(calc_id='test_entry_id_%d' % calc_id).first()
            assert calc is not None
            metadata = calc.metadata
            if edited:
                assert metadata.get('last_edit') is not None
            for key, value in kwargs.items():
                if metadata.get(key) != value:
                    return False
        return True

    def assert_elastic(self, *args, invert: bool = False, **kwargs):
        def assert_entry(get_entries):
            for arg in args:
                entry_id = 'test_entry_id_%d' % arg
                entries = list(get_entries(entry_id))
                assert len(entries) > 0, entry_id
                for entry in entries:
                    for key, value in kwargs.items():
                        if key in ['authors', 'owners']:
                            ids = [user['user_id'] for user in entry.get(key)]
                            if ids != value:
                                return False
                        else:
                            if entry.get(key) != value:
                                return False

            return True

        # test v0 data
        assert invert != assert_entry(
            lambda id: SearchRequest().search_parameters(calc_id=id).execute_scan())

        # test v1 data
        assert invert != assert_entry(
            lambda id: search(owner=None, query=dict(entry_id=id)).data)

    def test_edit_all_properties(self, test_user, other_test_user):
        edit_data = dict(
            comment='test_edit_props',
            references=['http://test', 'http://test2'],
            coauthors=[other_test_user.user_id],
            shared_with=[other_test_user.user_id])
        rv = self.perform_edit(**edit_data, query=dict(upload_id=['upload_1']))
        result = json.loads(rv.data)
        assert rv.status_code == 200
        actions = result.get('actions')
        for key in edit_data:
            assert key in actions
            quantity_actions = actions.get(key)
            if not isinstance(quantity_actions, list):
                quantity_actions = [quantity_actions]
            for quantity_action in quantity_actions:
                assert quantity_action['success']

        assert self.mongo(1, comment='test_edit_props')
        assert self.mongo(1, references=['http://test', 'http://test2'])
        assert self.mongo(1, coauthors=[other_test_user.user_id])
        assert self.mongo(1, shared_with=[other_test_user.user_id])

        self.assert_elastic(1, comment='test_edit_props')
        self.assert_elastic(1, references=['http://test', 'http://test2'])
        self.assert_elastic(1, authors=[test_user.user_id, other_test_user.user_id])
        self.assert_elastic(1, owners=[test_user.user_id, other_test_user.user_id])

        edit_data = dict(
            comment='',
            references=[],
            coauthors=[],
            shared_with=[])
        rv = self.perform_edit(**edit_data, query=dict(upload_id=['upload_1']))
        result = json.loads(rv.data)
        assert rv.status_code == 200
        actions = result.get('actions')
        for key in edit_data:
            assert key in actions
            quantity_actions = actions.get(key)
            if not isinstance(quantity_actions, list):
                quantity_actions = [quantity_actions]
            for quantity_action in quantity_actions:
                assert quantity_action['success']

        assert self.mongo(1, comment=None)
        assert self.mongo(1, references=[])
        assert self.mongo(1, coauthors=[])
        assert self.mongo(1, shared_with=[])

        self.assert_elastic(1, comment=None)
        self.assert_elastic(1, references=[])
        self.assert_elastic(1, authors=[test_user.user_id])
        self.assert_elastic(1, owners=[test_user.user_id])

    def test_edit_all(self):
        rv = self.perform_edit(comment='test_edit_all')
        self.assert_edit(rv, quantity='comment', success=True, message=False)
        assert self.mongo(1, 2, 3, comment='test_edit_all')
        self.assert_elastic(1, 2, 3, comment='test_edit_all')
        assert not self.mongo(4, comment='test_edit_all', edited=False)
        self.assert_elastic(4, comment='test_edit_all', edited=False, invert=True)

    def test_edit_multi(self):
        rv = self.perform_edit(comment='test_edit_multi', query=dict(upload_id=['upload_1', 'upload_2']))
        self.assert_edit(rv, quantity='comment', success=True, message=False)
        assert self.mongo(1, 2, 3, comment='test_edit_multi')
        self.assert_elastic(1, 2, 3, comment='test_edit_multi')
        assert not self.mongo(4, comment='test_edit_multi', edited=False)
        self.assert_elastic(4, comment='test_edit_multi', edited=False, invert=True)

    def test_edit_some(self):
        rv = self.perform_edit(comment='test_edit_some', query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, quantity='comment', success=True, message=False)
        assert self.mongo(1, comment='test_edit_some')
        self.assert_elastic(1, comment='test_edit_some')
        assert not self.mongo(2, 3, 4, comment='test_edit_some', edited=False)
        self.assert_elastic(2, 3, 4, comment='test_edit_some', edited=False, invert=True)

    def test_edit_verify(self):
        rv = self.perform_edit(
            comment='test_edit_verify', verify=True, query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, quantity='comment', success=True, message=False)
        assert not self.mongo(1, comment='test_edit_verify', edited=False)

    def test_edit_empty_list(self, other_test_user):
        rv = self.perform_edit(coauthors=[other_test_user.user_id], query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, quantity='coauthors', success=True, message=False)
        rv = self.perform_edit(coauthors=[], query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, quantity='coauthors', success=True, message=False)
        assert self.mongo(1, coauthors=[])

    def test_edit_duplicate_value(self, other_test_user):
        rv = self.perform_edit(coauthors=[other_test_user.user_id, other_test_user.user_id], query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, status_code=400, quantity='coauthors', success=False, message=True)

    def test_edit_uploader_as_coauthor(self, test_user):
        rv = self.perform_edit(coauthors=[test_user.user_id], query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, status_code=400, quantity='coauthors', success=False, message=True)

    def test_edit_ds(self):
        rv = self.perform_edit(
            datasets=[self.example_dataset.name], query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, quantity='datasets', success=True, message=False)
        assert self.mongo(1, datasets=[self.example_dataset.dataset_id])

    def test_edit_ds_remove_doi(self):
        rv = self.perform_edit(
            datasets=[self.example_dataset.name], query=dict(upload_id=['upload_1']))
        assert rv.status_code == 200
        rv = self.api.post('/datasets/%s' % self.example_dataset.name, headers=self.test_user_auth)
        assert rv.status_code == 200
        rv = self.perform_edit(datasets=[], query=dict(upload_id=['upload_1']))
        assert rv.status_code == 400
        data = json.loads(rv.data)
        assert not data['success']
        assert self.example_dataset.name in data['message']
        assert Dataset.m_def.a_mongo.get(dataset_id=self.example_dataset.dataset_id) is not None

    def test_edit_ds_remove(self):
        rv = self.perform_edit(
            datasets=[self.example_dataset.name], query=dict(upload_id=['upload_1']))
        assert rv.status_code == 200
        rv = self.perform_edit(datasets=[], query=dict(upload_id=['upload_1']))
        assert rv.status_code == 200
        with pytest.raises(KeyError):
            assert Dataset.m_def.a_mongo.get(dataset_id=self.example_dataset.dataset_id) is None

    def test_edit_ds_user_namespace(self, test_user):
        assert Dataset.m_def.a_mongo.objects(
            name=self.other_example_dataset.name).first() is not None

        rv = self.perform_edit(
            datasets=[self.other_example_dataset.name], query=dict(upload_id=['upload_1']))

        self.assert_edit(rv, quantity='datasets', success=True, message=True)
        new_dataset = Dataset.m_def.a_mongo.objects(
            name=self.other_example_dataset.name,
            user_id=test_user.user_id).first()
        assert new_dataset is not None
        assert self.mongo(1, datasets=[new_dataset.dataset_id])

    def test_edit_new_ds(self, test_user):
        rv = self.perform_edit(datasets=['new_dataset'], query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, quantity='datasets', success=True, message=True)
        new_dataset = Dataset.m_def.a_mongo.objects(name='new_dataset').first()
        assert new_dataset is not None
        assert new_dataset.user_id == test_user.user_id
        assert self.mongo(1, datasets=[new_dataset.dataset_id])

    def test_edit_bad_user(self):
        rv = self.perform_edit(coauthors=['bad_user'], query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, status_code=400, quantity='coauthors', success=False, message=True)

    def test_edit_user(self, other_test_user):
        rv = self.perform_edit(coauthors=[other_test_user.user_id], query=dict(upload_id=['upload_1']))
        self.assert_edit(rv, quantity='coauthors', success=True, message=False)

    def test_admin_only(self, other_test_user):
        rv = self.perform_edit(uploader=other_test_user.user_id)
        assert rv.status_code != 200


@pytest.mark.timeout(config.tests.default_timeout)
def test_edit_lift_embargo(api, published, other_test_user_auth, no_warn):
    example_calc = Calc.objects(upload_id=published.upload_id).first()
    assert example_calc.metadata['with_embargo']
    elastic_calc = next(
        search.SearchRequest().search_parameters(calc_id=example_calc.calc_id).execute_scan())
    assert elastic_calc['with_embargo'] is True
    with pytest.raises(files.Restricted):
        with files.UploadFiles.get(published.upload_id).read_archive(example_calc.calc_id) as archive:
            archive[example_calc.calc_id].to_dict()

    rv = api.post(
        '/repo/edit', headers=other_test_user_auth, content_type='application/json',
        data=json.dumps({
            'actions': {
                'with_embargo': {}
            }
        }))

    assert rv.status_code == 200, rv.data
    assert not Calc.objects(calc_id=example_calc.calc_id).first().metadata['with_embargo']

    Upload.get(published.upload_id).block_until_complete()

    elastic_calc = next(
        search.SearchRequest().search_parameters(calc_id=example_calc.calc_id).execute_scan())
    assert elastic_calc['with_embargo'] is False

    # should not raise Restricted anymore
    with files.UploadFiles.get(published.upload_id).read_archive(example_calc.calc_id) as archive:
        archive[example_calc.calc_id].to_dict()


@pytest.mark.timeout(config.tests.default_timeout)
def test_edit_lift_embargo_unnecessary(api, published_wo_user_metadata, other_test_user_auth):
    example_calc = Calc.objects(upload_id=published_wo_user_metadata.upload_id).first()
    assert not example_calc.metadata['with_embargo']
    rv = api.post(
        '/repo/edit', headers=other_test_user_auth, content_type='application/json',
        data=json.dumps({
            'actions': {
                'with_embargo': {
                    'value': 'lift'
                }
            }
        }))
    assert rv.status_code == 400
    data = json.loads(rv.data)
    assert not data['actions']['with_embargo']['success']


class TestRaw(UploadFilesBasedTests):

    def test_raw_file_from_calc(self, api, non_empty_processed, test_user_auth):
        calc = list(non_empty_processed.calcs)[0]
        url = '/raw/calc/%s/%s/%s' % (
            non_empty_processed.upload_id, calc.calc_id, os.path.basename(calc.mainfile))
        rv = api.get(url, headers=test_user_auth)
        assert rv.status_code == 200
        assert len(rv.data) > 0

        url = '/raw/calc/%s/%s/' % (non_empty_processed.upload_id, calc.calc_id)
        rv = api.get(url, headers=test_user_auth)
        assert rv.status_code == 200
        result = json.loads(rv.data)
        assert len(result['contents']) > 0

        url = '/raw/calc/%s/not_existing/' % (non_empty_processed.upload_id)
        rv = api.get(url, headers=test_user_auth)
        assert rv.status_code == 404

    @UploadFilesBasedTests.check_authorization
    def test_raw_file(self, api, upload, auth_headers):
        url = '/raw/%s/%s' % (upload, example_file_mainfile)
        rv = api.get(url, headers=auth_headers)
        assert rv.status_code == 200
        assert len(rv.data) > 0

    @UploadFilesBasedTests.check_authorization
    def test_raw_file_partial(self, api, upload, auth_headers):
        url = '/raw/%s/%s?offset=0&length=20' % (upload, example_file_mainfile)
        rv = api.get(url, headers=auth_headers)
        assert rv.status_code == 200
        start_data = rv.data
        assert len(start_data) == 20

        url = '/raw/%s/%s?offset=10&length=10' % (upload, example_file_mainfile)
        rv = api.get(url, headers=auth_headers)
        assert rv.status_code == 200
        next_data = rv.data
        assert len(rv.data) == 10
        assert start_data[10:] == next_data

    def test_raw_file_compressed(self, api, raw_files, admin_user_auth):
        upload = files.ArchiveBasedStagingUploadFiles(
            'upload_id', upload_path='tests/data/api/example_with_compressed.zip', create=True)
        upload.extract()
        for compression in ['gz', 'xz']:
            rv = api.get(
                'raw/upload_id/example_with_compressed/mainfile.%s?'
                'decompress=true&offset=5&length=3' % compression,
                headers=admin_user_auth)
            assert rv.status_code == 200
            assert rv.data == b'con'

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_signed(self, api, upload, _, test_user_signature_token):
        url = '/raw/%s/%s?signature_token=%s' % (upload, example_file_mainfile, test_user_signature_token)
        rv = api.get(url)
        assert rv.status_code == 200
        assert len(rv.data) > 0

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_missing_file(self, api, upload, auth_headers):
        url = '/raw/%s/does/not/exist' % upload
        rv = api.get(url, headers=auth_headers)
        assert rv.status_code == 404
        data = json.loads(rv.data)
        assert 'files' not in data

    @pytest.mark.parametrize('compress', [True, False])
    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_wildcard(self, api, upload, auth_headers, compress):
        url = '/raw/%s/examples*' % upload
        if compress:
            url = '%s?compress=1' % url
        rv = api.get(url, headers=auth_headers)

        assert rv.status_code == 200
        assert_zip_file(rv, files=len(example_file_contents))

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_wildcard_missing(self, api, upload, auth_headers):
        url = '/raw/%s/does/not/exist*' % upload
        rv = api.get(url, headers=auth_headers)
        assert rv.status_code == 404

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_missing_upload(self, api, upload, auth_headers):
        url = '/raw/doesnotexist/%s' % example_file_mainfile
        rv = api.get(url, headers=auth_headers)
        assert rv.status_code == 404

    @pytest.mark.parametrize('compress, strip', [(True, False), (False, False), (False, True)])
    @UploadFilesBasedTests.check_authorization
    def test_raw_files(self, api, upload, auth_headers, compress, strip):
        url = '/raw/%s?files=%s' % (
            upload, ','.join(example_file_contents))
        if compress:
            url = '%s&compress=1' % url
        if strip:
            url = '%s&strip=1' % url
        rv = api.get(url, headers=auth_headers)

        assert rv.status_code == 200
        assert_zip_file(rv, files=len(example_file_contents), basename=strip)

    @pytest.mark.parametrize('compress', [False, True])
    def test_raw_files_from_query_upload_id(self, api, non_empty_processed, test_user_auth, compress):
        url = '/raw/query?upload_id=%s&compress=%s' % (non_empty_processed.upload_id, 'true' if compress else 'false')
        rv = api.get(url, headers=test_user_auth)

        assert rv.status_code == 200
        assert_zip_file(rv, files=len(example_file_contents) + 1)

    @pytest.mark.parametrize('query_params', [
        {'atoms': 'Si'},
        {'authors': 'Sheldon Cooper'}
    ])
    def test_raw_files_from_query(self, api, processeds, test_user_auth, query_params):

        url = '/raw/query?%s' % urlencode(query_params)
        rv = api.get(url, headers=test_user_auth)

        assert rv.status_code == 200
        assert_zip_file(rv, files=len(example_file_contents) * len(processeds) + 1)
        with zipfile.ZipFile(io.BytesIO(rv.data)) as zip_file:
            with zip_file.open('manifest.json', 'r') as f:
                manifest = json.load(f)
                assert len(manifest) == len(processeds)

    def test_raw_files_from_empty_query(self, api, elastic):
        url = '/raw/query?upload_id=doesNotExist'
        rv = api.get(url)

        assert rv.status_code == 200
        assert_zip_file(rv, files=1)

    def test_raw_files_from_query_file_error(self, api, processeds, test_user_auth, monkeypatch):

        def raw_file(self, *args, **kwargs):
            raise Exception('test error')

        monkeypatch.setattr('nomad.files.StagingUploadFiles.raw_file', raw_file)
        url = '/raw/query?%s' % urlencode({'atoms': 'Si'})
        rv = api.get(url, headers=test_user_auth)

        assert rv.status_code == 200
        with zipfile.ZipFile(io.BytesIO(rv.data)) as zip_file:
            with zip_file.open('manifest.json', 'r') as f:
                manifest = json.load(f)
                assert 'errors' in manifest
                assert len(manifest['errors']) > 0

    @pytest.mark.parametrize('files, pattern, strip', [
        (1, '*.json', False),
        (1, '*.json', True),
        (5, ['*.json', '*.aux'], False)])
    def test_raw_query_pattern(self, api, non_empty_processed, test_user_auth, files, pattern, strip):
        params = dict(file_pattern=pattern)
        if strip:
            params.update(strip=True)
        url = '/raw/query?%s' % urlencode(params, doseq=True)
        rv = api.get(url, headers=test_user_auth)
        assert rv.status_code == 200
        assert_zip_file(rv, files=(files + 1), basename=strip)

    def test_post_raw_query(self, api, processeds, test_user_auth):
        data = {
            'query': {'$or': [{'atoms': ['Si']}, {'authors': ['Sheldon Cooper']}]},
            'file_pattern': ['*.json', '*.aux'], 'strip': False
        }
        rv = api.post(
            '/raw/query', content_type='application/json', data=json.dumps(data),
            headers=test_user_auth)

        assert rv.status_code == 200
        assert_zip_file(rv, files=len(example_file_contents) * len(processeds) + 1)

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_files_signed(self, api, upload, _, test_user_signature_token):
        url = '/raw/%s?files=%s&signature_token=%s' % (
            upload, ','.join(example_file_contents), test_user_signature_token)
        rv = api.get(url)

        assert rv.status_code == 200
        assert_zip_file(rv, files=len(example_file_contents))

    @pytest.mark.parametrize('compress', [True, False, None])
    @UploadFilesBasedTests.check_authorization
    def test_raw_files_post(self, api, upload, auth_headers, compress):
        url = '/raw/%s' % upload
        data = dict(files=example_file_contents)
        if compress is not None:
            data.update(compress=compress)
        rv = api.post(url, data=json.dumps(data), content_type='application/json', headers=auth_headers)

        assert rv.status_code == 200
        assert_zip_file(rv, files=len(example_file_contents))

    @pytest.mark.parametrize('compress', [True, False])
    @UploadFilesBasedTests.ignore_authorization
    def test_raw_files_missing_file(self, api, upload, auth_headers, compress):
        url = '/raw/%s?files=%s,missing/file.txt' % (upload, example_file_mainfile)
        if compress:
            url = '%s&compress=1' % url
        rv = api.get(url, headers=auth_headers)

        assert rv.status_code == 200
        assert_zip_file(rv, files=1)

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_files_missing_upload(self, api, upload, auth_headers):
        url = '/raw/doesnotexist?files=shoud/not/matter.txt'
        rv = api.get(url, headers=auth_headers)

        assert rv.status_code == 404

    @pytest.mark.parametrize('path', ['examples_template', 'examples_template/'])
    @UploadFilesBasedTests.ignore_authorization
    def test_raw_files_list(self, api, upload, auth_headers, path):
        url = '/raw/%s/%s' % (upload, path)
        rv = api.get(url, headers=auth_headers)
        assert rv.status_code == 200
        data = json.loads(rv.data)

        assert len(data['contents']) == 5
        assert data['upload_id'] == upload
        assert data['directory'] == 'examples_template'
        for content in data['contents']:
            assert content['name'] is not None
            assert content['size'] >= 0
        assert '1.aux' in list(content['name'] for content in data['contents'])

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_files_list_missing(self, api, upload, auth_headers):
        url = '/raw/%s/examples_' % upload
        rv = api.get(url, headers=auth_headers)
        assert rv.status_code == 404


class TestMirror:

    def test_upload(self, api, published, admin_user_auth, no_warn):
        url = '/mirror/%s' % published.upload_id
        rv = api.get(url, headers=admin_user_auth)
        assert rv.status_code == 200

        data = json.loads(rv.data)
        assert data['upload_id'] == published.upload_id
        assert json.loads(data['upload'])['_id'] == published.upload_id
        assert len(data['calcs']) == len(published.calcs)
        assert data['upload_files_path'] == published.upload_files.os_path

    def test_uploads(self, api, published, admin_user_auth, no_warn):
        rv = api.post(
            '/mirror/',
            content_type='application/json', data='{"query":{}}', headers=admin_user_auth)
        assert rv.status_code == 200, rv.data

        data = json.loads(rv.data)
        assert data[0]['upload_id'] == published.upload_id

    @pytest.mark.parametrize('with_doi', [False, True])
    def test_dataset(self, api, published_wo_user_metadata, admin_user_auth, test_user_auth, with_doi):
        rv = api.post(
            '/repo/edit', headers=test_user_auth, content_type='application/json',
            data=json.dumps({
                'actions': {
                    'datasets': [{
                        'value': 'test_dataset'
                    }]
                }
            }))
        assert rv.status_code == 200

        if with_doi:
            rv = api.post('/datasets/test_dataset', headers=test_user_auth)
            assert rv.status_code == 200

        rv = api.post(
            '/mirror/',
            content_type='application/json', data='{"query":{}}', headers=admin_user_auth)
        assert rv.status_code == 200, rv.data

        url = '/mirror/%s' % published_wo_user_metadata.upload_id
        rv = api.get(url, headers=admin_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert len(data['datasets']) == 1
        dataset = data['calcs'][0]['metadata']['datasets'][0]
        assert dataset in data['datasets']
        if with_doi:
            assert len(data['dois']) == 1
            assert data['datasets'][dataset]['doi'] is not None
            assert data['datasets'][dataset]['doi'] in data['dois']
        else:
            assert 'dois' not in data

    # TODO
    # - parametrize to also check raw
    # - compute the hex digest reference
    def test_files(self, api, published, admin_user_auth, no_warn):
        url = '/mirror/files/%s?prefix=archive' % published.upload_id
        rv = api.get(url, headers=admin_user_auth)
        assert rv.status_code == 200
        assert rv.data is not None
        assert md5(rv.data).hexdigest() == 'a50a980a4f1bd9892e95410936a36cdf'

    def test_users(self, api, published, admin_user_auth, no_warn):
        url = '/mirror/users'
        rv = api.get(url, headers=admin_user_auth)
        assert rv.status_code == 200
        users = json.loads(rv.data)
        assert users is not None
        assert len(users) == 3
        for user in users:
            assert 'email' not in user


class TestDataset:

    @pytest.fixture()
    def example_datasets(self, mongo, test_user):
        Dataset(dataset_id='1', user_id=test_user.user_id, name='ds1').a_mongo.create()
        Dataset(dataset_id='2', user_id=test_user.user_id, name='ds2', doi='test_doi').a_mongo.create()
        Dataset(dataset_id='3', user_id=test_user.user_id, name='weird+/*?& name').a_mongo.create()

    def assert_dataset(self, dataset, name: str = None, doi: bool = False):
        assert 'dataset_id' in dataset
        assert 'user_id' in dataset
        assert ('doi' in dataset) == doi
        assert dataset.get('name') is not None
        if name is not None:
            assert dataset.get('name') == name

    def assert_dataset_entry(self, api, dataset_id: str, exists: bool, with_doi: bool, **kwargs):
        rv = api.get('/repo/?dataset_id=%s' % dataset_id, **kwargs)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        total = data['pagination']['total']
        if exists:
            assert total > 0
        else:
            assert total == 0

        if exists:
            doi = data['results'][0]['datasets'][0]['doi']
            if with_doi:
                assert doi is not None
            else:
                assert doi is None

    def test_create_dataset(self, api, test_user_auth):
        rv = api.put(
            '/datasets/', headers=test_user_auth,
            data=json.dumps(dict(name='test_dataset')),
            content_type='application/json')
        assert rv.status_code == 200
        data = json.loads(rv.data)
        self.assert_dataset(data, 'test_dataset')

    @pytest.mark.parametrize('data', [
        dict(name='test_name', doi='something'),
        dict(name='test_name', dataset_id='something'),
        dict(name='test_name', user_id='something'),
        dict(name='test_name', unknown_key='something'),
        dict()])
    def test_create_dataset_bad_data(self, api, test_user_auth, data):
        rv = api.put(
            '/datasets/', headers=test_user_auth,
            data=json.dumps(data),
            content_type='application/json')
        assert rv.status_code >= 400

    def test_get_datasets(self, api, test_user_auth, example_datasets):
        rv = api.get('/datasets/', headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert 'pagination' in data
        assert data['pagination']['total'] == 3
        assert len(data['results']) == 3
        for dataset in data['results']:
            if dataset['name'] == 'ds2':
                self.assert_dataset(dataset, doi=True)
            else:
                self.assert_dataset(dataset)

    def test_get_datasets_prefix(self, api, test_user_auth, example_datasets):
        rv = api.get('/datasets/?prefix=ds1', headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert 'pagination' in data
        assert data['pagination']['total'] == 1
        assert len(data['results']) == 1
        assert data['results'][0]['name'] == 'ds1'

    def test_get_dataset(self, api, test_user_auth, example_datasets):
        rv = api.get('/datasets/ds1', headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        self.assert_dataset(data, name='ds1')

    def test_get_dataset_missing(self, api, other_test_user_auth, example_datasets):
        rv = api.get('/datasets/ds1', headers=other_test_user_auth)
        assert rv.status_code == 404

    @pytest.fixture()
    def example_dataset_with_entry(self, mongo, elastic, raw_files, example_datasets):
        example_data = ExampleData()
        example_entry = example_data.create_entry(
            calc_id='1', upload_id='1', published=True, with_embargo=False,
            datasets=['1'])
        example_data.save()
        from nomad import files
        assert files.UploadFiles.get('1') is not None
        assert example_entry.section_metadata.datasets[0].dataset_id == '1'

        # entry_archive = EntryArchive()
        # entry_metadata = entry_archive.m_create(
        #     EntryMetadata,
        #     domain='dft', calc_id='1', upload_id='1', published=True, with_embargo=False,
        #     datasets=['1'])
        # Calc(
        #     calc_id='1', upload_id='1', create_time=datetime.datetime.now(),
        #     metadata=entry_metadata.m_to_dict()).save()
        # upload_files = files.StagingUploadFiles(upload_id='1', create=True)
        # upload_files.write_archive('1', dict(section_metadata=entry_metadata.m_to_dict()))
        # upload_files.close()
        # index(entry_archive)

    def test_delete_dataset(self, api, test_user_auth, example_dataset_with_entry):
        rv = api.delete('/datasets/ds1', headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        self.assert_dataset(data, name='ds1')
        api.get('/datasets/ds1', headers=test_user_auth).status_code == 404
        self.assert_dataset_entry(api, '1', False, False, headers=test_user_auth)

    def test_get_dataset_with_doi(self, api, test_user_auth, example_datasets):
        rv = api.delete('/datasets/ds2', headers=test_user_auth)
        assert rv.status_code == 400

    @pytest.mark.parametrize('dataset_name,has_entry', [('ds1', True), ('weird+/*?& name', False)])
    def test_assign_doi(self, api, test_user_auth, example_dataset_with_entry, dataset_name, has_entry):
        import urllib.parse
        rv = api.post('/datasets/%s' % urllib.parse.quote(dataset_name), headers=test_user_auth)
        if has_entry:
            assert rv.status_code == 200
            data = json.loads(rv.data)
            self.assert_dataset(data, name=dataset_name, doi=True)
            self.assert_dataset_entry(api, '1', True, True, headers=test_user_auth)
        else:
            assert rv.status_code == 400

    def test_assign_doi_empty(self, api, test_user_auth, example_datasets):
        rv = api.post('/datasets/ds1', headers=test_user_auth)
        assert rv.status_code == 400

    def test_assign_doi_unpublished(self, api, test_user_auth, example_datasets):
        entry_metadata = EntryMetadata(
            domain='dft', calc_id='1', upload_id='1', published=False, with_embargo=False,
            datasets=['1'])
        Calc(
            calc_id='1', upload_id='1', create_time=datetime.datetime.now(),
            metadata=entry_metadata.m_to_dict()).save()
        rv = api.post('/datasets/ds1', headers=test_user_auth)
        assert rv.status_code == 400

    def test_resolve_doi(self, api, example_datasets):
        rv = api.get('/datasets/doi/test_doi')
        assert rv.status_code == 200
        data = json.loads(rv.data)
        self.assert_dataset(data, name='ds2', doi=True)
