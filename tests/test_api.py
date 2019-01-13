# Copyright 2018 Markus Scheidgen
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

import pytest
import time
import json
from mongoengine import connect
from mongoengine.connection import disconnect
import base64
import zipfile
import io
import inspect

from nomad import config
# for convinience we test the api without path prefix
services_config = config.services._asdict()
services_config.update(api_base_path='')
config.services = config.NomadServicesConfig(**services_config)

from nomad import api, coe_repo  # noqa
from nomad.files import UploadFiles, PublicUploadFiles  # noqa
from nomad.processing import Upload, Calc  # noqa
from nomad.coe_repo import User  # noqa

from tests.processing.test_data import example_files  # noqa
from tests.test_files import example_file, example_file_mainfile, example_file_contents  # noqa
from tests.test_files import create_staging_upload, create_public_upload  # noqa

# import fixtures
from tests.test_normalizing import normalized_template_example  # noqa pylint: disable=unused-import
from tests.test_parsing import parsed_template_example  # noqa pylint: disable=unused-import
# from tests.test_repo import example_elastic_calc  # noqa pylint: disable=unused-import
from tests.test_coe_repo import assert_coe_upload  # noqa


@pytest.fixture(scope='function')
def client(mockmongo):
    disconnect()
    connect('users_test', host=config.mongo.host, port=config.mongo.port, is_mock=True)

    api.app.config['TESTING'] = True
    client = api.app.test_client()

    yield client
    Upload._get_collection().drop()


def create_auth_headers(user):
    basic_auth_str = '%s:password' % user.email
    basic_auth_bytes = basic_auth_str.encode('utf-8')
    basic_auth_base64 = base64.b64encode(basic_auth_bytes).decode('utf-8')
    return {
        'Authorization': 'Basic %s' % basic_auth_base64
    }


@pytest.fixture(scope='session')
def test_user_auth(test_user: User):
    return create_auth_headers(test_user)


@pytest.fixture(scope='session')
def test_other_user_auth(other_test_user: User):
    return create_auth_headers(other_test_user)


@pytest.fixture(scope='session')
def admin_user_auth(admin_user: User):
    return create_auth_headers(admin_user)


class TestAdmin:

    @pytest.mark.timeout(10)
    def test_reset(self, client, admin_user_auth, repository_db):
        rv = client.post('/admin/reset', headers=admin_user_auth)
        assert rv.status_code == 200

    # TODO disabled as this will destroy the session repository_db beyond repair.
    # @pytest.mark.timeout(10)
    # def test_remove(self, client, admin_user_auth, repository_db):
    #     rv = client.post('/admin/remove', headers=admin_user_auth)
    #     assert rv.status_code == 200

    def test_doesnotexist(self, client, admin_user_auth):
        rv = client.post('/admin/doesnotexist', headers=admin_user_auth)
        assert rv.status_code == 404

    def test_only_admin(self, client, test_user_auth):
        rv = client.post('/admin/doesnotexist', headers=test_user_auth)
        assert rv.status_code == 401

    @pytest.fixture(scope='function')
    def disable_reset(self, monkeypatch):
        old_config = config.services
        new_config = config.NomadServicesConfig(
            config.services.api_host,
            config.services.api_port,
            config.services.api_base_path,
            config.services.api_secret,
            config.services.admin_password,
            True)
        monkeypatch.setattr(config, 'services', new_config)
        yield None
        monkeypatch.setattr(config, 'services', old_config)

    def test_disabled(self, client, admin_user_auth, disable_reset):
        rv = client.post('/admin/reset', headers=admin_user_auth)
        assert rv.status_code == 400


class TestAuth:
    def test_xtoken_auth(self, client, test_user: User, no_warn):
        rv = client.get('/uploads/', headers={
            'X-Token': test_user.email  # the test users have their email as tokens for convinience
        })

        assert rv.status_code == 200

    def test_xtoken_auth_denied(self, client, no_warn):
        rv = client.get('/uploads/', headers={
            'X-Token': 'invalid'
        })

        assert rv.status_code == 401

    def test_basic_auth(self, client, test_user_auth, no_warn):
        rv = client.get('/uploads/', headers=test_user_auth)
        assert rv.status_code == 200

    def test_basic_auth_denied(self, client, no_warn):
        basic_auth_base64 = base64.b64encode('invalid'.encode('utf-8')).decode('utf-8')
        rv = client.get('/uploads/', headers={
            'Authorization': 'Basic %s' % basic_auth_base64
        })
        assert rv.status_code == 401

    def test_get_token(self, client, test_user_auth, test_user: User, no_warn):
        rv = client.get('/auth/token', headers=test_user_auth)
        assert rv.status_code == 200
        assert rv.data.decode('utf-8') == test_user.get_auth_token().decode('utf-8')


class TestUploads:

    @pytest.fixture(scope='function')
    def proc_infra(self, repository_db, worker, no_warn):
        return dict(repository_db=repository_db)

    def assert_uploads(self, upload_json_str, count=0, **kwargs):
        data = json.loads(upload_json_str)
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

    def assert_processing(self, client, test_user_auth, upload_id):
        upload_endpoint = '/uploads/%s' % upload_id

        # poll until completed
        while True:
            time.sleep(0.1)
            rv = client.get(upload_endpoint, headers=test_user_auth)
            assert rv.status_code == 200
            upload = self.assert_upload(rv.data)
            assert 'upload_time' in upload
            if upload['completed']:
                break

        assert len(upload['tasks']) == 4
        assert upload['status'] == 'SUCCESS'
        assert upload['current_task'] == 'cleanup'
        upload_files = UploadFiles.get(upload['upload_id'])
        assert upload_files is not None
        calcs = upload['calcs']['results']
        for calc in calcs:
            assert calc['status'] == 'SUCCESS'
            assert calc['current_task'] == 'archiving'
            assert len(calc['tasks']) == 3
            assert client.get('/archive/logs/%s/%s' % (calc['upload_id'], calc['calc_id']), headers=test_user_auth).status_code == 200

        if upload['calcs']['pagination']['total'] > 1:
            rv = client.get('%s?page=2&per_page=1&order_by=status' % upload_endpoint, headers=test_user_auth)
            assert rv.status_code == 200
            upload = self.assert_upload(rv.data)
            assert len(upload['calcs']['results']) == 1

    def assert_unstage(self, client, test_user_auth, upload_id, proc_infra, meta_data={}):
        rv = client.get('/uploads/%s' % upload_id, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        empty_upload = upload['calcs']['pagination']['total'] == 0

        rv = client.post(
            '/uploads/%s' % upload_id,
            headers=test_user_auth,
            data=json.dumps(dict(operation='unstage', meta_data=meta_data)),
            content_type='application/json')
        assert rv.status_code == 200

        self.assert_upload_does_not_exist(client, upload_id, test_user_auth)
        assert_coe_upload(upload_id, empty=empty_upload, meta_data=meta_data)

    def assert_upload_does_not_exist(self, client, upload_id: str, test_user_auth):
        rv = client.get('/uploads/%s' % upload_id, headers=test_user_auth)
        assert rv.status_code == 404
        assert Upload.objects(upload_id=upload_id).first() is None
        assert Calc.objects(upload_id=upload_id).count() is 0
        upload_files = UploadFiles.get(upload_id)
        assert upload_files is None or isinstance(upload_files, PublicUploadFiles)

    def test_get_command(self, client, test_user_auth, no_warn):
        rv = client.get('/uploads/command', headers=test_user_auth)
        assert rv.status_code == 200
        data = json.loads(rv.data)
        assert 'upload_command' in data
        assert 'upload_url' in data

    def test_get_empty(self, client, test_user_auth, no_warn):
        rv = client.get('/uploads/', headers=test_user_auth)

        assert rv.status_code == 200
        self.assert_uploads(rv.data, count=0)

    def test_get_not_existing(self, client, test_user_auth, no_warn):
        rv = client.get('/uploads/123456789012123456789012', headers=test_user_auth)
        assert rv.status_code == 404

    @pytest.mark.timeout(30)
    @pytest.mark.parametrize('file', example_files)
    @pytest.mark.parametrize('mode', ['multipart', 'stream', 'local_path'])
    @pytest.mark.parametrize('name', [None, 'test_name'])
    def test_put(self, client, test_user_auth, proc_infra, file, mode, name):
        if name:
            url = '/uploads/?name=%s' % name
        else:
            url = '/uploads/'

        if mode == 'multipart':
            rv = client.put(
                url, data=dict(file=(open(file, 'rb'), 'file')), headers=test_user_auth)
        elif mode == 'stream':
            with open(file, 'rb') as f:
                rv = client.put(url, data=f.read(), headers=test_user_auth)
        elif mode == 'local_path':
            url += '&' if name else '?'
            url += 'local_path=%s' % file
            rv = client.put(url, headers=test_user_auth)
        else:
            assert False

        assert rv.status_code == 200
        if mode == 'local_path':
            upload = self.assert_upload(rv.data, local_path=file, name=name)
        else:
            upload = self.assert_upload(rv.data, name=name)

        self.assert_processing(client, test_user_auth, upload['upload_id'])

    def test_delete_not_existing(self, client, test_user_auth, no_warn):
        rv = client.delete('/uploads/123456789012123456789012', headers=test_user_auth)
        assert rv.status_code == 404

    def test_delete_during_processing(self, client, test_user_auth, proc_infra):
        rv = client.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        rv = client.delete('/uploads/%s' % upload['upload_id'], headers=test_user_auth)
        assert rv.status_code == 400
        self.assert_processing(client, test_user_auth, upload['upload_id'])

    def test_delete_unstaged(self, client, test_user_auth, proc_infra, clean_repository_db):
        rv = client.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(client, test_user_auth, upload['upload_id'])
        self.assert_unstage(client, test_user_auth, upload['upload_id'], proc_infra)
        rv = client.delete('/uploads/%s' % upload['upload_id'], headers=test_user_auth)
        assert rv.status_code == 404

    def test_delete(self, client, test_user_auth, proc_infra):
        rv = client.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(client, test_user_auth, upload['upload_id'])
        rv = client.delete('/uploads/%s' % upload['upload_id'], headers=test_user_auth)
        assert rv.status_code == 200
        self.assert_upload_does_not_exist(client, upload['upload_id'], test_user_auth)

    @pytest.mark.parametrize('example_file', example_files)
    def test_post(self, client, test_user_auth, example_file, proc_infra, clean_repository_db):
        rv = client.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(client, test_user_auth, upload['upload_id'])
        self.assert_unstage(client, test_user_auth, upload['upload_id'], proc_infra)

    def test_post_metadata(
            self, client, proc_infra, admin_user_auth, test_user_auth, test_user,
            other_test_user, clean_repository_db):
        rv = client.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(client, test_user_auth, upload['upload_id'])

        meta_data = dict(comment='test comment')
        self.assert_unstage(client, admin_user_auth, upload['upload_id'], proc_infra, meta_data)

    def test_post_metadata_forbidden(self, client, proc_infra, test_user_auth, clean_repository_db):
        rv = client.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
        upload = self.assert_upload(rv.data)
        self.assert_processing(client, test_user_auth, upload['upload_id'])
        rv = client.post(
            '/uploads/%s' % upload['upload_id'],
            headers=test_user_auth,
            data=json.dumps(dict(operation='unstage', meta_data=dict(_pid=256))),
            content_type='application/json')
        assert rv.status_code == 401

    # TODO validate metadata (or all input models in API for that matter)
    # def test_post_bad_metadata(self, client, proc_infra, test_user_auth, clean_repository_db):
    #     rv = client.put('/uploads/?local_path=%s' % example_file, headers=test_user_auth)
    #     upload = self.assert_upload(rv.data)
    #     self.assert_processing(client, test_user_auth, upload['upload_id'])
    #     rv = client.post(
    #         '/uploads/%s' % upload['upload_id'],
    #         headers=test_user_auth,
    #         data=json.dumps(dict(operation='unstage', meta_data=dict(doesnotexist='hi'))),
    #         content_type='application/json')
    #     assert rv.status_code == 400


class UploadFilesBasedTests:

    @staticmethod
    def fix_signature(func, wrapper):
        additional_args = list(inspect.signature(func).parameters.values())[4:]
        wrapper_sig = inspect.signature(wrapper)
        wrapper_args = list(wrapper_sig.parameters.values())[:3] + additional_args
        wrapper_sig = wrapper_sig.replace(parameters=tuple(wrapper_args))
        wrapper.__signature__ = wrapper_sig

    @staticmethod
    def check_authorizaton(func):
        @pytest.mark.parametrize('test_data', [
            [True, None, True],     # in staging for upload
            [True, None, False],    # in staging for different user
            [True, None, None],     # in staging for guest
            [False, True, True],    # in public, restricted for uploader
            [False, True, False],   # in public, restricted for different user
            [False, True, None],    # in public, restricted for guest
            [False, False, True],   # in public, public, for uploader
            [False, False, False],  # in public, public, for different user
            [False, False, None]    # in public, public, for guest
        ], indirect=True)
        def wrapper(self, client, test_data, *args, **kwargs):
            upload, authorized, auth_headers = test_data
            try:
                func(self, client, upload, auth_headers, *args, **kwargs)
            except AssertionError as assertion:
                assertion_str = str(assertion)
                if not authorized:
                    if '0 == 5' in assertion_str and 'ZipFile' in assertion_str:
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
        def wrapper(self, client, test_data, *args, **kwargs):
            upload, _, auth_headers = test_data
            func(self, client, upload, auth_headers, *args, **kwargs)
        UploadFilesBasedTests.fix_signature(func, wrapper)
        return wrapper

    @pytest.fixture(scope='function')
    def test_data(self, request, clean_repository_db, no_warn, test_user, other_test_user):
        # delete potential old test files
        for _ in [0, 1]:
            upload_files = UploadFiles.get('test_upload')
            if upload_files:
                upload_files.delete()

        in_staging, restricted, for_uploader = request.param

        if in_staging:
            authorized = for_uploader
        else:
            authorized = not restricted or for_uploader

        if for_uploader:
            auth_headers = create_auth_headers(test_user)
        elif for_uploader is False:
            auth_headers = create_auth_headers(other_test_user)
        else:
            auth_headers = None

        calc_specs = 'r' if restricted else 'p'
        if in_staging:
            Upload.create(user=test_user, upload_id='test_upload')
            upload_files = create_staging_upload('test_upload', calc_specs=calc_specs)
        else:
            upload_files = create_public_upload('test_upload', calc_specs=calc_specs)
            clean_repository_db.begin()
            coe_upload = coe_repo.Upload(
                upload_name='test_upload',
                user_id=test_user.user_id, is_processed=True)
            clean_repository_db.add(coe_upload)
            clean_repository_db.commit()

        yield 'test_upload', authorized, auth_headers

        upload_files.delete()


class TestArchive(UploadFilesBasedTests):
    @UploadFilesBasedTests.check_authorizaton
    def test_get(self, client, upload, auth_headers):
        rv = client.get('/archive/%s/0' % upload, headers=auth_headers)
        assert rv.status_code == 200
        assert json.loads(rv.data) is not None

    @UploadFilesBasedTests.check_authorizaton
    def test_get_calc_proc_log(self, client, upload, auth_headers):
        rv = client.get('/archive/logs/%s/0' % upload, headers=auth_headers)
        assert rv.status_code == 200
        assert len(rv.data) > 0

    @UploadFilesBasedTests.ignore_authorization
    def test_get_non_existing_archive(self, client, upload, auth_headers):
        rv = client.get('/archive/%s' % 'doesnt/exist', headers=auth_headers)
        assert rv.status_code == 404

    def test_get_metainfo(self, client):
        rv = client.get('/archive/metainfo/all.nomadmetainfo.json')
        assert rv.status_code == 200


class TestRepo(UploadFilesBasedTests):
    @UploadFilesBasedTests.ignore_authorization
    def test_calc(self, client, upload, auth_headers):
        rv = client.get('/repo/%s/0' % upload, headers=auth_headers)
        assert rv.status_code == 200

    @UploadFilesBasedTests.ignore_authorization
    def test_non_existing_calcs(self, client, upload, auth_headers):
        rv = client.get('/repo/doesnt/exist', headers=auth_headers)
        assert rv.status_code == 404

    # def test_calcs(self, client, example_elastic_calc, no_warn):
    #     rv = client.get('/repo/')
    #     assert rv.status_code == 200
    #     data = json.loads(rv.data)
    #     results = data.get('results', None)
    #     assert results is not None
    #     assert isinstance(results, list)
    #     assert len(results) >= 1

    # def test_calcs_pagination(self, client, example_elastic_calc, no_warn):
    #     rv = client.get('/repo/?page=1&per_page=1')
    #     assert rv.status_code == 200
    #     data = json.loads(rv.data)
    #     results = data.get('results', None)
    #     assert results is not None
    #     assert isinstance(results, list)
    #     assert len(results) == 1

    # def test_calcs_user(self, client, example_elastic_calc, test_user_auth, no_warn):
    #     rv = client.get('/repo/?owner=user', headers=test_user_auth)
    #     assert rv.status_code == 200
    #     data = json.loads(rv.data)
    #     results = data.get('results', None)
    #     assert results is not None
    #     assert len(results) >= 1

    # def test_calcs_user_authrequired(self, client, example_elastic_calc, no_warn):
    #     rv = client.get('/repo/?owner=user')
    #     assert rv.status_code == 401

    # def test_calcs_user_invisible(self, client, example_elastic_calc, test_other_user_auth, no_warn):
    #     rv = client.get('/repo/?owner=user', headers=test_other_user_auth)
    #     assert rv.status_code == 200
    #     data = json.loads(rv.data)
    #     results = data.get('results', None)
    #     assert results is not None
    #     assert len(results) == 0


class TestRaw(UploadFilesBasedTests):

    @UploadFilesBasedTests.check_authorizaton
    def test_raw_file(self, client, upload, auth_headers):
        url = '/raw/%s/%s' % (upload, example_file_mainfile)
        rv = client.get(url, headers=auth_headers)
        assert rv.status_code == 200
        assert len(rv.data) > 0

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_missing_file(self, client, upload, auth_headers):
        url = '/raw/%s/does/not/exist' % upload
        rv = client.get(url, headers=auth_headers)
        assert rv.status_code == 404
        data = json.loads(rv.data)
        assert 'files' not in data

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_listing(self, client, upload, auth_headers):
        url = '/raw/%s/examples' % upload
        rv = client.get(url, headers=auth_headers)
        assert rv.status_code == 404
        data = json.loads(rv.data)
        assert len(data['files']) == 5

    @pytest.mark.parametrize('compress', [True, False])
    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_wildcard(self, client, upload, auth_headers, compress):
        url = '/raw/%s/examples*' % upload
        if compress:
            url = '%s?compress=1' % url
        rv = client.get(url, headers=auth_headers)

        assert rv.status_code == 200
        assert len(rv.data) > 0
        with zipfile.ZipFile(io.BytesIO(rv.data)) as zip_file:
            assert zip_file.testzip() is None
            assert len(zip_file.namelist()) == len(example_file_contents)

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_wildcard_missing(self, client, upload, auth_headers):
        url = '/raw/%s/does/not/exist*' % upload
        rv = client.get(url, headers=auth_headers)
        assert rv.status_code == 404

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_file_missing_upload(self, client, upload, auth_headers):
        url = '/raw/doesnotexist/%s' % example_file_mainfile
        rv = client.get(url, headers=auth_headers)
        assert rv.status_code == 404

    @pytest.mark.parametrize('compress', [True, False])
    @UploadFilesBasedTests.check_authorizaton
    def test_raw_files(self, client, upload, auth_headers, compress):
        url = '/raw/%s?files=%s' % (
            upload, ','.join(example_file_contents))
        if compress:
            url = '%s&compress=1' % url
        rv = client.get(url, headers=auth_headers)

        assert rv.status_code == 200
        assert len(rv.data) > 0
        with zipfile.ZipFile(io.BytesIO(rv.data)) as zip_file:
            assert zip_file.testzip() is None
            assert len(zip_file.namelist()) == len(example_file_contents)

    @pytest.mark.parametrize('compress', [True, False, None])
    @UploadFilesBasedTests.check_authorizaton
    def test_raw_files_post(self, client, upload, auth_headers, compress):
        url = '/raw/%s' % upload
        data = dict(files=example_file_contents)
        if compress is not None:
            data.update(compress=compress)
        rv = client.post(url, data=json.dumps(data), content_type='application/json', headers=auth_headers)

        assert rv.status_code == 200
        assert len(rv.data) > 0
        with zipfile.ZipFile(io.BytesIO(rv.data)) as zip_file:
            assert zip_file.testzip() is None
            assert len(zip_file.namelist()) == len(example_file_contents)

    @pytest.mark.parametrize('compress', [True, False])
    @UploadFilesBasedTests.ignore_authorization
    def test_raw_files_missing_file(self, client, upload, auth_headers, compress):
        url = '/raw/%s?files=%s,missing/file.txt' % (upload, example_file_mainfile)
        if compress:
            url = '%s&compress=1' % url
        rv = client.get(url, headers=auth_headers)

        assert rv.status_code == 200
        assert len(rv.data) > 0
        with zipfile.ZipFile(io.BytesIO(rv.data)) as zip_file:
            assert zip_file.testzip() is None
            assert len(zip_file.namelist()) == 1

    @UploadFilesBasedTests.ignore_authorization
    def test_raw_files_missing_upload(self, client, upload, auth_headers):
        url = '/raw/doesnotexist?files=shoud/not/matter.txt'
        rv = client.get(url, headers=auth_headers)

        assert rv.status_code == 404


def test_docs(client):
    rv = client.get('/docs/index.html')
    rv = client.get('/docs/introduction.html')
    assert rv.status_code == 200
