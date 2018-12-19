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

import os
import os.path
import shutil
import pytest

from nomad import config
from nomad.uploads import Metadata, MetadataTimeout


class TestMetadata:
    @pytest.fixture(scope='function')
    def test_dir(self):
        path = os.path.join(config.fs.tmp, 'test_dir')
        os.makedirs(path)
        yield path
        shutil.rmtree(path)

    @pytest.fixture(scope='function')
    def md(self, test_dir):
        with Metadata(test_dir) as md:
            yield md

    def test_open_empty(self, test_dir):
        with Metadata(test_dir):
            pass

    def test_modify(self, test_dir):
        with Metadata(test_dir) as md:
            md.data['key'] = 'value'

        with Metadata(test_dir) as md:
            assert 'key' in md.data
            assert md.data['key'] == 'value'
            assert len(md.data) == 1

    def test_lock(self, test_dir):
        timeout = False
        with Metadata(test_dir):
            try:
                with Metadata(test_dir, lock_timeout=0.1):
                    pass
            except MetadataTimeout:
                timeout = True
        assert timeout

    def test_insert(self, md: Metadata):
        md.insert(dict(hash='0', data='test'))
        assert len(md.data) == 1
        assert '0' in md.data
        assert md.data['0']['data'] == 'test'

    def test_insert_fail(self, md: Metadata):
        failed = False
        md.insert(dict(hash='0', data='test'))
        try:
            md.insert(dict(hash='0', data='test'))
        except Exception:
            failed = True

        assert failed
        assert len(md.data) == 1

    def test_update(self, md: Metadata):
        md.insert(dict(hash='0', data='test'))
        md.update(dict(hash='0', data='updated'))
        assert len(md.data) == 1
        assert md.data['0']['data'] == 'updated'

    def test_update_fail(self, md: Metadata):
        failed = False
        try:
            md.update(dict(hash='0', data='updated'))
        except KeyError:
            failed = True
        assert failed
        assert len(md.data) == 0
