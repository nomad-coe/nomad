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

from nomad.metainfo import Package, MSection

m_package = Package(aliases=['nomad.datamodel.test_package'])


class TestSection(MSection):
    pass


def test_package_registry():
    assert m_package.name == 'tests.metainfo.test_package'
    assert m_package.aliases == ['nomad.datamodel.test_package']
    assert Package.registry['tests.metainfo.test_package'] is m_package
    assert Package.registry['nomad.datamodel.test_package'] is m_package


def test_resolve():
    json_obj = dict(m_def='tests.metainfo.test_package.TestSection')
    metainfo_obj = MSection.from_dict(json_obj)
    assert metainfo_obj.m_def is TestSection.m_def


def test_resolve_with_alias():
    json_obj = dict(m_def='nomad.datamodel.test_package.TestSection')
    metainfo_obj = MSection.from_dict(json_obj)
    assert metainfo_obj.m_def is TestSection.m_def


m_package.__init_metainfo__()
