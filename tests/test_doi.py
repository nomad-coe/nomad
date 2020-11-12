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

from nomad import config
from nomad.doi import DOI


def test_create(mongo, test_user, no_warn):
    doi = DOI.create('the_title', test_user)

    assert DOI.objects(doi=doi.doi).first() is not None
    assert doi.metadata_xml is not None


def test_create_doi_counter(mongo, test_user, no_warn):
    DOI.create('the_title', test_user)
    doi = DOI.create('the_title', test_user)
    assert doi.doi.endswith('-2')


def test_create_draft_doi(mongo, test_user, no_warn):
    if config.datacite.enabled:
        doi = DOI.create('the_title', test_user)
        doi.create_draft()
        doi.delete()

        assert DOI.objects(doi=doi.doi).first() is None
