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

from nomad import datamodel, search, processing, parsing
from nomad.search import Entry


def test_init_mapping(elastic):
    pass


def test_index_skeleton_calc(elastic):
    calc_with_metadata = datamodel.CalcWithMetadata(upload_id='test_upload', calc_id='test_calc')

    create_entry(calc_with_metadata)


def test_index_normalized_calc(elastic, normalized: parsing.LocalBackend):
    calc_with_metadata = normalized.to_calc_with_metadata()

    create_entry(calc_with_metadata)


def test_index_normalized_calc_with_metadata(
        elastic, normalized: parsing.LocalBackend, example_user_metadata: dict):

    calc_with_metadata = normalized.to_calc_with_metadata()
    calc_with_metadata.apply_user_metadata(example_user_metadata)

    create_entry(calc_with_metadata)


def test_index_upload(elastic, processed: processing.Upload):
    pass


def create_entry(calc_with_metadata: datamodel.CalcWithMetadata):
    search.Entry.from_calc_with_metadata(calc_with_metadata).save()
    assert_entry(calc_with_metadata.calc_id)


def assert_entry(calc_id):
    calc = Entry.get(calc_id)
    assert calc is not None
