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

from typing import List, Dict, Any, Union

from nomad import config, utils, infrastructure
from nomad.search import refresh, run_on_both_indexes


@run_on_both_indexes
def assert_search_upload(
        entries: Union[int, List] = -1,
        additional_keys: List[str] = [],
        upload_id: str = None,
        index: str = None,
        **kwargs):

    if isinstance(entries, list):
        size = len(entries)
    else:
        size = entries

    keys = ['calc_id', 'upload_id', 'mainfile']
    refresh(index=index)
    body: Dict[str, Any] = {}
    body.update(size=10)
    if upload_id is not None:
        body['query'] = dict(match=dict(upload_id=upload_id))

    search_results = infrastructure.elastic_client.search(index=index, body=body)['hits']

    if size != -1:
        assert search_results['total'] == size

    if search_results['total'] > 0:
        for hit in search_results['hits']:
            hit = utils.flat(hit['_source'])

            for key, value in kwargs.items():
                assert hit.get(key, None) == value

            if 'pid' in hit:
                assert int(hit.get('pid')) > 0

            for key in keys:
                assert key in hit, f'{key} is missing'

            for key in additional_keys:
                assert key in hit, f'{key} is missing'
                assert hit[key] != config.services.unavailable_value

            for coauthor in hit.get('coauthors', []):
                assert coauthor.get('name', None) is not None
