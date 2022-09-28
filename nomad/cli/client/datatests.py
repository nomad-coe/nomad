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

'''
A command that runs a series of data compatibility tests against this version
and another nomad@FAIRDI installation.
'''
from nomad.client import api
from nomad.datamodel.datamodel import EntryArchive
from nomad.metainfo import MEnum


def datatests(auth: api.Auth):
    # All enum values defined in results
    for key, values in get_enums(EntryArchive.results.sub_section, 'results').items():
        test_results_enums(key, values)


def get_enums(root, path='', es=True, results=None):
    '''Recursively navigates down the metainfo hierarchy to extract paths
    and possible values for enum quantities.

    Args:
        root: The root section definition to start from.
        path: Path prefix.
        es: Collect only if stored in ElasticSearch
        results: A dictionary where the results are stored.
    '''
    if results is None:
        results = {}

    subsections = root.sub_sections
    for quantity in root.quantities:
        dtype = quantity.type
        if es and 'elasticsearch' not in quantity.m_annotations:
            continue
        if isinstance(dtype, MEnum):
            values = dtype._values
            quantity_path = f'{path}.{quantity.name}'
            results[quantity_path] = values

    for subsection in subsections:
        section = subsection.sub_section
        get_enums(section, f'{path}.{subsection.name}', es, results)

    return results


def test_results_enums(name, enum_values, auth=None):
    '''Tests that the enum values defined in archive.results cover all of the
    available values in the public part of the data.
    '''
    max_size = 10000
    query = {
        "owner": "public",
        "query": {},
        "aggregations": {
            name: {
                "terms": {
                    "quantity": name,
                    "size": max_size
                },
            },
        },
        "pagination": {
            "page_size": 0,
        }
    }

    response = api.post('entries/query', json=query, auth=auth)
    status_code = response.status_code

    # 422 is allowed if it is about missing doc quantity: it means that this is
    # a new field that does not have any data yet in production, and is thus not
    # mapped in ES yet.
    if status_code == 422:
        error = response.json()['detail']
        assert len(error) == 1
        assert 'is not a doc quantity' in error[0]['msg']
    # If this field is available, test that we are really getting all of the
    # unique values and that the unique terms are a subset of the Enum values
    elif status_code == 200:
        agg_values = {x["value"] for x in response.json()['aggregations'][name]['terms']['data']}
        assert len(agg_values) < max_size
        assert agg_values.issubset(enum_values)
    else:
        raise ValueError(f'Invalid server status code: {status_code}')
