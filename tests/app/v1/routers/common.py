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
from typing import Set
import re
from devtools import debug
from urllib.parse import urlencode

from nomad.datamodel import results
from nomad.utils import deep_get

from tests.utils import assert_at_least, assert_url_query_args


def post_query_test_parameters(
        entity_id: str, total: int, material_prefix: str, entry_prefix: str):

    elements = f'{material_prefix}elements'
    program_name = f'{entry_prefix}results.method.simulation.program_name'
    method = f'{entry_prefix}results.method'
    properties = f'{entry_prefix}results.properties'
    upload_create_time = f'{entry_prefix}upload_create_time'

    return [
        pytest.param({}, 200, total, id='empty'),
        pytest.param('str', 422, -1, id='not-dict'),
        pytest.param({entity_id: 'id_01'}, 200, 1, id='match'),
        pytest.param({'mispelled': 'id_01'}, 422, -1, id='not-quantity'),
        pytest.param({entity_id: ['id_01', 'id_02']}, 200, 0, id='match-list-0'),
        pytest.param({entity_id: 'id_01', elements: ['H', 'O']}, 200, 1, id='match-list-1'),
        pytest.param({f'{entity_id}:any': ['id_01', 'id_02']}, 200, 2, id='any-short'),
        pytest.param({entity_id: {'any': ['id_01', 'id_02']}}, 200, 2, id='any'),
        pytest.param({entity_id: {'any': 'id_01'}}, 422, -1, id='any-not-list'),
        pytest.param({f'{entity_id}:any': 'id_01'}, 422, -1, id='any-short-not-list'),
        pytest.param({f'{entity_id}:gt': 'id_01'}, 200, total - 1, id='gt-short'),
        pytest.param({entity_id: {'gt': 'id_01'}}, 200, total - 1, id='gt'),
        pytest.param({entity_id: {'gt': ['id_01']}}, 422, total - 1, id='gt-list'),
        pytest.param({entity_id: {'missspelled': 'id_01'}}, 422, -1, id='not-op'),
        pytest.param({f'{entity_id}:lt': ['id_01']}, 422, -1, id='gt-shortlist'),
        pytest.param({f'{entity_id}:misspelled': 'id_01'}, 422, -1, id='not-op-short'),
        pytest.param({'or': [{entity_id: 'id_01'}, {entity_id: 'id_02'}]}, 200, 2, id='or'),
        pytest.param({'or': {entity_id: 'id_01', program_name: 'VASP'}}, 422, -1, id='or-not-list'),
        pytest.param({'and': [{entity_id: 'id_01'}, {entity_id: 'id_02'}]}, 200, 0, id='and'),
        pytest.param({'not': {entity_id: 'id_01'}}, 200, total - 1, id='not'),
        pytest.param({'not': [{entity_id: 'id_01'}]}, 422, -1, id='not-list'),
        pytest.param({'not': {'not': {entity_id: 'id_01'}}}, 200, 1, id='not-nested-not'),
        pytest.param({'not': {f'{entity_id}:any': ['id_01', 'id_02']}}, 200, total - 2, id='not-nested-any'),
        pytest.param({'and': [{f'{entity_id}:any': ['id_01', 'id_02']}, {f'{entity_id}:any': ['id_02', 'id_03']}]}, 200, 1, id='and-nested-any'),
        pytest.param({'and': [{'not': {entity_id: 'id_01'}}, {'not': {entity_id: 'id_02'}}]}, 200, total - 2, id='not-nested-not'),
        pytest.param({method: {'simulation.program_name': 'VASP'}}, 200, total, id='inner-object'),
        pytest.param({f'{properties}.electronic.dos_electronic.band_gap.type': 'direct'}, 200, 1, id='nested-implicit'),
        pytest.param({f'{properties}.electronic.dos_electronic.band_gap': {'type': 'direct'}}, 200, 1, id='nested-explicit'),
        pytest.param({properties: {'electronic.dos_electronic.band_gap': {'type': 'direct'}}}, 200, 1, id='nested-explicit-explicit'),
        pytest.param({f'{upload_create_time}:gt': '1970-01-01'}, 200, total, id='date-1'),
        pytest.param({f'{upload_create_time}:lt': '2099-01-01'}, 200, total, id='date-2'),
        pytest.param({f'{upload_create_time}:gt': '2099-01-01'}, 200, 0, id='date-3')
    ]


def get_query_test_parameters(
        entity_id: str, total: int, material_prefix: str, entry_prefix: str):

    elements = f'{material_prefix}elements'
    n_elements = f'{material_prefix}n_elements'
    upload_create_time = f'{entry_prefix}upload_create_time'

    return [
        pytest.param({}, 200, total, id='empty'),
        pytest.param({entity_id: 'id_01'}, 200, 1, id='match'),
        pytest.param({'mispelled': 'id_01'}, 200, total, id='not-quantity'),
        pytest.param({entity_id: ['id_01', 'id_02']}, 200, 2, id='match-many-or'),
        pytest.param({elements: ['H', 'O']}, 200, total, id='match-list-many-and-1'),
        pytest.param({elements: ['H', 'O', 'Zn']}, 200, 0, id='match-list-many-and-2'),
        pytest.param({n_elements: 2}, 200, total, id='match-int'),
        pytest.param({n_elements + '__gt': 2}, 200, 0, id='gt-int'),
        pytest.param({f'{entity_id}__any': ['id_01', 'id_02']}, 200, 2, id='any'),
        pytest.param({f'{entity_id}__any': 'id_01'}, 200, 1, id='any-not-list'),
        pytest.param({f'{entity_id}__gt': 'id_01'}, 200, total - 1, id='gt'),
        pytest.param({f'{entity_id}__gt': ['id_01', 'id_02']}, 422, -1, id='gt-list'),
        pytest.param({f'{entity_id}__missspelled': 'id_01'}, 422, -1, id='not-op-1'),
        pytest.param({n_elements + '__missspelled': 2}, 422, -1, id='not-op-2'),
        pytest.param({'q': f'{entity_id}__id_01'}, 200, 1, id='q-match'),
        pytest.param({'q': 'missspelled__id_01'}, 422, -1, id='q-bad-quantity'),
        pytest.param({'q': 'bad_encoded'}, 422, -1, id='q-bad-encode'),
        pytest.param({'q': f'{n_elements}__2'}, 200, total, id='q-match-int'),
        pytest.param({'q': f'{n_elements}__gt__2'}, 200, 0, id='q-gt'),
        pytest.param({'q': f'{entry_prefix}upload_create_time__gt__2014-01-01'}, 200, total, id='datetime'),
        pytest.param({'q': [elements + '__all__H', elements + '__all__O']}, 200, total, id='q-all'),
        pytest.param({'q': [elements + '__all__H', elements + '__all__X']}, 200, 0, id='q-all'),
        pytest.param({'q': f'{upload_create_time}__gt__1970-01-01'}, 200, total, id='date'),
        pytest.param({'json_query': f'{{"{elements}": ["H", "O"]}}'}, 200, total, id='json_query'),
        pytest.param({'json_query': f'{{"{elements}": ["H", "O"}}'}, 422, 0, id='invalid-json_query')
    ]


def owner_test_parameters():
    return [
        # (owner, user, status_code, total_entries, total_mainfiles, total_materials)
        pytest.param('user', None, 401, -1, -1, -1, id='user-wo-auth'),
        pytest.param('staging', None, 401, -1, -1, -1, id='staging-wo-auth'),
        pytest.param('visible', None, 200, 23, 23, 6, id='visible-wo-auth'),
        pytest.param('admin', None, 401, -1, -1, -1, id='admin-wo-auth'),
        pytest.param('shared', None, 401, -1, -1, -1, id='shared-wo-auth'),
        pytest.param('public', None, 200, 23, 23, 6, id='public-wo-auth'),

        pytest.param('user', 'test_user', 200, 32, 30, 13, id='user-test-user'),
        pytest.param('staging', 'test_user', 200, 6, 4, 4, id='staging-test-user'),
        pytest.param('visible', 'test_user', 200, 32, 30, 13, id='visible-test-user'),
        pytest.param('admin', 'test_user', 401, -1, -1, -1, id='admin-test-user'),
        pytest.param('shared', 'test_user', 200, 32, 30, 13, id='shared-test-user'),
        pytest.param('public', 'test_user', 200, 23, 23, 6, id='public-test-user'),

        pytest.param('user', 'other_test_user', 200, 0, 0, 0, id='user-other-test-user'),
        pytest.param('staging', 'other_test_user', 200, 2, 2, 2, id='staging-other-test-user'),
        pytest.param('visible', 'other_test_user', 200, 27, 27, 10, id='visible-other-test-user'),
        pytest.param('shared', 'other_test_user', 200, 4, 4, 4, id='shared-other-test-user'),
        pytest.param('public', 'other_test_user', 200, 23, 23, 6, id='public-other-test-user'),

        pytest.param('all', None, 200, 26, 26, 9, id='metadata-all-wo-auth'),
        pytest.param('all', 'test_user', 200, 32, 30, 13, id='metadata-all-test-user'),
        pytest.param('all', 'other_test_user', 200, 28, 28, 11, id='metadata-all-other-test-user'),

        pytest.param('admin', 'admin_user', 200, 32, 30, 13, id='admin-admin-user'),
        pytest.param('all', 'bad_user', 401, -1, -1, -1, id='bad-credentials')
    ]


def pagination_test_parameters(elements: str, n_elements: str, crystal_system: str, total: int):
    return [
        pytest.param({}, {'total': total, 'page_size': 10, 'next_page_after_value': 'id_10'}, 200, id='empty'),
        pytest.param({'page_size': 1}, {'total': total, 'page_size': 1, 'next_page_after_value': 'id_01'}, 200, id='size'),
        pytest.param({'page_size': 0}, {'total': total, 'page_size': 0}, 200, id='size-0'),
        pytest.param({'page_size': 1, 'page_after_value': 'id_01'}, {'page_after_value': 'id_01', 'next_page_after_value': 'id_02'}, 200, id='after'),
        pytest.param({'page_size': 1, 'page_after_value': 'id_02', 'order': 'desc'}, {'next_page_after_value': 'id_01'}, 200, id='after-desc'),
        pytest.param({'page_size': 10, 'page_after_value': 'id_22', 'order': 'asc'}, {'next_page_after_value': None}, 200, id='after-exhausted'),
        pytest.param({'page_size': 1, 'order_by': n_elements}, {'next_page_after_value': '2:id_01'}, 200, id='order-by-after-int'),
        pytest.param({'page_size': 1, 'order_by': crystal_system}, {'next_page_after_value': 'cubic:id_01'}, 200, id='order-by-after-nested'),
        pytest.param({'page_size': -1}, None, 422, id='bad-size'),
        pytest.param({'order': 'misspelled'}, None, 422, id='bad-order'),
        pytest.param({'order_by': 'misspelled'}, None, 422, id='bad-order-by'),
        pytest.param({'order_by': elements, 'page_after_value': 'H:id_01'}, None, 422, id='order-by-list'),
        pytest.param({'order_by': n_elements, 'page_after_value': 'some'}, None, 400, id='order-by-bad-after'),
        pytest.param({'page_offset': 0, 'page_size': 1}, {'total': total, 'next_page_after_value': 'id_01', 'page_offset': 0}, 200, id='page-offset-1'),
        pytest.param({'page_offset': 1, 'page_size': 1}, {'total': total, 'next_page_after_value': 'id_02', 'page_offset': 1}, 200, id='page-offset-2'),
        pytest.param({'page_offset': 9999}, None, 422, id='page-offset-too-large'),
        pytest.param({'page_offset': 9989}, None, 200, id='page-offset-just-small-enough'),
        pytest.param({'page': 1, 'page_size': 1}, {'total': total, 'page_size': 1, 'next_page_after_value': 'id_01', 'page': 1}, 200, id='page-1'),
        pytest.param({'page': 2, 'page_size': 1}, {'total': total, 'page_size': 1, 'next_page_after_value': 'id_02', 'page': 2}, 200, id='page-2'),
        pytest.param({'page': 1000, 'page_size': 10}, None, 422, id='page-too-large'),
        pytest.param({'page': 9999, 'page_size': 1}, None, 200, id='page-just-small-enough'),
        pytest.param({'page_offset': 1, 'page': 1}, None, 422, id='only-one-page-param'),
        pytest.param({'page_offset': 1, 'page_size': 0}, None, 422, id='page-param-only-with-page-size')
    ]


def aggregation_test_parameters(entity_id: str, material_prefix: str, entry_prefix: str, total: int):
    n_code_names = results.Simulation.program_name.a_elasticsearch[0].default_aggregation_size
    program_name = f'{entry_prefix}results.method.simulation.program_name'
    n_calculations = f'{entry_prefix}results.properties.n_calculations'
    upload_create_time = f'{entry_prefix}upload_create_time'

    return [
        pytest.param(
            {'statistics': {
                'metrics': ['n_entries', 'n_materials', 'n_uploads', 'n_calculations']
            }},
            3, 3, 200, 'test_user', id='statistics'
        ),
        pytest.param(
            {'terms': {'quantity': f'{entry_prefix}upload_id'}},
            8, 8, 200, 'test_user', id='default'),
        pytest.param(
            {
                'terms': {
                    'quantity': f'{entry_prefix}upload_id',
                    'pagination': {'order_by': f'{entry_prefix}main_author.user_id'}
                }
            },
            8, 8, 200, 'test_user', id='order-str'),
        pytest.param(
            {
                'terms': {
                    'quantity': f'{entry_prefix}upload_id',
                    'pagination': {'order_by': upload_create_time}
                }
            },
            8, 8, 200, 'test_user', id='order-date'),
        pytest.param(
            {
                'terms': {
                    'quantity': f'{entry_prefix}upload_id',
                    'pagination': {'order_by': f'{entry_prefix}results.properties.n_calculations'}
                }
            },
            8, 8, 200, 'test_user', id='order-int'),
        pytest.param(
            {'terms': {'quantity': f'{material_prefix}symmetry.structure_name'}},
            0, 0, 200, 'test_user', id='no-results'),
        pytest.param(
            {
                'terms': {
                    'quantity': f'{entry_prefix}upload_id',
                    'pagination': {'page_after_value': 'id_published'}
                }
            },
            8, 4, 200, 'test_user', id='after'),
        pytest.param(
            {
                'terms': {
                    'quantity': f'{entry_prefix}upload_id',
                    'pagination': {
                        'order_by': f'{entry_prefix}main_author.name',
                        'page_after_value': 'Sheldon Cooper:id_published'
                    }
                }
            },
            8, 4, 200, 'test_user', id='after-order'),
        pytest.param(
            {'terms': {'quantity': f'{entry_prefix}upload_id', 'entries': {'size': 10}}},
            8, 8, 200, 'test_user', id='entries'),
        pytest.param(
            {'terms': {'quantity': f'{entry_prefix}upload_id', 'entries': {'size': 1}}},
            8, 8, 200, 'test_user', id='entries-size'),
        pytest.param(
            {'terms': {'quantity': f'{entry_prefix}upload_id', 'entries': {'size': 0}}},
            -1, -1, 422, 'test_user', id='bad-entries'),
        pytest.param(
            {
                'terms': {
                    'quantity': f'{entry_prefix}upload_id',
                    'entries': {
                        'size': 10,
                        'required': {
                            'include': [f'{entry_prefix}entry_id', f'{entry_prefix}main_author.*']
                        }
                    }
                }
            },
            8, 8, 200, 'test_user', id='entries-include'),
        pytest.param(
            {'terms': {'quantity': program_name}},
            n_code_names, n_code_names, 200, None, id='fixed-values'),
        pytest.param(
            {'terms': {'quantity': program_name, 'metrics': ['n_uploads']}},
            n_code_names, n_code_names, 200, None, id='metrics'),
        pytest.param(
            {'terms': {'quantity': program_name, 'metrics': ['does not exist']}},
            -1, -1, 422, None, id='bad-metric'),
        pytest.param(
            {'terms': {'quantity': entity_id, 'size': 1000}},
            total, total, 200, None, id='size-to-large'),
        pytest.param(
            {'terms': {'quantity': entity_id, 'size': 5}},
            total, 5, 200, None, id='size'),
        pytest.param(
            {'terms': {'quantity': entity_id, 'size': -1}},
            -1, -1, 422, None, id='bad-size-1'),
        pytest.param(
            {'terms': {'quantity': entity_id, 'size': 0}},
            -1, -1, 422, None, id='bad-size-2'),
        pytest.param(
            {'terms': {'quantity': entity_id}},
            total, 10 if total > 10 else total, 200, None, id='size-default'),
        pytest.param(
            {
                'terms': {
                    'quantity': f'{entry_prefix}upload_id',
                    'pagination': {'order': 'asc'}
                }
            },
            8, 8, 200, 'test_user', id='order-direction'),
        pytest.param(
            {'terms': {'quantity': 'does not exist'}},
            -1, -1, 422, None, id='bad-quantity'),
        pytest.param(
            {'date_histogram': {'quantity': upload_create_time}},
            1, 1, 200, 'test-user', id='date-histogram'
        ),
        pytest.param(
            {'date_histogram': {'quantity': upload_create_time, 'metrics': ['n_uploads']}},
            1, 1, 200, 'test-user', id='date-histogram-metrics'
        ),
        pytest.param(
            {'date_histogram': {'quantity': upload_create_time, 'interval': '1s'}},
            1, 1, 200, 'test-user', id='date-histogram-interval'
        ),
        pytest.param(
            {'date_histogram': {'quantity': 'upload_id'}},
            -1, -1, 422, 'test-user', id='date-histogram-no-date'
        ),
        pytest.param(
            {'date_histogram': {'quantity': 'upload_id', 'interval': '1xy'}},
            -1, -1, 422, 'test-user', id='date-histogram-bad-interval'
        ),
        pytest.param(
            {'histogram': {'quantity': n_calculations, 'interval': 1}},
            1, 1, 200, None, id='histogram'
        ),
        pytest.param(
            {'histogram': {'quantity': n_calculations, 'interval': 1, 'metrics': ['n_uploads']}},
            1, 1, 200, None, id='histogram-metric'
        ),
        pytest.param(
            {'histogram': {'quantity': n_calculations}},
            -1, -1, 422, None, id='histogram-no-interval'
        ),
        pytest.param(
            {'histogram': {'quantity': 'upload_id'}},
            -1, -1, 422, None, id='histogram-no-number'
        ),
        pytest.param(
            {'min_max': {'quantity': n_calculations}},
            1, 1, 200, None, id='min-max'
        ),
        pytest.param(
            {'min_max': {'quantity': 'upload_id'}},
            -1, -1, 422, None, id='min-max-no-number'
        ),
    ]


def aggregation_exclude_from_search_test_parameters(resource: str, total_per_entity: int, total: int):
    if resource == "materials":
        entry_prefix = "entries."
        material_prefix = ""
    elif resource == "entries":
        entry_prefix = ""
        material_prefix = "results.material."
    else:
        raise ValueError("invalid resource")

    entry_id = f'{entry_prefix}entry_id'
    upload_id = f'{entry_prefix}upload_id'
    program_name = f'{entry_prefix}results.method.simulation.program_name'
    n_elements = f'{material_prefix}n_elements'
    band_gap = f'{entry_prefix}results.properties.electronic.band_structure_electronic.band_gap.value'

    def make_aggs(aggs):
        """Given a list of aggregation definitions, returns the API-compatible
        aggregations together with their types and expected result sizes.
        """
        aggs_api = {}
        types = []
        lengths = []
        for i, agg in enumerate(aggs):
            outer_agg = {f'{agg["type"]}': agg}
            aggs_api[f'agg_{i}'] = outer_agg
            types.append(agg["type"])
            lengths.append(agg["n_results"])
            del agg["n_results"]
            del agg["type"]
        return [aggs_api, types, lengths]

    return [
        pytest.param(
            {
                f'{entry_id}:any': ['id_01'],
                upload_id: 'id_published',
                program_name: 'VASP'
            },
            make_aggs([]), 1, 200,
            id='empty'
        ),
        pytest.param(
            {
                f'{entry_id}:any': ['id_01']
            },
            make_aggs([
                {
                    'type': 'terms',
                    'quantity': entry_id,
                    'exclude_from_search': True,
                    'n_results': 10,
                }
            ]),
            1, 200,
            id='exclude'
        ),
        pytest.param(
            {
                f'{entry_id}:any': ['id_01']
            },
            make_aggs([
                {
                    'type': 'terms',
                    'quantity': entry_id,
                    'exclude_from_search': False,
                    'n_results': total_per_entity,
                }
            ]),
            1, 200,
            id='dont-exclude'
        ),
        pytest.param(
            {
                f'{entry_id}:any': ['id_01'],
                upload_id: 'id_published',
                program_name: 'VASP'
            },
            make_aggs([
                {
                    'type': 'terms',
                    'quantity': entry_id,
                    'exclude_from_search': True,
                    'n_results': 10,
                },
                {
                    'type': 'terms',
                    'quantity': upload_id,
                    'exclude_from_search': True,
                    'n_results': 1,
                }
            ]),
            1, 200,
            id='two-aggs'
        ),
        pytest.param(
            {
                f'{entry_id}:any': ['id_01']
            },
            make_aggs([
                {
                    'type': 'terms',
                    'quantity': entry_id,
                    'exclude_from_search': True,
                    'n_results': 10,
                },
                {
                    'type': 'terms',
                    'quantity': entry_id,
                    'exclude_from_search': False,
                    'n_results': total_per_entity,
                }
            ]),
            1, 200,
            id='two-aggs-same-quantity'
        ),
        pytest.param(
            {},
            make_aggs([
                {
                    'type': 'terms',
                    'quantity': entry_id,
                    'exclude_from_search': True,
                    'n_results': 10,
                }
            ]),
            total, 200,
            id='not-in-query'
        ),
        pytest.param(
            {},
            make_aggs([
                {
                    'type': 'terms',
                    'quantity': entry_id,
                    'exclude_from_search': True,
                    'pagination': {
                        'page_size': 20
                    },
                    'n_results': 20,
                }
            ]),
            total, 422,
            id='with-pagination'
        ),
        pytest.param(
            {},
            make_aggs([
                {
                    'type': 'terms',
                    'quantity': entry_id,
                    'size': 20,
                    'exclude_from_search': True,
                    'n_results': 20,
                }
            ]),
            total, 200,
            id='with-size'
        ),
        pytest.param(
            {
                'or': [{entry_id: 'id_01'}, {entry_id: 'id_05'}]
            },
            make_aggs([
                {
                    'type': 'terms',
                    'quantity': entry_id,
                    'exclude_from_search': True,
                    'n_results': 10,
                }
            ]),
            2, 200,
            id='non-dict-query'
        ),
        pytest.param(
            {
                f'{n_elements}': {'gte': 0, 'lte': 10}
            },
            make_aggs([
                {
                    'type': 'min_max',
                    'quantity': n_elements,
                    'exclude_from_search': True,
                    'n_results': 2,
                }
            ]),
            total, 200,
            id='range_min_max'
        ),
        pytest.param(
            {
                f'{band_gap}': {'gte': 0, 'lte': 10}
            },
            make_aggs([
                {
                    'type': 'min_max',
                    'quantity': band_gap,
                    'exclude_from_search': True,
                    'n_results': 2,
                }
            ]),
            0, 200,
            id='nested_range_min_max'
        )
    ]


def assert_response(response, status_code=None):
    ''' General assertions for status_code and error messages '''
    if status_code and response.status_code != status_code:
        try:
            debug(response.json())
        except Exception:
            pass

    if status_code is not None:
        if response.status_code != status_code and response.status_code == 422:
            print(response.json()['detail'])
        assert response.status_code == status_code

    if status_code == 422:
        response_json = response.json()
        details = response_json['detail']
        assert len(details) > 0
        for detail in details:
            assert 'loc' in detail
            assert 'msg' in detail
        return

    if 400 <= status_code < 500:
        response_json = response.json()
        assert 'detail' in response_json


def assert_base_metadata_response(response, status_code=None):
    assert_response(response, status_code)

    if status_code != 200 or response.status_code != 200:
        return None

    response_json = response.json()
    assert 'es_query' not in response_json
    assert 'data' in response_json
    return response_json


def assert_metadata(response_json):
    if isinstance(response_json['data'], list):
        metadatas = response_json['data']
    else:
        metadatas = [response_json['data']]

    for metadata in metadatas:
        if 'required' not in response_json:
            assert 'license' in metadata

        if 'main_author' in metadata:
            assert 'email' not in metadata['main_author']


def assert_metadata_response(response, status_code=None):
    response_json = assert_base_metadata_response(response, status_code=status_code)
    if response_json is not None:
        assert_metadata(response_json)
    return response_json


def assert_required(data, required, default_key: str):
    # We flat out all keys in data and then make sure that the full qualified keys in the
    # data are consistent with the keys given in the required include and exclude.
    keys: Set[str] = set()

    def collect_keys(data, prefix=None):
        if isinstance(data, list):
            for item in data:
                collect_keys(item, prefix=prefix)

        elif isinstance(data, dict):
            for key, value in data.items():
                collect_keys(value, prefix=f'{prefix}.{key}' if prefix is not None else key)

        else:
            keys.add(prefix)

    collect_keys(data)

    if 'include' in required:
        for key in keys:
            found_include = False
            for include in required['include']:
                include_re = include.replace('.', r'\.').replace('*', r'[^\.\*]*')
                if key == default_key or re.match(include_re, key):
                    found_include = True

            assert found_include, key
    if 'exclude' in required:
        for exclude in required['exclude']:
            found_exclude = None
            for key in keys:
                exclude_re = exclude.replace('.', r'\.').replace('*', r'[^\.\*]*')
                if key != default_key and re.match(exclude_re, key):
                    found_exclude = key

            assert found_exclude is None, f'{exclude} excluded but found {found_exclude}'


def assert_aggregations(
        response_json, name, agg,
        total: int = -1, size: int = -1, default_key: str = None):
    assert 'aggregations' in response_json
    assert name in response_json['aggregations']
    agg_response_obj = response_json['aggregations'][name]
    assert len(agg_response_obj) == 1
    agg_type = next(iter(agg_response_obj.keys()))
    agg_response = agg_response_obj[agg_type]

    assert 'data' in agg_response
    if agg_type != 'statistics':
        assert 'quantity' in agg_response

    assert_at_least(agg, agg_response)

    data = agg_response['data']
    n_data = len(data)

    if 'pagination' in agg:
        assert agg_response['pagination']['total'] >= n_data
        if size >= 0:
            assert n_data == size
        if total >= 0:
            assert agg_response['pagination']['total'] == total

        assert_pagination(agg.get('pagination', {}), agg_response['pagination'], data, is_get=False)

    if agg_type == 'min_max':
        assert len(data) == 2
        assert isinstance(data[0], (float, int))
        assert isinstance(data[1], (float, int))
    elif agg_type == 'statistics':
        assert 'metrics' in agg_response
        for metric in agg.get('metrics', []):
            assert metric in data
            assert isinstance(data[metric], (float, int))
    else:
        assert total == -1 or total >= n_data
        assert size == -1 or size == n_data

        for bucket in data:
            assert 'value' in bucket
            if len(agg.get('metrics', [])) > 0:
                assert 'metrics' in bucket
            else:
                assert 'metrics' not in bucket
            assert 'count' in bucket

            value = bucket['value']
            if agg_type == 'date_histogram': assert re.match(r'\d{4}\-\d{2}\-\d{2}', value)
            elif agg_type == 'histogram': assert isinstance(value, (float, int))

            for metric in agg.get('metrics', []):
                assert metric in bucket['metrics']
                assert isinstance(bucket['metrics'][metric], (float, int))

    if 'entries' in agg:
        for bucket in data:
            assert 'entries' in bucket
            assert agg['entries'].get('size', 10) >= len(bucket['entries']) > 0
            if 'required' in agg['entries']:
                for entry in bucket['entries']:
                    assert_required(entry, agg['entries']['required'], default_key=default_key)


def assert_pagination(pagination, pagination_response, data, order_by=None, order=None, is_get=True):
    assert_at_least(pagination, pagination_response)
    assert len(data) <= pagination_response['page_size']
    assert len(data) <= pagination_response['total']

    if order is None:
        order = pagination_response.get('order', 'asc')
    if order_by is None:
        order_by = pagination_response.get('order_by')

    if order_by is not None:
        for index, item in enumerate(data):
            if index < len(data) - 1 and order_by in item:
                if order == 'desc':
                    assert item[order_by] >= data[index + 1][order_by]
                else:
                    assert item[order_by] <= data[index + 1][order_by]

    if is_get:
        page_size = pagination_response['page_size']
        page = pagination_response.get('page')
        page_url = pagination_response.get('page_url')
        first_page_url = pagination_response.get('first_page_url')
        prev_page_url = pagination_response.get('prev_page_url')
        next_page_url = pagination_response.get('next_page_url')
        next_page_after_value = pagination_response.get('next_page_after_value')

        assert page_url
        if page_size:
            assert first_page_url
            assert_url_query_args(first_page_url, page_after_value=None, page=None)
        if next_page_after_value:
            assert next_page_url
            assert_url_query_args(next_page_url, page_after_value=next_page_after_value, page=None)
        if page and page > 1:
            assert prev_page_url
            assert_url_query_args(prev_page_url, page=page - 1, page_after_value=None)


def perform_metadata_test(
        client, endpoint: str, owner=None, headers={}, status_code=200,
        total=None, http_method='get', **kwargs):

    if http_method == 'get':
        params = {}
        if owner is not None:
            params['owner'] = owner
        for value in kwargs.values():
            params.update(**value)
        response = client.get(
            f'{endpoint}?{urlencode(params, doseq=True)}', headers=headers)

    elif http_method == 'post':
        body = dict(**kwargs)
        if owner is not None:
            body['owner'] = owner
        response = client.post(f'{endpoint}/query', headers=headers, json=body)

    else:
        assert False

    response_json = assert_metadata_response(response, status_code=status_code)

    if response_json is None:
        return

    assert 'pagination' in response_json

    if total is not None and total >= 0:
        assert response_json['pagination']['total'] == total, response_json['pagination']['total']

    return response_json


def perform_owner_test(
        client, test_user_auth, other_test_user_auth, admin_user_auth,
        owner, user, status_code, total, http_method, test_method):

    headers = None
    if user == 'test_user':
        headers = test_user_auth
    elif user == 'other_test_user':
        headers = other_test_user_auth
    elif user == 'admin_user':
        headers = admin_user_auth
    elif user == 'bad_user':
        headers = {'Authorization': 'Bearer NOTATOKEN'}

    test_method(
        client, headers=headers, owner=owner, status_code=status_code, total=total,
        http_method=http_method)


def perform_quantity_search_test(quantity, resource, search, result, client):
    if resource == "materials":
        if quantity.startswith("results.material"):
            quantity = quantity[len("results.material."):]
        else:
            quantity = f"entries.{quantity}"

    body = {"query": {f"{quantity}": search}}
    response = client.post(f"{resource}/query", json=body, headers={})
    assert_response(response, 200)
    response = response.json()
    api_result = deep_get(response["data"][0], *quantity.split("."))
    assert api_result == result
