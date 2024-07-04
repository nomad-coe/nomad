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

import numpy as np
import pytest

from nomad.graph.graph_reader import DefinitionReader
from nomad.metainfo import (
    Package,
    MSection,
    Quantity,
    Reference,
    SectionProxy,
    SubSection,
)

m_package = Package()


class Inner(MSection):
    n_impurities = Quantity(type=np.int32)


class Base(MSection):
    dimensionality = Quantity(type=np.int32)

    n_points = Quantity(type=np.int32)


class Derived(Base):
    weights = Quantity(type=np.float64, shape=['*'])

    inner = Quantity(type=Reference(SectionProxy('Inner')))


class Holder(MSection):
    base = Quantity(type=Reference(SectionProxy('Base')), shape=[])

    derived = Quantity(type=Reference(SectionProxy('Derived')), shape=[])

    derived_section = SubSection(sub_section=Derived.m_def)


m_package.init_metainfo()

m_def = Holder.m_def

prefix = 'metainfo/tests.graph.test_definition_reader/section_definitions'


def remove_cache(result):
    if '__CACHE__' in result:
        del result['__CACHE__']
    return result


def assert_list(l1, l2):
    assert len(l1) == len(l2)
    for i, j in zip(l1, l2):
        if isinstance(i, dict):
            assert_dict(i, j)
        elif isinstance(i, list):
            assert_list(i, j)
        else:
            assert i == j


def assert_dict(d1, d2):
    assert set(d1.keys()) == set(d2.keys())
    for k, v in d1.items():
        if isinstance(v, dict):
            assert_dict(v, d2[k])
        elif isinstance(v, list):
            assert_list(v, d2[k])
        else:
            assert v == d2[k]


@pytest.mark.parametrize(
    'query,result',
    [
        # plain get default quantities
        # the references are not resolved
        pytest.param(
            {'m_request': {'directive': 'plain'}},
            {
                'm_def': f'{prefix}/3',
                'metainfo': {
                    'tests.graph.test_definition_reader': {
                        'section_definitions': [
                            None,
                            None,
                            None,
                            {
                                'name': 'Holder',
                                'quantities': [
                                    {
                                        'name': 'base',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/1',
                                        },
                                        'shape': [],
                                    },
                                    {
                                        'name': 'derived',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/2',
                                        },
                                        'shape': [],
                                    },
                                ],
                                'sub_sections': [
                                    {
                                        'name': 'derived_section',
                                        'sub_section': f'{prefix}/2',
                                    }
                                ],
                            },
                        ]
                    }
                },
            },
            id='plain-retrieval',
        ),
        # now resolve all referenced quantities and sections
        pytest.param(
            {'m_request': {'directive': 'resolved'}},
            {
                'm_def': f'{prefix}/3',
                'metainfo': {
                    'tests.graph.test_definition_reader': {
                        'section_definitions': [
                            {
                                'name': 'Inner',
                                'quantities': [
                                    {
                                        'name': 'n_impurities',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    }
                                ],
                            },
                            {
                                'name': 'Base',
                                'quantities': [
                                    {
                                        'name': 'dimensionality',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    },
                                    {
                                        'name': 'n_points',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    },
                                ],
                            },
                            {
                                'name': 'Derived',
                                'base_sections': [f'{prefix}/1'],
                                'quantities': [
                                    {
                                        'name': 'weights',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'float64',
                                        },
                                        'shape': ['*'],
                                    },
                                    {
                                        'name': 'inner',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/0',
                                        },
                                    },
                                ],
                            },
                            {
                                'name': 'Holder',
                                'quantities': [
                                    {
                                        'name': 'base',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/1',
                                        },
                                        'shape': [],
                                    },
                                    {
                                        'name': 'derived',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/2',
                                        },
                                        'shape': [],
                                    },
                                ],
                                'sub_sections': [f'{prefix}/2'],
                            },
                        ]
                    }
                },
            },
            id='resolve-all',
        ),
        # limit resolve depth via depth
        # as a result, inner is not resolved
        pytest.param(
            {'m_request': {'directive': 'resolved', 'depth': 1}},
            {
                'm_def': f'{prefix}/3',
                'metainfo': {
                    'tests.graph.test_definition_reader': {
                        'section_definitions': [
                            None,
                            {
                                'name': 'Base',
                                'quantities': [
                                    {
                                        'name': 'dimensionality',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    },
                                    {
                                        'name': 'n_points',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    },
                                ],
                            },
                            {
                                'name': 'Derived',
                                'base_sections': [f'{prefix}/1'],
                                'quantities': [
                                    {
                                        'name': 'weights',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'float64',
                                        },
                                        'shape': ['*'],
                                    },
                                    {
                                        'name': 'inner',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/0',
                                        },
                                    },
                                ],
                            },
                            {
                                'name': 'Holder',
                                'quantities': [
                                    {
                                        'name': 'base',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/1',
                                        },
                                        'shape': [],
                                    },
                                    {
                                        'name': 'derived',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/2',
                                        },
                                        'shape': [],
                                    },
                                ],
                                'sub_sections': [f'{prefix}/2'],
                            },
                        ]
                    }
                },
            },
            id='resolve-with-depth',
        ),
        pytest.param(
            {
                'all_quantities': {'m_request': {'directive': 'plain'}},
                'inherited_sections': {'m_request': {'directive': 'plain'}},
                'all_base_sections': {'m_request': {'directive': 'plain'}},
            },
            {
                'm_def': f'{prefix}/3',
                'metainfo': {
                    'tests.graph.test_definition_reader': {
                        'section_definitions': [
                            None,
                            None,
                            None,
                            {
                                'all_quantities': {
                                    'base': f'{prefix}/3/quantities/0',
                                    'derived': f'{prefix}/3/quantities/1',
                                },
                                'inherited_sections': [
                                    'metainfo/tests.graph.test_definition_reader/section_definitions/3'
                                ],
                                'all_base_sections': [],
                            },
                        ]
                    }
                },
            },
            id='get-derived-only',
        ),
        # resolve all quantities
        pytest.param(
            {'all_quantities': {'m_request': {'directive': 'resolved'}}},
            {
                'm_def': f'{prefix}/3',
                'metainfo': {
                    'tests.graph.test_definition_reader': {
                        'section_definitions': [
                            {
                                'name': 'Inner',
                                'quantities': [
                                    {
                                        'name': 'n_impurities',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    }
                                ],
                            },
                            {
                                'name': 'Base',
                                'quantities': [
                                    {
                                        'name': 'dimensionality',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    },
                                    {
                                        'name': 'n_points',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    },
                                ],
                            },
                            {
                                'name': 'Derived',
                                'base_sections': [f'{prefix}/1'],
                                'quantities': [
                                    {
                                        'name': 'weights',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'float64',
                                        },
                                        'shape': ['*'],
                                    },
                                    {
                                        'name': 'inner',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/0',
                                        },
                                    },
                                ],
                            },
                            {
                                'all_quantities': {
                                    'base': f'{prefix}/3/quantities/0',
                                    'derived': f'{prefix}/3/quantities/1',
                                },
                                'quantities': [
                                    {
                                        'name': 'base',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/1',
                                        },
                                        'shape': [],
                                    },
                                    {
                                        'name': 'derived',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/2',
                                        },
                                        'shape': [],
                                    },
                                ],
                            },
                        ]
                    }
                },
            },
            id='get-derived-resolved',
        ),
        pytest.param(
            {'all_quantities': {'m_request': {'directive': 'resolved', 'depth': 1}}},
            {
                'm_def': f'{prefix}/3',
                'metainfo': {
                    'tests.graph.test_definition_reader': {
                        'section_definitions': [
                            None,
                            {
                                'name': 'Base',
                                'quantities': [
                                    {
                                        'name': 'dimensionality',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    },
                                    {
                                        'name': 'n_points',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'int32',
                                        },
                                    },
                                ],
                            },
                            {
                                'name': 'Derived',
                                'base_sections': [f'{prefix}/1'],
                                'quantities': [
                                    {
                                        'name': 'weights',
                                        'type': {
                                            'type_kind': 'numpy',
                                            'type_data': 'float64',
                                        },
                                        'shape': ['*'],
                                    },
                                    {
                                        'name': 'inner',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/0',
                                        },
                                    },
                                ],
                            },
                            {
                                'all_quantities': {
                                    'base': f'{prefix}/3/quantities/0',
                                    'derived': f'{prefix}/3/quantities/1',
                                },
                                'quantities': [
                                    {
                                        'name': 'base',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/1',
                                        },
                                        'shape': [],
                                    },
                                    {
                                        'name': 'derived',
                                        'type': {
                                            'type_kind': 'reference',
                                            'type_data': f'{prefix}/2',
                                        },
                                        'shape': [],
                                    },
                                ],
                            },
                        ]
                    }
                },
            },
            id='get-derived-resolved-with-depth',
        ),
    ],
)
def test_definition_reader(query, result):
    with DefinitionReader(query) as reader:
        response = remove_cache(reader.sync_read(m_def))
    assert_dict(response, result)
