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
from pydantic import ValidationError

from nomad.metainfo import Quantity
from nomad.datamodel.metainfo.annotations import PlotAnnotation, ELNAnnotation


@pytest.mark.parametrize('quantity, annotation, result, error', [
    pytest.param(Quantity(type=str), {}, {}, False, id='empty'),
    pytest.param(Quantity(type=str), {'component': 'StringEditQuantity'}, None, False, id='plain'),
    pytest.param(Quantity(type=str), {'component': 'NumberEditQuantity'}, None, True, id='wrong-type'),
    pytest.param(Quantity(type=str, shape=[1, 2]), {'component': 'StringEditQuantity'}, None, True, id='bad-shape')
])
def test_eln_validation(quantity, annotation, result, error):
    if error:
        with pytest.raises(ValidationError):
            annotation_model = ELNAnnotation(**annotation)
            annotation_model.m_definition = quantity
    else:
        annotation_model = ELNAnnotation(**annotation)
        annotation_model.m_definition = quantity
        assert annotation_model.dict(exclude_none=True) == (result or annotation)


@pytest.mark.parametrize('annotation, result, error', [
    pytest.param({'data': {'x': '#x', 'y': '#y'}}, None, False, id='plain'),
    pytest.param({
        'data': {'x': [1, 2], 'y': [1, 2]},
        'layout': {},
        'config': {}
    }, None, False, id='raw-data'),
    pytest.param({
        'data': {'x': [1, 2], 'y': [1, 2], 'z': '#a/b/c'},
        'layout': {},
        'config': {}
    }, None, False, id='mixture'),
    pytest.param({
        'data': {'x': ['#1', '#2'], 'y': ['#1', '#2']},
        'layout': {},
        'config': {}
    }, None, True, id='error-list-of-references'),
    pytest.param({
        'data': {'x': [1, 2], 'y': [1, '2']},
        'layout': {},
        'config': {}
    }, None, True, id='invalid-numbers'),
    pytest.param({
        'data': {'x': ['#1', '#2'], 'y': ['#1']},
        'label': 'test',
        'lines': [{}, {}],
        'layout': {},
        'config': {}
    }, None, True, id='x-y-mismatch'),
    pytest.param({'data': {}}, None, True, id='data-required'),
    pytest.param({'data': {'x': '#./a/0/b', 'y': '#a/b'}}, None, False, id='x-y-pattern-ok'),
    pytest.param({'data': {'x': '#./a', 'y': '#a/../b'}}, None, True, id='x-y-pattern-fail'),
    pytest.param({'data': {'x': '#a/:/b', 'y': '#a/:/b'}}, None, False, id='slice-all'),
    pytest.param({'data': {'x': '#a/1:2/b', 'y': '#a/1:2/b'}}, None, False, id='slice-1:2'),
    pytest.param({'data': {'x': '#a/:/a/b', 'y': '#a/:/a/b'}}, None, False, id='slice-multiple-after'),
    pytest.param({'data': {'x': '#a/-2:-1/b', 'y': '#a/-2:-1/b'}}, None, False, id='slice-negative'),
    pytest.param({'data': {'x': '#a/:/a/:/b', 'y': '#a/:/a/:/b'}}, None, False, id='slice-multiple'),
    pytest.param({'data': {'x': '#:/::/b', 'y': '#:/::/b'}}, None, True, id='slice-too-many'),
    pytest.param({'data': {'x': '#a/:a/b', 'y': '#a/:a/b'}}, None, True, id='slice-invalid-surrounding'),
    pytest.param({'data': {'x': '#a/:/:/b', 'y': '#a/:/:/b'}}, None, True, id='slice-two-in-row'),
    pytest.param({'data': {'x': '#a/:', 'y': '#a/:'}}, None, True, id='slice-last'),
    pytest.param({'data': {'x': '#a/1:2.1/b', 'y': '#a/1:2.1/b'}}, None, True, id='slice-non-integer'),
    pytest.param({'data': {'x': '#a/1:2:5/b', 'y': '#a/1:2:5/b'}}, None, True, id='slice-step-unsupported'),
    pytest.param({'data': [{'x': '#x', 'y': '#y'}, {'x': '#x', 'y': '#y'}]}, None, False, id='list'),
])
def test_plot_validation(annotation, result, error):
    if error:
        with pytest.raises(ValidationError):
            PlotAnnotation(**annotation)
    else:
        assert PlotAnnotation(**annotation).dict(exclude_none=True) == (result or annotation)
