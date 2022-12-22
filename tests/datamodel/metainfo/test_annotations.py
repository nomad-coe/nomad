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

from nomad.datamodel.metainfo.annotations import PlotAnnotation


@pytest.mark.parametrize('annotation, result, error', [
    pytest.param({'x': 'x', 'y': 'y'}, None, False, id='plain'),
    pytest.param({
        'x': ['1', '2'],
        'y': ['1', '2'],
        'label': 'test',
        'lines': [{}, {}],
        'layout': {},
        'config': {}
    }, None, False, id='full'),
    pytest.param({
        'x': ['1', '2'],
        'y': ['1'],
        'label': 'test',
        'lines': [{}, {}],
        'layout': {},
        'config': {}
    }, None, True, id='x-y-mismatch'),
    pytest.param({'y': ['1']}, None, True, id='x-and-y-required-y'),
    pytest.param({'x': ['1']}, None, True, id='x-and-y-required-x'),
    pytest.param({}, None, True, id='x-and-y-required'),
    pytest.param({'x_axis': '1', 'y_axis': '1'}, {'x': '1', 'y': '1'}, False, id='deprecated-aliases'),
    pytest.param({'x': './a/0/b', 'y': 'a/b'}, None, False, id='x-y-pattern-ok'),
    pytest.param({'x': './a', 'y': 'a/../b'}, None, True, id='x-y-pattern-fail'),
])
def test_plot_validation(annotation, result, error):
    if error:
        with pytest.raises(ValidationError):
            PlotAnnotation(**annotation)
    else:
        assert PlotAnnotation(**annotation).dict(exclude_none=True) == (result or annotation)
