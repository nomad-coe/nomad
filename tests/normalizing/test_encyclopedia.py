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

from nomad.metainfo.encyclopedia import Encyclopedia
from tests.normalizing.conftest import geometry_optimization, molecular_dynamics   # pylint: disable=unused-import


def test_geometry_optimization(geometry_optimization: Encyclopedia):
    """Tests that geometry optimizations are correctly processed."
    """
    run_type = geometry_optimization.calculation.run_type
    assert run_type == "geometry optimization"


# def test_molecular_dynamics(molecular_dynamics: Encyclopedia):
    # """Tests that geometry optimizations are correctly processed."
    # """
    # run_type = molecular_dynamics.calculation.run_type
    # assert run_type == "molecular dynamics"


def test_system_type(geometry_optimization: Encyclopedia):
    """Tests that geometry optimizations are correctly processed."
    """
    system_type = geometry_optimization.material.system_type
    assert system_type == "bulk"
