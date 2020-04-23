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

import numpy as np
import pytest

from tests.normalizing.conftest import (  # pylint: disable=unused-import
    phonon,
    bands_unpolarized_no_gap,
    bands_polarized_no_gap,
    bands_unpolarized_gap_indirect,
    bands_polarized_gap_indirect,
    band_path_fcc,
)
from pint import UnitRegistry
ureg = UnitRegistry()


def test_band_gaps(bands_unpolarized_no_gap, bands_polarized_no_gap, bands_unpolarized_gap_indirect, bands_polarized_gap_indirect):
    """Tests that band gaps are correctly identified for different cases.
    """
    def test_generic(bs, n_channels):
        """Generic tests for band structure data."""
        assert bs.brillouin_zone is not None
        assert bs.reciprocal_cell.shape == (3, 3)

    # Unpolarized, no gaps
    bs = bands_unpolarized_no_gap.entry_archive.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    test_generic(bs, n_channels=1)
    assert bs.section_band_gap is None
    assert bs.section_band_gap_spin_up is None
    assert bs.section_band_gap_spin_down is None

    # Polarized, no gaps
    bs = bands_polarized_no_gap.entry_archive.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    test_generic(bs, n_channels=2)
    assert bs.section_band_gap is None
    assert bs.section_band_gap_spin_up is None
    assert bs.section_band_gap_spin_down is None

    # Unpolarized, finite gap, indirect
    bs = bands_unpolarized_gap_indirect.entry_archive.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    test_generic(bs, n_channels=1)
    gap_ev = (bs.section_band_gap.value * ureg.J).to(ureg.eV).magnitude
    assert gap_ev == pytest.approx(0.62, 0.01)
    assert bs.section_band_gap.type == "indirect"
    assert bs.section_band_gap_spin_up is None
    assert bs.section_band_gap_spin_down is None

    # Polarized, finite gap, indirect
    bs = bands_polarized_gap_indirect.entry_archive.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    test_generic(bs, n_channels=2)
    gap = bs.section_band_gap
    gap_up = bs.section_band_gap_spin_up
    gap_down = bs.section_band_gap_spin_down
    gap_ev = (gap.value * ureg.J).to(ureg.eV).magnitude
    gap_up_ev = (gap_up.value * ureg.J).to(ureg.eV).magnitude
    gap_down_ev = (gap_down.value * ureg.J).to(ureg.eV).magnitude
    assert gap_up.type == "indirect"
    assert gap_down.type == "indirect"
    assert gap_up_ev != gap_down_ev
    assert gap_up_ev == gap_ev
    assert gap_up_ev == pytest.approx(0.956, 0.01)
    assert gap_down_ev == pytest.approx(1.230, 0.01)


def test_band_paths(band_path_fcc):
    """Tests that the path labeling workds for different lattices.
    """
    # FCC
    scc = band_path_fcc.entry_archive.section_run[0].section_single_configuration_calculation[0]
    system = scc.single_configuration_calculation_to_system_ref
    space_group_number = system.section_symmetry[0].space_group_number
    assert space_group_number == 227
    assumed_path = np.array([
        ["G", "X"],
        ["X", "W"],
        ["W", "K"],
        ["K", "G"],
        ["G", "L"],
        ["L", "U"],
        ["U", "W"],
        ["W", "L"],
        ["L", "K"],
        ["U", "X"],
    ])
    bs = band_path_fcc.entry_archive.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    for i, segment in enumerate(bs.section_k_band_segment):
        labels = segment.band_path_labels
        assert np.array_equal(labels, assumed_path[i, :])


def test_phonon_band(phonon):
    """Ensures that phonon bands are not touched.
    """
    bs = phonon.entry_archive.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    assert bs.is_standard_path is None
    assert bs.section_band_gap is None
    assert bs.section_band_gap_spin_up is None
    assert bs.section_band_gap_spin_down is None
