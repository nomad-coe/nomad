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

from tests.normalizing.conftest import (  # pylint: disable=unused-import
    phonon,
    bands_unpolarized_no_gap,
    bands_polarized_no_gap,
    bands_unpolarized_gap_indirect,
    bands_polarized_gap_indirect,
    band_path_cF,
    band_path_tP,
    band_path_hP,
    band_path_mP_nonstandard,
    band_path_cF_nonstandard,
)

from nomad.units import ureg


def test_band_gaps(bands_unpolarized_no_gap, bands_polarized_no_gap, bands_unpolarized_gap_indirect, bands_polarized_gap_indirect):
    """Tests that band gaps are correctly identified for different cases.
    """
    def test_generic(bs):
        """Generic tests for band structure data."""
        assert bs.brillouin_zone is not None
        assert bs.reciprocal_cell.shape == (3, 3)

    # Unpolarized, no gaps
    bs = bands_unpolarized_no_gap.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    test_generic(bs)
    assert len(bs.section_band_gap) == 1
    assert bs.section_band_gap[0].value == 0

    # Polarized, no gaps
    bs = bands_polarized_no_gap.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    test_generic(bs)
    assert len(bs.section_band_gap) == 2
    assert bs.section_band_gap[0].value == 0
    assert bs.section_band_gap[1].value == 0

    # Unpolarized, finite gap, indirect
    bs = bands_unpolarized_gap_indirect.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    test_generic(bs)
    assert len(bs.section_band_gap) == 1
    gap = bs.section_band_gap[0]
    gap_ev = (gap.value * ureg.J).to(ureg.eV).magnitude
    assert gap_ev == pytest.approx(0.62, 0.01)
    assert gap.type == "indirect"

    # TODO: AL I cannot find a polarized example with band gap! Previous parser got the band gap wrong.
    # Polarized, finite gap, indirect
    # bs = bands_polarized_gap_indirect.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    # test_generic(bs)
    # assert len(bs.section_band_gap) == 2
    # gap_up = bs.section_band_gap[0]
    # gap_down = bs.section_band_gap[1]
    # gap_up_ev = (gap_up.value * ureg.J).to(ureg.eV).magnitude
    # gap_down_ev = (gap_down.value * ureg.J).to(ureg.eV).magnitude
    # assert gap_up.type == "indirect"
    # assert gap_down.type == "indirect"
    # assert gap_up_ev == pytest.approx(0.956, 0.01)
    # assert gap_down_ev == pytest.approx(1.230, 0.01)


def test_paths(band_path_cF, band_path_tP, band_path_hP):
    """Tests that the paths are labeled correctly.
    """
    # Face-centered cubic (FCC, cF)
    assumed_labels = np.array([
        ["Γ", "X"],
        ["X", "W"],
        ["W", "K"],
        ["K", "Γ"],
        ["Γ", "L"],
        ["L", "U"],
        ["U", "W"],
        ["W", "L"],
        ["L", "K"],
        ["U", "X"],
    ])
    bs = band_path_cF.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    for i, segment in enumerate(bs.section_k_band_segment):
        labels = segment.band_segm_labels
        assert np.array_equal(labels, assumed_labels[i, :])

    # Tetragonal (TET, tP)
    assumed_labels = np.array([
        ["Γ", "X"],
        ["X", "M"],
        ["M", "Γ"],
        ["Γ", "Z"],
        ["Z", "R"],
        ["R", "A"],
        ["A", "Z"],
        ["X", "R"],
        ["M", "A"],
    ])
    bs = band_path_tP.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    for i, segment in enumerate(bs.section_k_band_segment):
        labels = segment.band_segm_labels
        assert np.array_equal(labels, assumed_labels[i, :])

    # Hexagonal (HEX, hP)
    assumed_labels = np.array([
        ["Γ", "M"],
        ["M", "K"],
        ["K", "Γ"],
        ["Γ", "A"],
        ["A", "L"],
        ["L", "H"],
        ["H", "A"],
        ["L", "M"],
        ["K", "H"],
    ])
    bs = band_path_hP.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    for i, segment in enumerate(bs.section_k_band_segment):
        labels = segment.band_segm_labels
        assert np.array_equal(labels, assumed_labels[i, :])


def test_non_standard(band_path_mP_nonstandard, band_path_cF_nonstandard):
    """Tests for lattice that do not follow the Setyawan/Curtarolo standard.
    """
    # The ordering of the lattice does not follow the standard: a, b <= c. Not labels defined.
    bs = band_path_mP_nonstandard.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    for segment in bs.section_k_band_segment:
        labels = segment.band_segm_labels
        assert labels is None

    # Existing labels should not be overridden
    assumed_labels = np.array([
        ["W", "L"],
        ["L", "Γ"],
        ["Γ", "X"],
        ["X", "W"],
        ["W", "K"],
    ])
    bs = band_path_cF_nonstandard.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    for i, segment in enumerate(bs.section_k_band_segment):
        labels = segment.band_segm_labels
        assert np.array_equal(labels, assumed_labels[i, :])


def test_phonon_band(phonon):
    """Ensures that band gaps are not added to phonon bands.
    """
    bs = phonon.section_run[0].section_single_configuration_calculation[0].section_k_band[0]
    assert bs.is_standard_path is None
    assert len(bs.section_band_gap) == 0
