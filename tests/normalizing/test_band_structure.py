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
        assert bs.reciprocal_cell.shape == (3, 3)

    # Unpolarized, no gaps
    bs = bands_unpolarized_no_gap.run[0].calculation[0].band_structure_electronic[0]
    test_generic(bs)
    assert len(bs.band_gap) == 1
    assert bs.band_gap[0].value == 0

    # Polarized, no gaps
    bs = bands_polarized_no_gap.run[0].calculation[0].band_structure_electronic[0]
    test_generic(bs)
    assert len(bs.band_gap) == 2
    assert bs.band_gap[0].value == 0
    assert bs.band_gap[1].value == 0

    # Unpolarized, finite gap, indirect
    bs = bands_unpolarized_gap_indirect.run[0].calculation[0].band_structure_electronic[0]
    test_generic(bs)
    assert len(bs.band_gap) == 1
    info = bs.band_gap[0]
    gap_joule = info.value
    gap_ev = gap_joule.to(ureg.eV).magnitude
    assert gap_ev == pytest.approx(1, 0.001)
    assert info.type == "indirect"

    # Polarized, finite gap, indirect
    bs = bands_polarized_gap_indirect.run[0].calculation[0].band_structure_electronic[0]
    test_generic(bs)
    assert len(bs.band_gap) == 2
    channel_up = bs.band_gap[0]
    channel_down = bs.band_gap[1]
    gap_up_ev = channel_up.value.to(ureg.eV).magnitude
    gap_down_ev = channel_down.value.to(ureg.eV).magnitude
    assert channel_up.type == "indirect"
    assert channel_down.type == "indirect"
    assert gap_up_ev == pytest.approx(1, 0.01)
    assert gap_down_ev == pytest.approx(0.8, 0.01)


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
    bs = band_path_cF.run[0].calculation[0].band_structure_electronic[0]
    for i, segment in enumerate(bs.segment):
        labels = segment.endpoints_labels
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
    bs = band_path_tP.run[0].calculation[0].band_structure_electronic[0]
    for i, segment in enumerate(bs.segment):
        labels = segment.endpoints_labels
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
    bs = band_path_hP.run[0].calculation[0].band_structure_electronic[0]
    for i, segment in enumerate(bs.segment):
        labels = segment.endpoints_labels
        assert np.array_equal(labels, assumed_labels[i, :])


def test_non_standard(band_path_mP_nonstandard, band_path_cF_nonstandard):
    """Tests for lattice that do not follow the Setyawan/Curtarolo standard.
    """
    # The ordering of the lattice does not follow the standard: a, b <= c. No labels defined.
    bs = band_path_mP_nonstandard.run[0].calculation[0].band_structure_electronic[0]
    for segment in bs.segment:
        labels = segment.endpoints_labels
        assert labels is None

    # Existing labels should not be overridden
    assumed_labels = np.array([
        ["W", "L"],
        ["L", "Γ"],
        ["Γ", "X"],
        ["X", "W"],
        ["W", "K"],
    ])
    bs = band_path_cF_nonstandard.run[0].calculation[0].band_structure_electronic[0]
    for i, segment in enumerate(bs.segment):
        labels = segment.endpoints_labels
        assert np.array_equal(labels, assumed_labels[i, :])


def test_phonon_band(phonon):
    """Ensures that band gaps are not added to phonon bands.
    """
    bs = phonon.run[0].calculation[0].band_structure_phonon[0]
    assert bs.path_standard is None
    assert len(bs.band_gap) == 0
