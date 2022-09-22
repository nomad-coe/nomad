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
    get_template_band_structure,
    band_path_cF,
    band_path_tP,
    band_path_hP,
    band_path_mP_nonstandard,
    band_path_cF_nonstandard,
)

from nomad.units import ureg


@pytest.mark.parametrize('gaps,has_reciprocal_cell', [
    pytest.param([None], True, id="unpolarized, no gap"),
    pytest.param([(1, 'indirect')], True, id="unpolarized, finite gap"),
    pytest.param([None, None], True, id="polarized, no gap"),
    pytest.param([(1, 'indirect'), (0.8, 'indirect')], True, id="polarized, different gaps"),
    pytest.param([(1, 'indirect')], False, id="no reciprocal cell"),
])
def test_band_gaps(gaps, has_reciprocal_cell):
    """Tests that band gaps are correctly identified for different cases.
    """
    bs = get_template_band_structure(gaps, has_reciprocal_cell=has_reciprocal_cell).run[0].calculation[0].band_structure_electronic[0]
    if has_reciprocal_cell:
        assert bs.reciprocal_cell.shape == (3, 3)
    else:
        assert bs.reciprocal_cell is None
    assert len(bs.band_gap) == len(gaps)
    for index, gap in enumerate(gaps):
        channel_info = bs.band_gap[index]
        value = channel_info.value
        if gap is not None and has_reciprocal_cell:
            assert channel_info.type == gap[1]
        else:
            assert channel_info.type is None
        if gap is None:
            assert value == 0
        else:
            gap_ev = value.to(ureg.eV).magnitude
            assert gap_ev == pytest.approx(gap[0], 0.001)
        eho_ev = channel_info.energy_highest_occupied.to(ureg.eV).magnitude
        assert eho_ev == pytest.approx(1 if gap else 0, 0.001)


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
