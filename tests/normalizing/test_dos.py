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

from nptyping import NDArray
import numpy as np
import pytest

from nomad.units import ureg
from nomad.datamodel import EntryArchive
from tests.normalizing.conftest import (  # pylint: disable=unused-import
    dos_si_fhiaims,
    dos_si_exciting,
    dos_si_vasp,
    get_template_dos,
    run_normalize,
)

from nomad_dos_fingerprints import DOSFingerprint


def test_fingerprint(dos_si_vasp):
    # Check if DOS fingerprint was created
    dos_fingerprint_dict = dos_si_vasp.m_xpath(
        '''
        run[*].calculation[*].dos_electronic[*].fingerprint
        ''')[-1][-1][0]
    dos_fingerprint = DOSFingerprint().from_dict(dos_fingerprint_dict)
    assert dos_fingerprint.get_similarity(dos_fingerprint) == 1
    assert dos_fingerprint.filling_factor != 0
    assert dos_fingerprint.filling_factor != 1


@pytest.mark.parametrize(
    "ranges, highest, lowest, fermi, expected_highest, expected_lowest, n",
    [
        # Explicit highest/lowest occupied given by parser: The current
        # behaviour is to override these values based on the data that is
        # actually stored in the DOS if there is a minor difference.
        ([[[0, 1], [2, 3]]], [1.04], [1.94], None, [1], [1.9], 101),
        # Fermi energy in the middle of a gap, inaccuracy due to discretization.
        ([[[0, 1], [2, 3]]], None, None, [1.5], [1.0], [1.9], 101),
        # Fermi energy near the top of the valence band. Since Fermi energy
        # is close enough to the zero value, gap is detected.
        ([[[0, 1], [2, 3]]], None, None, [0.99], [1.0], [1.9], 101),
        # Fermi energy near the top of the valence band, but still too far away
        # for a gap to be detected.
        ([[[0, 1], [2, 3]]], None, None, [0.89], [0.9], [0.9], 101),
        # Fermi energy near the bottom of the conduction band. Since Fermi energy
        # is close enough to the zero value, gap is detected.
        ([[[0, 1], [2, 3]]], None, None, [1.91], [1.0], [1.9], 101),
        # Fermi energy near the bottom of the conduction band, but still too
        # far away for a gap to be detected.
        ([[[0, 1], [2, 3]]], None, None, [2.01], [2.0], [2.0], 101),
        # Fermi energy at the center of a tiny peak.
        ([[[1, 1.1]]], None, None, [1], [1], [1], 101),
    ]
)
def test_energy_reference_detection(ranges, highest, lowest, fermi, expected_highest, expected_lowest, n):
    """Checks that the energy reference detection for DOS works in different
    scenarios.
    """
    fermi = fermi[0] if fermi else fermi
    lowest = lowest[0] if lowest else lowest
    highest = highest[0] if highest else highest
    archive = get_template_dos(ranges, fermi, highest, lowest, n)
    dos = archive.run[0].calculation[0].dos_electronic[0]
    n_channels = len(dos.total)
    for i_channel in range(n_channels):
        gap = dos.band_gap[i_channel]
        assert gap.energy_highest_occupied.to(ureg.electron_volt).magnitude == pytest.approx(
            expected_highest[i_channel]
        )
        assert gap.energy_lowest_unoccupied.to(ureg.electron_volt).magnitude == pytest.approx(
            expected_lowest[i_channel]
        )


@pytest.mark.skip('Temporarly disabled, but needs to be fixed.')
def test_dos_magnitude(dos_si_vasp: EntryArchive, dos_si_exciting: EntryArchive, dos_si_fhiaims: EntryArchive):
    """
    Ensure the DOS normalizer acted on the DOS values. The order of magnitude
    for normalized DOS values in VASP, exciting and FHIAims are currently
    tested.
    """
    def get_dos_values_normalized(archive):
        total_dos = archive.run[0].calculation[-1].dos_electronic[-1].total
        return np.array([d.value.to('1/eV') * d.normalization_factor for d in total_dos])

    dos_vasp = get_dos_values_normalized(dos_si_vasp)
    dos_exciting = get_dos_values_normalized(dos_si_exciting)
    dos_fhiaims = get_dos_values_normalized(dos_si_fhiaims)

    dos_vasp_mean = mean_nonzero(dos_vasp)
    dos_exciting_mean = mean_nonzero(dos_exciting)
    dos_fhiaims_mean = mean_nonzero(dos_fhiaims)

    assert is_same_magnitude(dos_vasp_mean, dos_exciting_mean, dos_fhiaims_mean)


def mean_nonzero(dos: NDArray):
    """Returns the mean value of all nonzero elements in the given array.
    """
    return dos[np.nonzero(dos)].mean()


def is_same_magnitude(*args):
    """Used to test that all given floating point numbers are of the expected
    order of magnitude.
    """
    correct_magnitude = 1e-2
    tolerance = 10
    values = np.array(args)
    values_normalized = values / correct_magnitude
    return ((values_normalized <= tolerance) & (values_normalized >= 1 / tolerance)).all()
