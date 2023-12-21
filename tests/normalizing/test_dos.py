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

from nomad.units import ureg
from nomad.datamodel import EntryArchive
from tests.normalizing.conftest import (  # pylint: disable=unused-import
    dos_si_fhiaims,
    dos_si_exciting,
    dos_si_vasp,
    get_template_dos,
    run_normalize,
)

from nomad_dos_fingerprints import DOSFingerprint  # pylint: disable=import-error
from nomad.normalizing.dos_integrator import integrate_dos


def approx(value, abs=0, rel=1e-1):
    return pytest.approx(value, abs=abs, rel=rel)


def test_fingerprint(dos_si_vasp):
    # Check if DOS fingerprint was created
    dos_fingerprint_dict = dos_si_vasp.m_xpath(
        """
        run[*].calculation[*].dos_electronic[*].fingerprint
        """
    )[-1][-1][0]
    dos_fingerprint = DOSFingerprint().from_dict(dos_fingerprint_dict)
    assert dos_fingerprint.get_similarity(dos_fingerprint) == 1
    assert dos_fingerprint.filling_factor != 0
    assert dos_fingerprint.filling_factor != 1


@pytest.mark.parametrize(
    'ranges, highest, lowest, fermi, expected_highest, expected_lowest, n',
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
    ],
)
def test_energy_reference_detection(
    ranges, highest, lowest, fermi, expected_highest, expected_lowest, n
):
    """Checks that the energy reference detection for DOS works in different
    scenarios.
    """
    fermi = fermi[0] if fermi else fermi
    lowest = lowest[0] if lowest else lowest
    highest = highest[0] if highest else highest
    archive = get_template_dos(ranges, fermi, highest, lowest, n)
    assert len(archive.run[0].calculation[0].dos_electronic) == 1
    dos = archive.run[0].calculation[0].dos_electronic[0]
    gap = dos.band_gap[0]
    assert gap.energy_highest_occupied.to(
        ureg.electron_volt
    ).magnitude == pytest.approx(expected_highest[0])
    assert gap.energy_lowest_unoccupied.to(
        ureg.electron_volt
    ).magnitude == pytest.approx(expected_lowest[0])


def test_dos_magnitude(
    dos_si_vasp: EntryArchive,
    dos_si_exciting: EntryArchive,
    dos_si_fhiaims: EntryArchive,
):
    """
    Verify that the raw DOS extracted from similar systems describes the same number of
    electrons. Testing for VASP, exciting and FHI-aims DOS Si2 parsing.
    """
    codes_to_check = list(locals().values())  # TODO add test for normalized DOS
    dos_ints = [
        integrate_dos(dos_si.run[0].calculation[-1].dos_electronic)
        for dos_si in codes_to_check
    ]

    # Compare each DOS with its neighbor
    for index in range(len(dos_ints))[:-1]:
        assert approx(dos_ints[index]) == approx(dos_ints[index + 1])
