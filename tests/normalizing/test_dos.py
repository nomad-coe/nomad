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
# import matplotlib.pyplot as mpl

from nomad.datamodel import EntryArchive
from tests.normalizing.conftest import (  # pylint: disable=unused-import
    dos_si_fhiaims,
    dos_si_exciting,
    dos_si_vasp,
)

from nomad_dos_fingerprints import DOSFingerprint


def test_fingerprint(dos_si_vasp):
    # Check if DOS fingerprint was created
    dos_fingerprint_dict = dos_si_vasp.m_xpath(
        '''
        section_run[*].section_single_configuration_calculation[*].
        section_dos[*].section_dos_fingerprint
        ''')[0][0][0]
    dos_fingerprint = DOSFingerprint().from_dict(dos_fingerprint_dict)
    assert dos_fingerprint.get_similarity(dos_fingerprint) == 1
    assert dos_fingerprint.filling_factor != 0
    assert dos_fingerprint.filling_factor != 1


def test_dos_magnitude(dos_si_vasp: EntryArchive, dos_si_exciting: EntryArchive, dos_si_fhiaims: EntryArchive):
    """
    Ensure the DOS normalizer acted on the DOS values. The order of magnitude
    for normalized DOS values in VASP, exciting and FHIAims are currently
    tested.
    """
    def get_dos_values_normalized(archive):
        return archive.section_run[0].m_xpath(
            'section_single_configuration_calculation[*].section_dos[*].dos_values_normalized')[0][0]

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
    correct_magnitude = 1e18
    tolerance = 10
    values = np.array(args)
    values_normalized = values / correct_magnitude
    return ((values_normalized <= tolerance) & (values_normalized >= 1 / tolerance)).all()
