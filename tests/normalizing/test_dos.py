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
# import matplotlib.pyplot as mpl

from nomad.parsing.legacy import Backend
from tests.normalizing.conftest import (  # pylint: disable=unused-import
    dos_si_fhiaims,
    dos_si_exciting,
    dos_si_vasp,
)

from nomad_dos_fingerprints import DOSFingerprint


def test_fingerprint(dos_unpolarized_vasp):
    # Check if DOS fingerprint was created
    backend_dos_fingerprint = dos_unpolarized_vasp.get_value('section_dos_fingerprint', 0)
    dos_fingerprint = DOSFingerprint().from_dict(dict(
        bins=backend_dos_fingerprint.bins,
        indices=backend_dos_fingerprint.indices,
        grid_id=backend_dos_fingerprint.grid_id,
        stepsize=backend_dos_fingerprint.stepsize,
        filling_factor=backend_dos_fingerprint.filling_factor))
    assert dos_fingerprint.get_similarity(dos_fingerprint) == 1
    assert dos_fingerprint.filling_factor != 0
    assert dos_fingerprint.filling_factor != 1


# def test_dos_energies(dos_si_vasp: Backend, dos_si_exciting: Backend, dos_si_fhiaims: Backend):
#  """For debugging.
#  """
#  x_exciting = dos_si_exciting.get_value('dos_energies_normalized', 0)
#  y_exciting = dos_si_exciting.get_value('dos_values_normalized', 0)
#  x_vasp = dos_si_vasp.get_value('dos_energies_normalized', 0)
#  y_vasp = dos_si_vasp.get_value('dos_values_normalized', 0)
#  x_fhiaims = dos_si_fhiaims.get_value('dos_energies_normalized', 0)
#  y_fhiaims = dos_si_fhiaims.get_value('dos_values_normalized', 0)
#  mpl.plot(x_vasp, y_vasp[0], label="VASP")
#  mpl.plot(x_exciting, y_exciting[0], label="exciting")
#  mpl.plot(x_fhiaims, y_fhiaims[0], label="FHI-aims")
#  mpl.legend()
#  mpl.show()


def test_dos_magnitude(dos_si_vasp: Backend, dos_si_exciting: Backend, dos_si_fhiaims: Backend):
    """
    Ensure the DOS normalizer acted on the DOS values. The order of magnitude
    for normalized DOS values in VASP, exciting and FHIAims are currently
    tested.
    """
    dos_vasp = dos_si_vasp.get_value('dos_values_normalized', 0)
    dos_exciting = dos_si_exciting.get_value('dos_values_normalized', 0)
    dos_fhiaims = dos_si_fhiaims.get_value('dos_values_normalized', 0)

    dos_vasp_mean = mean_nonzero(dos_vasp)
    dos_exciting_mean = mean_nonzero(dos_exciting)
    dos_fhiaims_mean = mean_nonzero(dos_fhiaims)

    assert is_same_magnitude(dos_vasp_mean, dos_exciting_mean, dos_fhiaims_mean)


def mean_nonzero(dos: np.array):
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
