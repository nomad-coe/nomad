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

from nomad.datamodel.metainfo.simulation.calculation import Dos
from nomad.utils import get_logger
import numpy as np
from typing import List, Union


# TODO add tester for all functionalities

logger = get_logger(__name__)  # set logger


def get_fermi_energy(dos_object, efermi):
    """
    Handle extraction of the Fermi level
    """
    if efermi is not None:
        e_fermi = efermi
    else:
        try:
            e_fermi = dos_object.energy_fermi
        except AttributeError:
            logger.warning('Missing Fermi energy')
            return

    return e_fermi


def get_energy_index(dos_object, energy_level):
    """
    Obtain the closest index that contains all electron states up to `energy_level`. If the
    closest energy state index cannot be resolved, the function returns:
        - the maximum of energies if energy_level > max(dos_energies)
        - the minimum of energies if energy_level < min(dos_energies)
    """
    if hasattr(energy_level, 'magnitude'):
        dos_energies = dos_object.energies.to(energy_level.units).magnitude  # ensure that the units are correct
        energy_level = energy_level.magnitude  # convert to float if necessary
    else:
        dos_energies = dos_object.energies.magnitude  # now it's up to the user to ensure correct units

    closest_index = np.where(dos_energies >= energy_level)
    if len(closest_index[0]) > 0:
        return closest_index[0][0]
    else:
        logger.warning('Could not find closest_index for the energy_level.')
        if energy_level > np.max(dos_energies):
            return len(dos_energies)
        else:
            return 0


def integrate_dos(
        dos_object: Union[List[Dos], None],
        spin_channels: List[int] = [0, 1],
        limits: List[str] = ['min', 'fermi'],
        efermi: Union[float, None] = None):
    """
    Integrate a NOMAD run `dos_object` over the stated `spin_channels`. In non-normalized cases, this simply yields the number of band electrons.
    - `limits`: 2-object array determining the integration range in energy units. Outside of explicit values, one can also choose `min`, `max` and `fermi`.
    - `efermi`: explicitly passed Fermi level. To be used when the DOS object does not contain any Fermi level itself.
    """
    if len(limits) != 2:
        logger.warning(f'Expected a list of length 2, but got {len(limits)}')
        return

    dos_integrated = 0.
    for spin_channel in spin_channels:
        try:
            dos_spin = dos_object[spin_channel]
        except Exception:
            continue
        # Setting the integral limits
        limit_keywords_map = {
            'min': 0,
            'max': -1,
            'fermi': get_energy_index(dos_spin, get_fermi_energy(dos_spin, efermi))}
        try:
            mapped_limits = [limit_keywords_map[limit] for limit in limits]
        except KeyError:
            mapped_limits = [get_energy_index(dos_spin, limit) for limit in limits]

        # Extract energies and DOS values to perform the integration
        sel_energies = dos_spin.energies.magnitude[mapped_limits[0]:mapped_limits[1]]
        try:
            dos_values = dos_spin.total[0].value.magnitude[mapped_limits[0]:mapped_limits[1]]
            dos_integrated += np.trapz(x=sel_energies, y=dos_values)
        except IndexError:
            continue

    return dos_integrated
