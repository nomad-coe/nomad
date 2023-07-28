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

from nomad.datamodel.metainfo.simulation.calculation import (
    Spectra, ElectronicStructureProvenance
)
from nomad.normalizing.normalizer import Normalizer


class SpectraNormalizer(Normalizer):
    """Normalizer for Spectra in run.calculation with the following responsibilities:

      - Sets provenance method section.
      - Normalizes intensities to their maximum value.
    """
    def normalize(self, logger=None) -> None:
        # Setup logger
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # SinglePoint
        if self.entry_archive.m_xpath('run[-1].calculation[-1].spectra'):
            calc_section = self.entry_archive.run[-1].calculation[-1]
            for spectra in calc_section.spectra:
                if self.is_valid_spectra(spectra):
                    # Find photon method provenance
                    if not calc_section.method_ref:
                        return
                    if calc_section.method_ref.photon and not spectra.provenance:
                        provenance = ElectronicStructureProvenance(methodology=calc_section.method_ref, label='photon')
                        spectra.m_add_sub_section(Spectra.provenance, provenance)
                    # Normalizing intensities to their maximum value.
                    # spectra.intensities = spectra.intensities / max(spectra.intensities)
                    # TODO set same excitations_energies reference for OCEAN (BSE), exciting (BSE),
                    # and XSPECTRA (Fermi rule) codes. Generalize for others

    def is_valid_spectra(self, spectra: Spectra) -> bool:
        """Used to check that a spectra has all required information for normalization.
        """
        if spectra.excitation_energies is not None and spectra.intensities is not None:
            n_energies = len(spectra.excitation_energies)
            n_intensities = len(spectra.intensities)
            if n_energies > 0 and n_intensities > 0 and n_energies == n_intensities:
                spectra.n_energies = n_energies
            else:
                self.logger.warning(
                    'Empty arrays or size of arrays do not coincide: could not validate spectra.')
                return False
            if (spectra.intensities < 0.0).any():
                self.logger.warning('Invalid negative intensities found: could not validate spectra.')
                return False
            return True
        else:
            self.logger.warning(
                'Parsed spectra do not contain excitation_energies and intensitites to be normalized.')
            return False
