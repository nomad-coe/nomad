# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the 'License');
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an'AS IS' BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import Any, Dict
import numpy as np
import re
import ase.data

from nomad.normalizing.normalizer import SystemBasedNormalizer
from nomad.metainfo import units
from nomad.metainfo.optimade import OptimadeEntry, Species

species_re = re.compile(r'^([A-Z][a-z]?)(\d*)$')


class OptimadeNormalizer(SystemBasedNormalizer):

    """
    This normalizer performs all produces a section all data necessary for the Optimade API.
    It assumes that the :class:`SystemNormalizer` was run before.
    """
    def __init__(self, backend):
        super().__init__(backend, only_representatives=True)

    def get_optimade_data(self, index) -> OptimadeEntry:
        """
        The 'main' method of this :class:`SystemBasedNormalizer`.
        Normalizes the section with the given `index`.
        Normalizes geometry, classifies, system_type, and runs symmetry analysis.
        """
        optimade = OptimadeEntry()

        def get_value(key: str, default: Any = None, numpy: bool = False) -> Any:
            try:
                value = self._backend.get_value(key, index)
                if type(value) == np.ndarray and not numpy:
                    return value.tolist()
                if isinstance(value, list) and numpy:
                    return np.array(value)

                return value
            except KeyError:
                return default

        from nomad.normalizing.system import normalized_atom_labels

        nomad_species = get_value('atom_labels')

        # elements
        atoms = normalized_atom_labels(nomad_species)
        atom_count = len(atoms)
        atom_counts: Dict[str, int] = {}
        for atom in atoms:
            current = atom_counts.setdefault(atom, 0)
            current += 1
            atom_counts[atom] = current

        optimade.elements = list(set(atoms))
        optimade.elements.sort()
        optimade.nelements = len(optimade.elements)
        optimade.elements_ratios = [
            atom_counts[element] / atom_count
            for element in optimade.elements]

        # formulas
        optimade.chemical_formula_reduced = get_value('chemical_composition_reduced')
        optimade.chemical_formula_hill = get_value('chemical_composition_bulk_reduced')
        optimade.chemical_formula_descriptive = optimade.chemical_formula_hill
        optimade.chemical_formula_anonymous = ''.join([
            '%s' % element + (str(atom_counts[element]) if atom_counts[element] > 1 else '')
            for element in optimade.elements])

        # sites
        optimade.nsites = len(nomad_species)
        optimade.species_at_sites = nomad_species
        optimade.lattice_vectors = (get_value('lattice_vectors', numpy=True) * units.m).to(units.angstrom).magnitude
        optimade.cartesian_site_positions = (get_value('atom_positions', numpy=True) * units.m).to(units.angstrom).magnitude
        optimade.dimension_types = [
            1 if value else 0
            for value in get_value('configuration_periodic_dimensions')]

        # species
        for species_label in set(nomad_species):
            match = re.match(species_re, species_label)

            element_label, index = match.groups(1), match.groups(2)

            species = optimade.m_create(Species)
            species.name = species_label
            if element_label in ase.data.chemical_symbols:
                chemical_symbol = element_label
            else:
                chemical_symbol = 'x'
            species.chemical_symbols = [chemical_symbol]
            species.concentration = [1.0]

        optimade.structure_features = []

        return optimade

    def normalize_system(self, index, is_representative):
        if not is_representative:
            return False

        try:
            optimade = self.get_optimade_data(index)
            self._backend.add_mi2_section(optimade)
        except Exception as e:
            import traceback
            traceback.print_exc()
            self.logger.warn('could not acquire optimade data', exc_info=e)
