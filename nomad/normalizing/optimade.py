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
from string import ascii_uppercase
import pint.quantity

from nomad.normalizing.normalizer import SystemBasedNormalizer
from nomad.units import ureg
from nomad.datamodel import OptimadeEntry, Species, DFTMetadata, EntryMetadata
from nomad.datamodel.metainfo.public import section_system

species_re = re.compile(r'^([A-Z][a-z]?)(\d*)$')


class OptimadeNormalizer(SystemBasedNormalizer):

    '''
    This normalizer performs all produces a section all data necessary for the Optimade API.
    It assumes that the :class:`SystemNormalizer` was run before.
    '''
    def __init__(self, backend):
        super().__init__(backend, only_representatives=True)

    def add_optimade_data(self, index) -> OptimadeEntry:
        '''
        The 'main' method of this :class:`SystemBasedNormalizer`.
        Normalizes the section with the given `index`.
        Normalizes geometry, classifies, system_type, and runs symmetry analysis.
        '''

        if self.entry_archive.section_metadata is None:
            self.entry_archive.m_create(EntryMetadata)
        if self.entry_archive.section_metadata.dft is None:
            self.entry_archive.section_metadata.m_create(DFTMetadata)
        optimade = self.entry_archive.section_metadata.dft.m_create(OptimadeEntry)

        def get_value(quantity_def, default: Any = None, numpy: bool = False, unit=None) -> Any:
            try:
                value = self.section_run.section_system[-1].m_get(quantity_def)
                if type(value) == np.ndarray and not numpy:
                    return value.tolist()
                if isinstance(value, list) and numpy:
                    return np.array(value)

                if numpy and unit is not None:
                    if isinstance(value, pint.quantity._Quantity):
                        value = value.to(unit)
                    else:
                        value = value * unit

                return value
            except KeyError:
                return default

        from nomad.normalizing.system import normalized_atom_labels

        nomad_species = get_value(section_system.atom_labels)

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
        optimade.chemical_formula_reduced = get_value(section_system.chemical_composition_reduced)
        optimade.chemical_formula_hill = get_value(section_system.chemical_composition_bulk_reduced)
        optimade.chemical_formula_descriptive = optimade.chemical_formula_hill
        optimade.chemical_formula_anonymous = ''
        for i in range(len(optimade.elements)):
            part = '%s' % ascii_uppercase[i % len(ascii_uppercase)]
            if atom_counts[optimade.elements[i]] > 1:
                part += str(atom_counts[optimade.elements[i]])
            optimade.chemical_formula_anonymous += part

        # sites
        optimade.nsites = len(nomad_species)
        optimade.species_at_sites = nomad_species
        optimade.lattice_vectors = get_value(section_system.lattice_vectors, numpy=True, unit=ureg.m)
        optimade.cartesian_site_positions = get_value(section_system.atom_positions, numpy=True, unit=ureg.m)
        optimade.dimension_types = [
            1 if value else 0
            for value in get_value(section_system.configuration_periodic_dimensions)]

        # species
        for species_label in set(nomad_species):
            match = re.match(species_re, species_label)

            element_label = match.group(1)

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

    def normalize_system(self, system, is_representative):
        if not is_representative:
            return False

        try:
            self.add_optimade_data(system.m_parent_index)
            return True

        except Exception as e:
            self.logger.warn('could not acquire optimade data', exc_info=e)
