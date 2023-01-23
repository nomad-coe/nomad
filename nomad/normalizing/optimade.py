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

from typing import Any, Dict
import numpy as np
import re
import ase.data
import ase.formula
import pint.quantity

from nomad.atomutils import Formula
from nomad.normalizing.normalizer import SystemBasedNormalizer
from nomad.units import ureg
from nomad.datamodel import OptimadeEntry, Species, EntryMetadata
from nomad.datamodel.metainfo.simulation.system import Atoms, System


species_re = re.compile(r'^([A-Z][a-z]?)(\d*)$')


# TODO this should be the default and not necessary
def transform_to_v1(entry: EntryMetadata) -> EntryMetadata:
    '''
    Transformation function to use during re-indexing of entries with outdated optimade
    format. Fixes formulas and periodic dimensions, removed entries with X in formula.
    '''
    optimade = entry.optimade
    if optimade is None:
        return entry

    if optimade.chemical_formula_reduced is None or 'X' in optimade.chemical_formula_reduced:
        entry.m_remove_sub_section(EntryMetadata.optimade, -1)
        return entry

    if optimade.chemical_formula_hill is not None:
        formula = Formula(optimade.chemical_formula_hill)
        optimade.chemical_formula_reduced = formula.format('reduced')
        optimade.chemical_formula_hill = formula.format('hill')
        optimade.chemical_formula_anonymous = formula.format('anonymous')
        optimade.chemical_formula_descriptive = optimade.chemical_formula_hill

    dimension_types = optimade.dimension_types
    if dimension_types is None:
        optimade.dimension_types = [0, 0, 0]
    elif isinstance(dimension_types, int):
        optimade.dimension_types = [1] * dimension_types + [0] * (3 - dimension_types)

    return entry


class OptimadeNormalizer(SystemBasedNormalizer):

    '''
    This normalizer performs all produces a section all data necessary for the Optimade API.
    It assumes that the :class:`SystemNormalizer` was run before.
    '''
    def __init__(self, archive):
        super().__init__(archive, only_representatives=True)

    def add_optimade_data(self, index) -> OptimadeEntry:
        '''
        The 'main' method of this :class:`SystemBasedNormalizer`.
        Normalizes the section with the given `index`.
        Normalizes geometry, classifies, system_type, and runs symmetry analysis.
        '''
        if self.entry_archive.metadata is None:
            self.entry_archive.m_create(EntryMetadata)
        optimade = self.entry_archive.metadata.m_create(OptimadeEntry)

        def get_value(quantity_def, default: Any = None, numpy: bool = False, unit=None, source: Any = None) -> Any:
            try:
                source = source if source is not None else self.section_run.system[-1]
                value = source.m_get(quantity_def)
                if value is None:
                    return
                if type(value) == np.ndarray and not numpy:
                    return value.tolist()
                if isinstance(value, list) and numpy:
                    return np.array(value)

                if numpy and unit is not None:
                    if isinstance(value, pint.quantity._Quantity):
                        value = value.to(unit)
                    elif value is not None:
                        value = value * unit

                return value
            except KeyError:
                return default

        from nomad.normalizing.system import normalized_atom_labels

        system = self.section_run.system[-1] if self.section_run.system else None
        if system is None:
            return optimade

        nomad_species = get_value(Atoms.labels, source=system.atoms)
        nomad_species = [] if nomad_species is None else nomad_species

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
        original_formula = get_value(System.chemical_composition_hill, source=system)
        if original_formula is not None:
            formula = Formula(original_formula)
            optimade.chemical_formula_reduced = formula.format('reduced')
            optimade.chemical_formula_hill = formula.format('hill')
            optimade.chemical_formula_anonymous = formula.format('anonymous')
            optimade.chemical_formula_descriptive = optimade.chemical_formula_hill

        # sites
        optimade.nsites = len(nomad_species)
        optimade.species_at_sites = nomad_species
        optimade.lattice_vectors = get_value(Atoms.lattice_vectors, numpy=True, unit=ureg.m, source=system.atoms)
        optimade.cartesian_site_positions = get_value(Atoms.positions, numpy=True, unit=ureg.m, source=system.atoms)
        pbc = get_value(Atoms.periodic, source=system.atoms)
        if pbc is not None:
            optimade.dimension_types = [1 if value else 0 for value in pbc]

        # species
        for species_label in set(nomad_species):
            match = re.match(species_re, species_label)

            element_label = match.group(1) if match else species_label

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
