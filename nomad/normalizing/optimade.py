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

from nomad.datamodel import EntryArchive
from nomad.atomutils import Formula
from nomad.normalizing.normalizer import SystemBasedNormalizer
from nomad.units import ureg
from nomad.datamodel import OptimadeEntry, Species, EntryMetadata


species_re = re.compile(r'^([A-Z][a-z]?)(\d*)$')
atom_label_re = re.compile(
    '|'.join(sorted(ase.data.chemical_symbols, key=lambda x: len(x), reverse=True))
)


def normalized_atom_labels(atom_labels):
    """
    Normalizes the given atom labels: they either are labels right away, or contain
    additional numbers (to distinguish same species but different labels, see meta-info),
    or we replace them with ase placeholder atom for unknown elements 'X'.
    """
    return [
        ase.data.chemical_symbols[0] if match is None else match.group(0)
        for match in [
            re.search(atom_label_re, atom_label) for atom_label in atom_labels
        ]
    ]


# TODO this should be the default and not necessary
def transform_to_v1(entry: EntryMetadata) -> EntryMetadata:
    """
    Transformation function to use during re-indexing of entries with outdated optimade
    format. Fixes formulas and periodic dimensions, removed entries with X in formula.
    """
    optimade = entry.optimade
    if optimade is None:
        return entry

    if (
        optimade.chemical_formula_reduced is None
        or 'X' in optimade.chemical_formula_reduced
    ):
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
    """
    This normalizer performs all produces a section all data necessary for the Optimade API.
    It assumes that the :class:`SystemNormalizer` was run before.
    """

    normalizer_level = 1

    def __init__(self):
        super().__init__(only_representatives=True)

    def add_optimade_data(self, archive: EntryArchive) -> OptimadeEntry:
        """
        The 'main' method of this :class:`SystemBasedNormalizer`.
        Normalizes the section with the given `index`.
        Normalizes geometry, classifies, system_type, and runs symmetry analysis.
        """
        if archive.metadata is None:
            archive.m_create(EntryMetadata)
        optimade = archive.metadata.m_create(OptimadeEntry)

        def get_value(
            quantity_def,
            default: Any = None,
            numpy: bool = False,
            unit=None,
            source: Any = None,
        ) -> Any:
            try:
                source = source if source is not None else archive.run[0].system[-1]
                value = source.m_get(quantity_def)
                if value is None:
                    return
                if type(value) is np.ndarray and not numpy:
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

        system = archive.run[0].system[-1] if archive.run[0].system else None
        if system is None:
            return optimade

        atoms_cls = system.m_def.all_sub_sections['atoms'].sub_section.section_cls
        nomad_species = get_value(atoms_cls.labels, source=system.atoms)
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
            atom_counts[element] / atom_count for element in optimade.elements
        ]

        # formulas
        system_cls = (
            archive.run[0].m_def.all_sub_sections['system'].sub_section.section_cls
        )
        original_formula = get_value(
            system_cls.chemical_composition_hill, source=system
        )
        if original_formula is not None:
            formula = Formula(original_formula)
            optimade.chemical_formula_reduced = formula.format('reduced')
            optimade.chemical_formula_hill = formula.format('hill')
            optimade.chemical_formula_anonymous = formula.format('anonymous')
            optimade.chemical_formula_descriptive = optimade.chemical_formula_hill

        # sites
        optimade.nsites = len(nomad_species)
        optimade.species_at_sites = nomad_species
        optimade.lattice_vectors = get_value(
            atoms_cls.lattice_vectors, numpy=True, unit=ureg.m, source=system.atoms
        )
        optimade.cartesian_site_positions = get_value(
            atoms_cls.positions, numpy=True, unit=ureg.m, source=system.atoms
        )
        pbc = get_value(atoms_cls.periodic, source=system.atoms)
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

    def normalize_system(self, archive: EntryArchive, system, is_representative):
        if not is_representative:
            return False

        try:
            self.add_optimade_data(archive)
            return True

        except Exception as e:
            self.logger.warn('could not acquire optimade data', exc_info=e)
