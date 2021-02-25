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

import json
import numpy as np
from typing import Dict, List

from nomad import atomutils
from nomad.normalizing.normalizer import Normalizer
from nomad.datamodel.encyclopedia import EncyclopediaMetadata
from nomad.datamodel.optimade import OptimadeEntry, Species
from nomad.datamodel.metainfo.common_dft import section_symmetry, section_system
from nomad.datamodel.results import (
    Results,
    Material,
    Symmetry,
    Structure,
    LatticeParameters,
    WyckoffSet,
)


class ResultsNormalizer(Normalizer):
    """Populates the results section in the metainfo.
    """
    def normalize(self, logger=None) -> None:
        # Setup logger
        if logger is not None:
            self.logger = logger.bind(normalizer=self.__class__.__name__)

        # Do nothing if section_run is not present
        if self.section_run is None:
            return

        # Fetch different information resources from which data is gathered
        repr_sys = None
        for section in self.section_run.section_system:
            if section.is_representative:
                repr_sys = section
                break
        try:
            encyclopedia = self.entry_archive.section_metadata.encyclopedia
        except Exception:
            encyclopedia = None
        try:
            optimade = self.entry_archive.section_metadata.dft.optimade
        except Exception:
            optimade = None

        # Create the section and populate the subsections
        results = Results()
        results.material = self.material(repr_sys, encyclopedia, optimade)
        self.entry_archive.results = results

    def material(
            self,
            repr_sys: section_system,
            encyclopedia: EncyclopediaMetadata,
            optimade: OptimadeEntry) -> Material:
        """Returns a populated Material subsection."""
        material = Material()
        symm = repr_sys.section_symmetry if repr_sys else None
        symm = symm[0] if symm else None

        if repr_sys:
            material.type_structural = repr_sys.system_type
            names, counts = atomutils.get_hill_decomposition(repr_sys.atom_labels, reduced=True)
            material.chemical_formula_reduced_fragments = [
                "{}{}".format(n, int(c) if c != 1 else "") for n, c in zip(names, counts)
            ]
            struct_orig = self.structure_original(repr_sys)
            if struct_orig:
                material.structure_original = struct_orig
        if encyclopedia:
            material.material_id = encyclopedia.material.material_id
            classes = encyclopedia.material.material_classification
            if classes:
                classifications = json.loads(classes)
                material.type_functional = classifications.get("material_class_springer")
                material.type_compound = classifications.get("compound_class_springer")
            material.material_name = encyclopedia.material.material_name
        if optimade:
            material.elements = optimade.elements
            material.nelements = optimade.nelements
            material.chemical_formula_descriptive = optimade.chemical_formula_descriptive
            material.chemical_formula_reduced = optimade.chemical_formula_reduced
            material.chemical_formula_hill = optimade.chemical_formula_hill
            material.chemical_formula_anonymous = optimade.chemical_formula_anonymous

        if symm:
            symmetry = self.symmetry(repr_sys, symm, encyclopedia)
            if symmetry:
                material.symmetry = symmetry
            struct_prim = self.structure_primitive(symm)
            if struct_prim:
                material.structure_primitive = struct_prim
            struct_conv = self.structure_conventional(symm)
            if struct_conv:
                material.structure_conventional = struct_conv

        return material

    def symmetry(
            self,
            repr_sys: section_system,
            symmetry: section_symmetry,
            encyclopedia: EncyclopediaMetadata) -> Symmetry:
        """Returns a populated Symmetry subsection."""
        result = Symmetry()
        filled = False

        if symmetry:
            result.hall_number = symmetry.hall_number
            result.hall_symbol = symmetry.hall_symbol
            result.bravais_lattice = symmetry.bravais_lattice
            result.crystal_system = symmetry.crystal_system
            result.space_group_number = symmetry.space_group_number
            result.space_group_symbol = symmetry.international_short_symbol
            result.point_group = symmetry.point_group
            filled = True

        if encyclopedia and encyclopedia.material.bulk:
            result.strukturbericht_designation = encyclopedia.material.bulk.strukturbericht_designation
            result.structure_name = encyclopedia.material.bulk.structure_type
            result.prototype_formula = encyclopedia.material.bulk.structure_prototype
            filled = True

        proto = repr_sys.section_prototype if repr_sys else None
        proto = proto[0] if proto else None
        if proto:
            result.prototype_aflow_id = proto.prototype_aflow_id
            filled = True

        if filled:
            return result
        return None

    def structure_original(self, repr_sys: section_system) -> Structure:
        """Returns a populated Structure subsection for the original
        structure.
        """
        if repr_sys:
            struct = Structure()
            struct.cartesian_site_positions = repr_sys.atom_positions
            struct.species_at_sites = repr_sys.atom_labels
            self.species(struct.species_at_sites, struct)
            if atomutils.is_valid_basis(repr_sys.lattice_vectors):
                struct.dimension_types = np.array(repr_sys.configuration_periodic_dimensions).astype(int)
                struct.lattice_vectors = repr_sys.lattice_vectors
                struct.cell_volume = atomutils.get_volume(repr_sys.lattice_vectors.magnitude),
                param_values = atomutils.cell_to_cellpar(repr_sys.lattice_vectors.magnitude)
                params = LatticeParameters()
                params.a = float(param_values[0])
                params.b = float(param_values[1])
                params.c = float(param_values[2])
                params.alpha = float(param_values[3])
                params.beta = float(param_values[4])
                params.gamma = float(param_values[5])
                struct.lattice_parameters = params
            return struct

        return None

    def structure_primitive(self, symmetry: section_symmetry) -> Structure:
        """Returns a populated Structure subsection for the primitive
        structure.
        """
        if symmetry:
            struct = Structure()
            prim_sys = symmetry.section_primitive_system[0]
            struct.cartesian_site_positions = prim_sys.atom_positions_primitive
            struct.species_at_sites = atomutils.chemical_symbols(prim_sys.atomic_numbers_primitive)
            self.species(struct.species_at_sites, struct)
            if atomutils.is_valid_basis(prim_sys.lattice_vectors_primitive):
                struct.dimension_types = [1, 1, 1]
                struct.lattice_vectors = prim_sys.lattice_vectors_primitive
                struct.cell_volume = atomutils.get_volume(prim_sys.lattice_vectors_primitive.magnitude),
                param_values = atomutils.cell_to_cellpar(prim_sys.lattice_vectors_primitive.magnitude)
                params = LatticeParameters()
                params.a = float(param_values[0])
                params.b = float(param_values[1])
                params.c = float(param_values[2])
                params.alpha = float(param_values[3])
                params.beta = float(param_values[4])
                params.gamma = float(param_values[5])
                struct.lattice_parameters = params
            return struct

        return None

    def structure_conventional(self, symmetry: section_symmetry) -> Structure:
        """Returns a populated Structure subsection for the conventional
        structure.
        """
        if symmetry:
            struct = Structure()
            conv_sys = symmetry.section_std_system[0]
            struct.cartesian_site_positions = conv_sys.atom_positions_std
            struct.species_at_sites = atomutils.chemical_symbols(conv_sys.atomic_numbers_std)
            self.species(struct.species_at_sites, struct)
            if atomutils.is_valid_basis(conv_sys.lattice_vectors_std):
                struct.dimension_types = [1, 1, 1]
                struct.lattice_vectors = conv_sys.lattice_vectors_std
                struct.cell_volume = atomutils.get_volume(conv_sys.lattice_vectors_std.magnitude),
                param_values = atomutils.cell_to_cellpar(conv_sys.lattice_vectors_std.magnitude)
                params = LatticeParameters()
                params.a = float(param_values[0])
                params.b = float(param_values[1])
                params.c = float(param_values[2])
                params.alpha = float(param_values[3])
                params.beta = float(param_values[4])
                params.gamma = float(param_values[5])
                struct.lattice_parameters = params
                analyzer = symmetry.m_cache["symmetry_analyzer"]
                sets = analyzer.get_wyckoff_sets_conventional(return_parameters=True)
                self.wyckoff_sets(struct, sets)
            return struct

        return None

    def wyckoff_sets(self, struct: Structure, wyckoff_sets: Dict) -> None:
        """Populates the Wyckoff sets in the given structure.
        """
        for group in wyckoff_sets:
            wset = struct.m_create(WyckoffSet)
            if group.x is not None or group.y is not None or group.z is not None:
                if group.x is not None:
                    wset.x = float(group.x)
                if group.y is not None:
                    wset.y = float(group.y)
                if group.z is not None:
                    wset.z = float(group.z)
            wset.indices = group.indices
            wset.element = group.element
            wset.wyckoff_letter = group.wyckoff_letter

    def species(self, labels: List[str], struct: Structure) -> None:
        """Given a list of species labels, returns the
        corresponding Species instance.
        """
        unique_labels = sorted(list(set(labels)))
        for label in unique_labels:
            i_species = struct.m_create(Species)
            i_species.name = label
            i_species.chemical_elements = [label]
            i_species.concentration = [1.0]
