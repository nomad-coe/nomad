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
from typing import Dict, List, Union

from nomad import config
from nomad import atomutils
from nomad.normalizing.normalizer import Normalizer
from nomad.datamodel.encyclopedia import EncyclopediaMetadata
from nomad.datamodel.optimade import OptimadeEntry, Species
from nomad.datamodel.metainfo.common_dft import section_symmetry, section_system, GeometryOptimization
from nomad.datamodel.results import (
    Results,
    Material,
    Method,
    Properties,
    Symmetry,
    Structure,
    StructureOriginal,
    StructurePrimitive,
    StructureConventional,
    LatticeParameters,
    WyckoffSet,
    Simulation,
    DFT,
    GW,
    xc_treatments,
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

        symmetry = None
        if repr_sys and repr_sys.section_symmetry:
            symmetry = repr_sys.section_symmetry[0]

        # Create the section and populate the subsections
        results = Results()
        results.material = self.material(repr_sys, symmetry, encyclopedia, optimade)
        results.method = self.method(encyclopedia)
        results.properties = self.properties(repr_sys, symmetry, encyclopedia)
        self.entry_archive.results = results

    def material(
            self,
            repr_sys: section_system,
            symmetry: section_symmetry,
            encyclopedia: EncyclopediaMetadata,
            optimade: OptimadeEntry) -> Material:
        """Returns a populated Material subsection."""
        material = Material()

        if repr_sys:
            material.type_structural = repr_sys.system_type
            names, counts = atomutils.get_hill_decomposition(repr_sys.atom_labels, reduced=True)
            material.chemical_formula_reduced_fragments = [
                "{}{}".format(n, int(c) if c != 1 else "") for n, c in zip(names, counts)
            ]
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

        if symmetry:
            symm = self.symmetry(repr_sys, symmetry, encyclopedia)
            if symm:
                material.symmetry = symm

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

    def wyckoff_sets(self, struct: StructureConventional, wyckoff_sets: Dict) -> None:
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
        """Given a list of species labels, creates the corresponding Species
        sections in the given structure.
        """
        if labels is None:
            return
        unique_labels = sorted(list(set(labels)))
        for label in unique_labels:
            i_species = struct.m_create(Species)
            i_species.name = label
            i_species.chemical_elements = [label]
            i_species.concentration = [1.0]

    def basis_set_type(self) -> Union[str, None]:
        name = self.section_run.program_basis_set_type
        if name:
            key = name.replace('_', '').replace('-', '').replace(' ', '').lower()
            name_mapping = {
                'gaussians': 'gaussians',
                'realspacegrid': 'real-space grid',
                'planewaves': 'plane waves'
            }
            name = name_mapping.get(key, name)
        return name

    def core_electron_treatment(self) -> str:
        treatment = config.services.unavailable_value
        code_name = self.section_run.program_name
        if code_name is not None:
            core_electron_treatments = {
                'VASP': 'pseudopotential',
                'FHI-aims': 'full all electron',
                'exciting': 'full all electron',
                'quantum espresso': 'pseudopotential'
            }
            treatment = core_electron_treatments.get(code_name, config.services.unavailable_value)
        return treatment

    def xc_functional_names(self, repr_method) -> Union[List[str], None]:
        if repr_method:
            functionals = repr_method.section_XC_functionals
            if functionals:
                names = []
                for functional in functionals:
                    name = functional.XC_functional_name
                    if name:
                        names.append(name)
                return sorted(names)
        return None

    def xc_functional_type(self, xc_functionals) -> str:
        if xc_functionals:
            name = xc_functionals[0]
            return xc_treatments.get(name[:3].lower(), config.services.unavailable_value)
        else:
            return config.services.unavailable_value

    def method(
            self,
            encyclopedia: EncyclopediaMetadata) -> Method:
        """Returns a populated Method subsection."""
        method = Method()
        simulation = Simulation()
        repr_method = None
        method_name = config.services.unavailable_value
        methods = self.section_run.section_method
        n_methods = len(methods)

        if n_methods == 1:
            repr_method = methods[0]
            method_name = repr_method.electronic_structure_method
            if method_name is None:
                method_name = config.services.unavailable_value
        elif n_methods > 1:
            for sec_method in methods:
                # GW
                electronic_structure_method = sec_method.electronic_structure_method
                if electronic_structure_method in {"G0W0", "scGW"}:
                    repr_method = sec_method
                    method_name = electronic_structure_method
                    break

                # Methods linked to each other through references. Get all
                # linked methods, try to get electronic_structure_method from
                # each.
                try:
                    refs = sec_method.section_method_to_method_refs
                except KeyError:
                    pass
                else:
                    linked_methods = [sec_method]
                    for ref in refs:
                        method_to_method_kind = ref.method_to_method_kind
                        method_to_method_ref = ref.method_to_method_ref
                        if method_to_method_kind == "core_settings":
                            linked_methods.append(method_to_method_ref)

                    for i_method in linked_methods:
                        electronic_structure_method = i_method.electronic_structure_method
                        if electronic_structure_method is not None:
                            repr_method = sec_method
                            method_name = electronic_structure_method

        if method_name in {"G0W0", "scGW"}:
            method.method_name = "GW"
            gw = GW()
            gw.gw_type = repr_method.gw_type
            gw.starting_point = repr_method.gw_starting_point.split()
            simulation.gw = gw
        elif method_name in {"DFT", "DFT+U"}:
            method.method_name = "DFT"
            dft = DFT()
            dft.basis_set_type = self.basis_set_type()
            dft.core_electron_treatment = self.core_electron_treatment()
            dft.smearing_kind = repr_method.smearing_kind
            dft.smearing_width = repr_method.smearing_width
            if repr_method.number_of_spin_channels:
                dft.spin_polarized = repr_method.number_of_spin_channels > 1
            dft.xc_functional_names = self.xc_functional_names(repr_method)
            dft.xc_functional_type = self.xc_functional_type(dft.xc_functional_names)
            dft.scf_threshold_energy_change = repr_method.scf_threshold_energy_change
            dft.van_der_Waals_method = repr_method.van_der_Waals_method
            dft.relativity_method = repr_method.relativity_method
            simulation.dft = dft

        if encyclopedia and encyclopedia.method:
            method.method_id = encyclopedia.method.method_id
        simulation.program_name = self.section_run.program_name
        simulation.program_version = self.section_run.program_version
        method.simulation = simulation

        return method

    def geometry_optimization(self) -> Union[GeometryOptimization, None]:
        """Returns a reference to a workflow result for geometry optimization.
        """
        workflow = self.entry_archive.section_workflow
        if workflow is not None:
            geo_opt = workflow.section_geometry_optimization
            if geo_opt:
                return geo_opt
        return None

    def properties(
            self,
            repr_sys: section_system,
            symmetry: section_symmetry,
            encyclopedia: EncyclopediaMetadata) -> Properties:
        """Returns a populated Properties subsection."""
        properties = Properties()

        # Structures
        struct_orig = self.structure_original(repr_sys)
        if struct_orig:
            properties.structure_original = struct_orig
        struct_prim = self.structure_primitive(symmetry)
        if struct_prim:
            properties.structure_primitive = struct_prim
        struct_conv = self.structure_conventional(symmetry)
        if struct_conv:
            properties.structure_conventional = struct_conv

        geo_opt = self.geometry_optimization()
        if geo_opt:
            properties.geometry_optimization

        return properties

    def structure_original(self, repr_sys: section_system) -> StructureOriginal:
        """Returns a populated Structure subsection for the original
        structure.
        """
        if repr_sys:
            struct = StructureOriginal()
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

    def structure_primitive(self, symmetry: section_symmetry) -> StructurePrimitive:
        """Returns a populated Structure subsection for the primitive
        structure.
        """
        if symmetry:
            struct = StructurePrimitive()
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

    def structure_conventional(self, symmetry: section_symmetry) -> StructureConventional:
        """Returns a populated Structure subsection for the conventional
        structure.
        """
        if symmetry:
            struct = StructureConventional()
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
