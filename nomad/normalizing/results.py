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
import re
import numpy as np
from typing import Dict, List, Union, Any
import ase.data

from nomad import config
from nomad import atomutils
from nomad.normalizing.normalizer import Normalizer
from nomad.datamodel.encyclopedia import EncyclopediaMetadata
from nomad.datamodel.optimade import OptimadeEntry, Species
from nomad.datamodel.metainfo.common_dft import (
    section_symmetry,
    section_system,
)
from nomad.datamodel.results import (
    Results,
    Material,
    Method,
    GeometryOptimizationProperties,
    GeometryOptimizationMethod,
    Properties,
    Symmetry,
    Structures,
    Structure,
    StructureOriginal,
    StructurePrimitive,
    StructureConventional,
    StructureOptimized,
    LatticeParameters,
    WyckoffSet,
    Simulation,
    DFT,
    GW,
    xc_treatments,
    ElectronicProperties,
    VibrationalProperties,
    BandStructureElectronic,
    BandStructurePhonon,
    DOSElectronic,
    DOSPhonon,
    EnergyFreeHelmholtz,
    HeatCapacityConstantVolume,
)

re_label = re.compile("^([a-zA-Z][a-zA-Z]?)[^a-zA-Z]*")
elements = set(ase.data.chemical_symbols)


def label_to_chemical_symbol(label: str) -> Union[str, None]:
    """Tries to extract a valid chemical symbol from a label. Currently can
    only handle labels that correspond to chemical names (no matter what case)
    and that possbily have a non-alphabetic postfix. Raises an error if no
    match can be made.
    """
    match = re_label.match(label)
    symbol = None
    if match:
        test = match.group(1).capitalize()
        if test in elements:
            symbol = test
    if symbol is None:
        raise ValueError(
            "Could not identify a chemical element from the given label: {}".format(label)
        )
    return symbol


def valid_array(array: Any) -> bool:
    """Checks if the given variable is a non-empty array.
    """
    return array is not None and len(array) > 0


def isint(value: Any) -> bool:
    """Checks if the given variable can be interpreted as an integer.
    """
    try:
        int(value)
        return True
    except ValueError:
        return False


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
        self.entry_archive.results = results
        results.material = self.material(repr_sys, symmetry, encyclopedia, optimade)
        results.method = self.method(encyclopedia)
        results.properties = self.properties(repr_sys, symmetry, encyclopedia)
        self.geometry_optimization(results.method, results.properties)

        # Add the present quantities. The full path will be saved for each
        # property, and if the leaf quantity/section name is unambiguous, also
        # it will be saved. Each repeating section will be investigated as well
        # to find all present properties.
        available_properties = set()
        for section, m_def, _ in results.properties.m_traverse():
            parent_path = section.m_path()
            path = "{}/{}".format(parent_path if parent_path != "/" else "", m_def.name)[20:]
            parts = filter(lambda x: not isint(x), path.split("/"))
            path = ".".join(parts)
            available_properties.add(path)
        shorthand_prop = list()
        for prop in available_properties:
            name = prop.rsplit(".", 1)[-1]
            shorthand_prop.append(name)
        u, c = np.unique(shorthand_prop, return_counts=True)
        shorthand_prop = u[c == 1]
        available_properties |= set([str(x) for x in shorthand_prop])
        results.properties.available_properties = list(available_properties)

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
            try:
                symbol = label_to_chemical_symbol(label)
            except ValueError:
                self.logger.info("could not identify chemical symbol from the label: {}".format(label))
            else:
                i_species.chemical_symbols = [symbol]
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

    def band_structure_electronic(self) -> Union[BandStructureElectronic, None]:
        """Returns a new section containing an electronic band structure. In
        the case of multiple valid band structures, only the latest one is
        considered.

       Band structure is reported only under the following conditions:
          - There is a non-empty array of band_k_points.
          - There is a non-empty array of band_energies.
          - The reported band_structure_kind is not "vibrational".
        """
        path = ["section_run", "section_single_configuration_calculation", "section_k_band"]
        for bs in self.traverse_reversed(path):
            kind = bs.band_structure_kind
            if kind == "vibrational" or not bs.section_k_band_segment:
                continue
            valid = True
            for segment in bs.section_k_band_segment:
                energies = segment.band_energies
                k_points = segment.band_k_points
                if not valid_array(energies) or not valid_array(k_points):
                    valid = False
                    break
            if valid:
                # Fill band structure data to the newer, improved data layout
                bs_new = BandStructureElectronic()
                bs_new.reciprocal_cell = bs
                bs_new.segments = bs.section_k_band_segment
                bs_new.spin_polarized = bs_new.segments[0].band_energies.shape[0] > 1
                bs_new.energy_fermi = bs.m_parent.energy_reference_fermi
                bs_new.energy_references = bs.energy_references
                if bs.section_band_gap:
                    energy_highest_occupied = []
                    energy_lowest_unoccupied = []
                    band_gap = []
                    band_gap_type = []
                    for gap in bs.section_band_gap:
                        band_gap.append(gap.value)
                        band_gap_type.append("no_gap" if gap.type is None else gap.type)
                        energy_highest_occupied.append(gap.valence_band_max_energy)
                        energy_lowest_unoccupied.append(gap.conduction_band_min_energy)
                    bs_new.band_gap = band_gap
                    bs_new.band_gap_type = band_gap_type
                    if bs.m_parent.energy_reference_highest_occupied is not None:
                        bs_new.energy_highest_occupied = bs.m_parent.energy_reference_highest_occupied
                    else:
                        bs_new.energy_highest_occupied = energy_highest_occupied
                    if bs.m_parent.energy_reference_lowest_unoccupied is not None:
                        bs_new.energy_lowest_unoccupied = bs.m_parent.energy_reference_lowest_unoccupied
                    else:
                        bs_new.energy_lowest_unoccupied = energy_lowest_unoccupied
                return bs_new

        return None

    def dos_electronic(self) -> Union[DOSElectronic, None]:
        """Returns a reference to the section containing an electronic dos. In
        the case of multiple valid DOSes, only the latest one is reported.

       DOS is reported only under the following conditions:
          - There is a non-empty array of dos_values_normalized.
          - There is a non-empty array of dos_energies.
          - The reported dos_kind is not "vibrational".
        """
        path = ["section_run", "section_single_configuration_calculation", "section_dos"]
        for dos in self.traverse_reversed(path):
            kind = dos.dos_kind
            energies = dos.dos_energies
            values = dos.dos_values_normalized
            if kind != "vibrational" and valid_array(energies) and valid_array(values):
                dos_new = DOSElectronic()
                dos_new.energies = dos
                dos_new.densities = dos
                n_channels = values.shape[0]
                dos_new.spin_polarized = n_channels > 1
                dos_new.energy_references = dos.energy_references
                return dos_new

        return None

    def band_structure_phonon(self) -> Union[BandStructurePhonon, None]:
        """Returns a new section containing a phonon band structure. In
        the case of multiple valid band structures, only the latest one is
        considered.

       Band structure is reported only under the following conditions:
          - There is a non-empty array of band_k_points.
          - There is a non-empty array of band_energies.
          - The reported band_structure_kind is "vibrational".
        """
        path = ["section_run", "section_single_configuration_calculation", "section_k_band"]
        for bs in self.traverse_reversed(path):
            kind = bs.band_structure_kind
            if kind != "vibrational" or not bs.section_k_band_segment:
                continue
            valid = True
            for segment in bs.section_k_band_segment:
                energies = segment.band_energies
                k_points = segment.band_k_points
                if not valid_array(energies) or not valid_array(k_points):
                    valid = False
                    break
            if valid:
                # Fill band structure data to the newer, improved data layout
                bs_new = BandStructurePhonon()
                bs_new.segments = bs.section_k_band_segment
                return bs_new

        return None

    def dos_phonon(self) -> Union[DOSPhonon, None]:
        """Returns a section containing phonon dos data. In the case of
        multiple valid data sources, only the latest one is reported.

       DOS is reported only under the following conditions:
          - There is a non-empty array of dos_values_normalized.
          - There is a non-empty array of dos_energies.
          - The reported dos_kind is "vibrational".
        """
        path = ["section_run", "section_single_configuration_calculation", "section_dos"]
        for dos in self.traverse_reversed(path):
            kind = dos.dos_kind
            energies = dos.dos_energies
            values = dos.dos_values_normalized
            if kind == "vibrational" and valid_array(energies) and valid_array(values):
                # Fill dos data to the newer, improved data layout
                dos_new = DOSPhonon()
                dos_new.energies = dos
                dos_new.densities = dos
                return dos_new

        return None

    def energy_free_helmholtz(self) -> Union[EnergyFreeHelmholtz, None]:
        """Returns a section Helmholtz free energy data. In the case of
        multiple valid data sources, only the latest one is reported.

       Helmholtz free energy is reported only under the following conditions:
          - There is a non-empty array of temperatures.
          - There is a non-empty array of energies.
        """
        path = ["section_run", "section_frame_sequence", "section_thermodynamical_properties"]
        for thermo_prop in self.traverse_reversed(path):
            temperatures = thermo_prop.thermodynamical_property_temperature
            energies = thermo_prop.vibrational_free_energy_at_constant_volume
            if valid_array(temperatures) and valid_array(energies):
                energy_free = EnergyFreeHelmholtz()
                energy_free.energies = thermo_prop
                energy_free.temperatures = thermo_prop
                return energy_free

        return None

    def heat_capacity_constant_volume(self) -> Union[HeatCapacityConstantVolume, None]:
        """Returns a section containing heat capacity data. In the case of
        multiple valid data sources, only the latest one is reported.

       Heat capacity is reported only under the following conditions:
          - There is a non-empty array of temperatures.
          - There is a non-empty array of energies.
        """
        path = ["section_run", "section_frame_sequence", "section_thermodynamical_properties"]
        for thermo_prop in self.traverse_reversed(path):
            temperatures = thermo_prop.thermodynamical_property_temperature
            heat_capacities = thermo_prop.specific_heat_capacity
            if valid_array(temperatures) and valid_array(heat_capacities):
                heat_cap = HeatCapacityConstantVolume()
                heat_cap.heat_capacities = thermo_prop
                heat_cap.temperatures = thermo_prop
                return heat_cap

        return None

    def geometry_optimization(self, method: Method, properties: Properties) -> None:
        """Populates both geometry optimization methodology and calculated
        properties based on the first found geometry optimization workflow.
        """
        path = ["section_workflow"]
        for workflow in self.traverse_reversed(path):
            # Check validity
            if workflow.workflow_type == "geometry_optimization" and workflow.calculations_ref:

                # Method
                geo_opt_wf = workflow.section_geometry_optimization
                if geo_opt_wf is not None:
                    geo_opt_meth = GeometryOptimizationMethod()
                    geo_opt_meth.geometry_optimization_type = geo_opt_wf.geometry_optimization_type
                    method.simulation.geometry_optimization = geo_opt_meth

                # Properties
                geo_opt_prop = GeometryOptimizationProperties()
                geo_opt_prop.trajectory = workflow.calculations_ref
                structure_optimized = self.structure_optimized(
                    workflow.calculation_result_ref.single_configuration_calculation_to_system_ref
                )
                if structure_optimized:
                    geo_opt_prop.structure_optimized = structure_optimized
                if geo_opt_wf is not None:
                    if geo_opt_wf.energies is not None:
                        geo_opt_prop.energies = geo_opt_wf
                    geo_opt_prop.final_energy_difference = geo_opt_wf.final_energy_difference
                    geo_opt_prop.final_force_maximum = geo_opt_wf.final_force_maximum
                properties.geometry_optimization = geo_opt_prop

                return

    def properties(
            self,
            repr_sys: section_system,
            symmetry: section_symmetry,
            encyclopedia: EncyclopediaMetadata) -> Properties:
        """Returns a populated Properties subsection."""
        properties = Properties()

        # Structures
        struct_orig = self.structure_original(repr_sys)
        struct_prim = self.structure_primitive(symmetry)
        struct_conv = self.structure_conventional(symmetry)
        if struct_orig or struct_prim or struct_conv:
            structures = Structures()
            if struct_conv:
                structures.structure_conventional = struct_conv
            if struct_prim:
                structures.structure_primitive = struct_prim
            if struct_orig:
                structures.structure_original = struct_orig
            properties.structures = structures

        # Electronic
        bs_electronic = self.band_structure_electronic()
        dos_electronic = self.dos_electronic()
        if bs_electronic or dos_electronic:
            electronic = ElectronicProperties()
            if bs_electronic:
                electronic.band_structure_electronic = bs_electronic
            if dos_electronic:
                electronic.dos_electronic = dos_electronic
            properties.electronic = electronic

        # Vibrational
        bs_phonon = self.band_structure_phonon()
        dos_phonon = self.dos_phonon()
        energy_free = self.energy_free_helmholtz()
        heat_cap = self.heat_capacity_constant_volume()
        if bs_phonon or dos_phonon or energy_free or heat_cap:
            vibrational = VibrationalProperties()
            if dos_phonon:
                vibrational.dos_phonon = dos_phonon
            if bs_phonon:
                vibrational.band_structure_phonon = bs_phonon
            if energy_free:
                vibrational.energy_free_helmholtz = energy_free
            if heat_cap:
                vibrational.heat_capacity_constant_volume = heat_cap
            properties.vibrational = vibrational

        try:
            n_calc = len(self.section_run.section_single_configuration_calculation)
        except Exception:
            n_calc = 0
        properties.n_calculations = n_calc

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
                struct.cell_volume = atomutils.get_volume(repr_sys.lattice_vectors.magnitude)
                struct.lattice_parameters = self.lattice_parameters(repr_sys.lattice_vectors)
            return struct

        return None

    def structure_primitive(self, symmetry: section_symmetry) -> StructurePrimitive:
        """Returns a populated Structure subsection for the primitive
        structure.
        """
        if symmetry:
            struct = StructurePrimitive()
            prim_sys = symmetry.section_primitive_system[0]
            struct.species_at_sites = atomutils.chemical_symbols(prim_sys.atomic_numbers_primitive)
            self.species(struct.species_at_sites, struct)
            lattice_vectors = prim_sys.lattice_vectors_primitive
            if atomutils.is_valid_basis(lattice_vectors):
                struct.cartesian_site_positions = atomutils.to_cartesian(prim_sys.atom_positions_primitive, lattice_vectors)
                struct.dimension_types = [1, 1, 1]
                struct.lattice_vectors = lattice_vectors
                struct.cell_volume = atomutils.get_volume(lattice_vectors.magnitude)
                struct.lattice_parameters = self.lattice_parameters(lattice_vectors)
            return struct

        return None

    def structure_conventional(self, symmetry: section_symmetry) -> StructureConventional:
        """Returns a populated Structure subsection for the conventional
        structure.
        """
        if symmetry:
            struct = StructureConventional()
            conv_sys = symmetry.section_std_system[0]
            struct.species_at_sites = atomutils.chemical_symbols(conv_sys.atomic_numbers_std)
            self.species(struct.species_at_sites, struct)
            lattice_vectors = conv_sys.lattice_vectors_std
            if atomutils.is_valid_basis(lattice_vectors):
                struct.cartesian_site_positions = atomutils.to_cartesian(conv_sys.atom_positions_std, lattice_vectors)
                struct.dimension_types = [1, 1, 1]
                struct.lattice_vectors = lattice_vectors
                struct.cell_volume = atomutils.get_volume(lattice_vectors.magnitude)
                struct.lattice_parameters = self.lattice_parameters(lattice_vectors)
                analyzer = symmetry.m_cache["symmetry_analyzer"]
                sets = analyzer.get_wyckoff_sets_conventional(return_parameters=True)
                self.wyckoff_sets(struct, sets)
            return struct

        return None

    def structure_optimized(self, system: section_system) -> Structure:
        """Returns a populated Structure subsection for the optimized
        structure.
        """
        if system:
            struct = StructureOptimized()
            struct.cartesian_site_positions = system.atom_positions
            struct.species_at_sites = system.atom_labels
            self.species(struct.species_at_sites, struct)
            if atomutils.is_valid_basis(system.lattice_vectors):
                struct.dimension_types = np.array(system.configuration_periodic_dimensions).astype(int)
                struct.lattice_vectors = system.lattice_vectors
                struct.cell_volume = atomutils.get_volume(system.lattice_vectors.magnitude)
                struct.lattice_parameters = self.lattice_parameters(system.lattice_vectors)
            return struct

        return None

    def lattice_parameters(self, lattice_vectors) -> LatticeParameters:
        """Converts the given cell into LatticeParameters"""
        param_values = atomutils.cell_to_cellpar(lattice_vectors.magnitude)
        params = LatticeParameters()
        params.a = float(param_values[0])
        params.b = float(param_values[1])
        params.c = float(param_values[2])
        params.alpha = float(param_values[3])
        params.beta = float(param_values[4])
        params.gamma = float(param_values[5])
        return params

    def traverse_reversed(self, path: List[str]) -> Any:
        """Traverses the given metainfo path in reverse order. Useful in
        finding the latest reported section or value.
        """
        def traverse(root, path, i):
            sections = getattr(root, path[i])
            if isinstance(sections, list):
                for section in reversed(sections):
                    if i == len(path) - 1:
                        yield section
                    else:
                        for s in traverse(section, path, i + 1):
                            yield s
            else:
                if i == len(path) - 1:
                    yield sections
                else:
                    for s in traverse(sections, path, i + 1):
                        yield s
        for t in traverse(self.entry_archive, path, 0):
            if t is not None:
                yield t
