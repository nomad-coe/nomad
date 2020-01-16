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

from hashlib import sha512
from typing import Dict
from abc import abstractmethod
import math
import ase
import numpy as np

from nomad.normalizing.normalizer import Normalizer, s_scc, s_system, s_frame_sequence, r_frame_sequence_to_sampling, s_sampling_method, r_frame_sequence_local_frames, r_scc_to_system
from nomad.metainfo.encyclopedia import Encyclopedia, Material, Calculation
from nomad.normalizing import structure
from nomad import config


class EncyclopediaNormalizer(Normalizer):
    """
    This normalizer emulates the functionality of the old Encyclopedia backend.
    The data used by the encyclopedia have been assigned under new metainfo
    within section_encyclopedia. In the future these separate metainfos could
    be absorbed into the existing metainfo hiearchy.
    """
    def __init__(self, backend):
        super().__init__(backend)

    # NOTE: Band structure normalizer
    def band_gap(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def band_gap_position(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def band_gap_type(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def band_structure(self) -> None:
        pass

    # NOTE: Method normalizer
    def basis_set_short_name(self) -> None:
        pass

    # NOTE: Method normalizer
    def basis_set_type(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def brillouin_zone(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def brillouin_zone_viewer(self) -> None:
        pass

    def calculation(self) -> None:
        pass

    def calculation_pid(self) -> None:
        pass

    # NOTE: Parser
    def code_name(self) -> None:
        pass

    # NOTE: Parser
    def code_version(self) -> None:
        pass

    # NOTE: Repo
    def contributor_first_name(self) -> None:
        pass

    # NOTE: Repo
    def contributor_last_name(self) -> None:
        pass

    # NOTE: Repo
    def contributor_type(self) -> None:
        pass

    # NOTE: Repo
    def contributors(self) -> None:
        pass

    # NOTE: Method normalizer
    def core_electron_treatment(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def dos(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def elastic_constants_matrix(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def elastic_deformation_energies(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def elastic_fitting_parameters(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def elastic_moduli(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def elastic_properties(self) -> None:
        pass

    def energies(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def fermi_surface(self) -> None:
        pass

    # NOTE: Method normalizer
    def functional_long_name(self) -> None:
        pass

    # NOTE: Method normalizer
    def functional_type(self) -> None:
        pass

    # TODO: ??
    def group_e_min(self) -> None:
        pass

    # TODO: ??
    def group_type(self) -> None:
        pass

    # TODO: Method normalizer
    def gw_starting_point(self) -> None:
        pass

    # TODO: Method normalizer
    def gw_type(self) -> None:
        pass

    # NOTE: Enc specific
    def has_bs(self) -> None:
        pass

    # NOTE: Enc specific
    def has_dos(self) -> None:
        pass

    # NOTE: Enc specific
    def has_fermi_surface(self) -> None:
        pass

    # NOTE: Enc specific
    def has_thermal_properties(self) -> None:
        pass

    def helmholtz_free_energy(self) -> None:
        pass

    def k_point_grid_description(self) -> None:
        pass

    def mainfile_uri(self) -> None:
        pass

    # NOTE: Postprocessing
    def number_of_calculation(self) -> None:
        pass

    # NOTE: Phonon normalizer
    def phonon_dispersion(self) -> None:
        pass

    # NOTE: Phonon normalizer
    def phonon_dos(self) -> None:
        pass

    # NOTE: Method normalizer
    def pseudopotential_type(self) -> None:
        pass

    # NOTE: Repo
    def repository_dowload_uri(self) -> None:
        pass

    # NOTE: Repo
    def repository_upload_comment(self) -> None:
        pass

    # NOTE: Repo
    def repository_uri(self) -> None:
        pass

    # NOTE: Enc specific
    def run_type(self, calculation) -> str:
        """Decides what type of calculation this is: single_point, md,
        geometry_optimization, etc.
        """
        run_enums = Calculation.run_type.type
        run_type = run_enums.unavailable

        try:
            sccs = self._backend[s_scc]
        except Exception:
            sccs = []
        try:
            frame_sequences = self._backend[s_frame_sequence]
        except Exception:
            frame_sequences = []

        n_scc = len(sccs)
        n_frame_seq = len(frame_sequences)

        # Only one system, no sequences
        if n_scc == 1 and n_frame_seq == 0:
            program_name = self._backend["program_name"]
            if program_name == "elastic":
                # TODO move to taylor expansion as soon as data is correct in archive
                run_type = run_enums.elastic_constants
            else:
                run_type = run_enums.single_point
        # One sequence. Currently calculations with multiple sequences are
        # unsupported.
        elif n_frame_seq == 1:
            frame_seq = frame_sequences[0]

            # See if sampling_method is present
            try:
                i_sampling_method = frame_seq[r_frame_sequence_to_sampling]
            except KeyError:
                self.logger.info(
                    "Cannot determine encyclopedia run type because missing "
                    "value for frame_sequence_to_sampling_ref"
                )
                return run_type

            # See if local frames are present
            try:
                frames = frame_seq[r_frame_sequence_local_frames]
            except KeyError:
                self.logger.info(
                    "section_frame_sequence_local_frames not found although a "
                    "frame_sequence exists"
                )
                return run_type
            if len(frames) == 0:
                self.logger.info("No frames referenced in section_frame_sequence_local_frames")
                return run_type

            section_sampling_method = self._backend[s_sampling_method][i_sampling_method]
            sampling_method = section_sampling_method["sampling_method"]

            if sampling_method == "molecular_dynamics":
                run_type = run_enums.molecular_dynamics
            if sampling_method == "geometry_optimization":
                run_type = run_enums.geometry_optimization
            if sampling_method == "taylor_expansion":
                run_type = run_enums.phonon_calculation

        calculation.run_type = run_type
        return run_type

    def scf_threshold(self) -> None:
        pass

    # NOTE: Enc specific
    def similar_materials(self) -> None:
        pass

    # NOTE: Method normalizer
    def smearing_kind(self) -> None:
        pass

    # NOTE: Method normalizer
    def smearing_parameters(self) -> None:
        pass

    # NOTE: System normalizer
    def space_group(self) -> None:
        pass

    # NOTE: System normalizer
    def space_group_international_short_symbol(self) -> None:
        pass

    # NOTE: System normalizer
    def space_group_number(self) -> None:
        pass

    # NOTE: Phonon normalizer
    def specific_heat_cv(self) -> None:
        pass

    # NOTE: System normalizer
    def springer_classification(self) -> None:
        pass

    # NOTE: System normalizer
    def springer_compound_class(self) -> None:
        pass

    # NOTE: System normalizer
    def springer_prototype(self) -> None:
        pass

    # NOTE: System normalizer
    def structure_prototype(self) -> None:
        pass

    # NOTE: System normalizer
    def structure_type(self) -> None:
        pass

    # NOTE: System normalizer
    def strukturbericht_designation(self) -> None:
        pass

    # NOTE: System normalizer
    def system_type(self, material: Material, calculation: Calculation) -> tuple:
        # Select the representative system from which system type is retrieved.
        # For geometry optimizations system type is analyzed from last relaxed
        # frame. For phonon calculations system type is analyzed from first
        # undistorted frame. For molecular dynamics system type is analyzed
        # from first frame
        system_type = config.services.unavailable_value
        system = None
        run_enums = Calculation.run_type.type
        system_enums = Material.system_type.type
        if calculation.run_type in {
                run_enums.geometry_optimization,
                run_enums.molecular_dynamics,
                run_enums.phonon_calculation}:
            frame_seqs = self._backend[s_frame_sequence]
            frame_seq = frame_seqs[0]
            frames = frame_seq[r_frame_sequence_local_frames]
            sccs = self._backend[s_scc]
            systems = self._backend[s_system]
            if calculation.run_type == run_enums.geometry_optimization:
                idx = -1
            elif calculation.run_type == run_enums.phonon_calculation:
                idx = 0
            elif calculation.run_type == run_enums.molecular_dynamics:
                idx = 0
            scc = sccs[frames[idx]]
            r_system = scc[r_scc_to_system]
            system = systems[r_system]
        elif calculation.run_type == run_enums.single_point:
            system = self._backend[s_system][0]

        # Try to find system type information from backend for the selected system.
        try:
            stype = system["system_type"]
        except KeyError:
            self.logger.info("System type information not available for encyclopedia")
        else:
            if stype == "2D / surface":
                system_type = system_enums.two_d
            if stype == system_enums.bulk or stype == system_enums.one_d:
                system_type = stype

        material.system_type = system_type
        return system, system_type

    def template(self) -> None:
        pass

    # NOTE: System normalizer
    def wyckoff_groups(self) -> None:
        pass

    def fill(self, run_type, system_type, representative_system):
        # Fill structure related meta
        if system_type == Material.system_type.type.bulk:
            structure_worker = StructureBulk()
            structure_worker.fill(self._backend, representative_system)

    def normalize(self, logger=None) -> None:
        super().normalize(logger)
        system_enums = Material.system_type.type

        # Initialise metainfo structure
        sec_enc = Encyclopedia()
        material = sec_enc.m_create(Material)
        calculation = sec_enc.m_create(Calculation)

        # Determine run type, stop if unknown
        run_type = self.run_type(calculation)
        if run_type == config.services.unavailable_value:
            self.logger.info("unknown run type for encyclopedia")
            return

        # Get the system type, stop if unknown
        representative_system, system_type = self.system_type(material, calculation)
        if system_type != system_enums.bulk and system_type != system_enums.two_d and system_type != system_enums.one_d:
            self.logger.info("unknown system type for encyclopedia")
            return

        # Get the method type, stop if unknown
        # TODO

        # Put the encyclopedia section into backend
        self._backend.add_mi2_section(sec_enc)
        self.fill(run_type, system_type, representative_system)


class Structure():
    """A base class that is used for processing structure related information
    in the Encylopedia.
    """
    @abstractmethod
    def atom_labels(self, material: Material, std_atoms: ase.Atoms) -> None:
        pass

    @abstractmethod
    def atom_positions(self, material: Material, std_atoms: ase.Atoms) -> None:
        pass

    @abstractmethod
    def atomic_density(self, calculation: Calculation, repr_system: ase.Atoms) -> None:
        pass

    @abstractmethod
    def bravais_lattice(self, material: Material, section_system: Dict) -> None:
        pass

    @abstractmethod
    def cell_angles_string(self, calculation: Calculation) -> None:
        pass

    @abstractmethod
    def cell_normalized(self, material: Material, std_atoms: ase.Atoms) -> None:
        pass

    @abstractmethod
    def cell_primitive(self, material: Material, std_atoms: ase.Atoms) -> None:
        pass

    @abstractmethod
    def cell_volume(self, calculation: Calculation, std_atoms: ase.Atoms) -> None:
        pass

    @abstractmethod
    def crystal_system(self, material: Material, section_system: Dict) -> None:
        pass

    @abstractmethod
    def formula(self, material: Material, prim_sys: ase.Atoms) -> None:
        pass

    @abstractmethod
    def formula_reduced(self, material: Material, prim_sys: ase.Atoms) -> None:
        pass

    # def free_wyckoff_parameters(self) -> None:
        # pass

    @abstractmethod
    def lattice_parameters(self, calculation: Calculation, std_atoms: ase.Atoms) -> None:
        pass

    @abstractmethod
    def mass_density(self, calculation: Calculation, repr_system: ase.Atoms) -> None:
        pass

    @abstractmethod
    def material_hash(self, material: Material, section_system: Dict) -> None:
        pass

    # def material_name(self) -> None:
        # pass

    @abstractmethod
    def number_of_atoms(self, material: Material, std_atoms: ase.Atoms) -> None:
        pass

    @abstractmethod
    def periodicity(self, material: Material) -> None:
        pass

    @abstractmethod
    def point_group(self, material: Material, section_symmetry: Dict) -> None:
        pass

    def fill(self, backend, representative_system: Dict) -> None:
        # Fetch resources
        sec_enc = backend.get_mi2_section(Encyclopedia.m_def)
        material = sec_enc.material
        calculation = sec_enc.calculation
        sec_symmetry = representative_system["section_symmetry"][0]
        std_atoms = sec_symmetry["section_std_system"][0].tmp["std_atoms"]  # Temporary value stored by SystemNormalizer
        prim_atoms = sec_symmetry["section_primitive_system"][0].tmp["prim_atoms"]  # Temporary value stored by SystemNormalizer
        repr_atoms = sec_symmetry["section_original_system"][0].tmp["orig_atoms"]  # Temporary value stored by SystemNormalizer

        # Fill structural information
        self.mass_density(calculation, repr_atoms)
        self.material_hash(material, sec_symmetry)
        self.number_of_atoms(material, std_atoms)
        self.atom_labels(material, std_atoms)
        self.atomic_density(calculation, repr_atoms)
        self.bravais_lattice(material, sec_symmetry)
        self.cell_normalized(material, std_atoms)
        self.cell_volume(calculation, std_atoms)
        self.crystal_system(material, sec_symmetry)
        self.cell_primitive(material, prim_atoms)
        self.formula(material, prim_atoms)
        self.formula_reduced(material, prim_atoms)
        self.lattice_parameters(calculation, std_atoms)
        self.cell_angles_string(calculation)
        self.periodicity(material)
        self.point_group(material, sec_symmetry)


class StructureBulk(Structure):
    """Processes structure related metainfo for Encyclopedia bulk structures.
    """
    def atom_labels(self, material: Material, std_atoms: ase.Atoms) -> None:
        material.atom_labels = std_atoms.get_chemical_symbols()

    def atom_positions(self, material: Material, std_atoms: ase.Atoms) -> None:
        material.atom_positions = std_atoms.get_scaled_positions(wrap=False)

    def atomic_density(self, calculation: Calculation, repr_system: ase.Atoms) -> None:
        orig_n_atoms = len(repr_system)
        orig_volume = repr_system.get_volume() * (1e-10)**3
        calculation.atomic_density = float(orig_n_atoms / orig_volume)

    def mass_density(self, calculation: Calculation, repr_system: ase.Atoms) -> None:
        orig_volume = repr_system.get_volume() * (1e-10)**3
        mass = structure.get_summed_atomic_mass(repr_system.get_atomic_numbers())
        calculation.mass_density = float(mass / orig_volume)

    def material_hash(self, material: Material, section_symmetry: Dict) -> None:
        # Get symmetry information from the section
        space_group_number = section_symmetry["space_group_number"]
        section_std_system = section_symmetry["section_std_system"][0]
        wyckoff_sets = section_std_system.tmp["wyckoff_sets"]

        # Create and store hash based on SHA512
        norm_hash_string = structure.create_symmetry_string(space_group_number, wyckoff_sets)
        material.material_hash = sha512(norm_hash_string.encode('utf-8')).hexdigest()

    def number_of_atoms(self, material: Material, std_atoms: ase.Atoms) -> None:
        material.number_of_atoms = len(std_atoms)

    def bravais_lattice(self, material: Material, section_symmetry: Dict) -> None:
        bravais_lattice = section_symmetry["bravais_lattice"]
        material.bravais_lattice = bravais_lattice

    def cell_angles_string(self, calculation: Calculation) -> None:
        angles = calculation.lattice_parameters[3:6]
        angles_rounded = []
        for angle in angles:
            angle_deg = math.degrees(angle)
            angles_rounded.append(
                round(angle_deg / config.normalize.angle_rounding) * config.normalize.angle_rounding)
        calculation.cell_angles_string = "/".join([str(angle) for angle in angles_rounded])

    def cell_normalized(self, material: Material, std_atoms: ase.Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        material.cell_normalized = cell_normalized

    def cell_primitive(self, material: Material, prim_atoms: ase.Atoms) -> None:
        cell_prim = prim_atoms.get_cell()
        cell_prim *= 1e-10
        material.cell_primitive = cell_prim

    def cell_volume(self, calculation: Calculation, std_atoms: ase.Atoms) -> None:
        calculation.cell_volume = float(std_atoms.get_volume() * 1e-10**3)

    def crystal_system(self, material: Material, section_symmetry: Dict) -> None:
        material.crystal_system = section_symmetry["crystal_system"]

    def formula(self, material: Material, prim_sys: ase.Atoms) -> None:
        chemical_symbols = prim_sys.get_chemical_symbols()
        formula = structure.get_hill_formula(chemical_symbols, reduced=False)
        material.formula = formula

    def formula_reduced(self, material: Material, prim_sys: ase.Atoms) -> None:
        chemical_symbols = prim_sys.get_chemical_symbols()
        formula = structure.get_hill_formula(chemical_symbols, reduced=True)
        material.formula_reduced = formula

    def lattice_parameters(self, calculation: Calculation, std_atoms: ase.Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        calculation.lattice_parameters = structure.get_lattice_parameters(cell_normalized)

    def periodicity(self, material: Material) -> None:
        material.periodicity = np.array([0, 1, 2], dtype=np.int8)

    def point_group(self, material: Material, section_symmetry: Dict) -> None:
        point_group = section_symmetry["point_group"]
        material.point_group = point_group
