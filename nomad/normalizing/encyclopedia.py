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
from typing import Dict, List
from math import gcd as gcd
from functools import reduce
from abc import abstractmethod
import json
import ase
import ase.data
from ase import Atoms
import numpy as np
from matid import SymmetryAnalyzer
import matid.geometry

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

    # NOTE: Phonon normalizer
    def specific_heat_cv(self) -> None:
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
            if stype == system_enums.bulk or stype == system_enums.one_d or stype == system_enums.two_d:
                system_type = stype

        material.system_type = system_type
        return system, system_type

    def template(self) -> None:
        pass

    def fill(self, run_type, system_type, representative_system):
        # Fill structure related metainfo
        struct = None
        if system_type == Material.system_type.type.bulk:
            struct = StructureBulk(self._backend, self.logger)
        elif system_type == Material.system_type.type.two_d:
            struct = Structure2D(self._backend, self.logger)
        elif system_type == Material.system_type.type.one_d:
            struct = Structure1D(self._backend, self.logger)
        if struct is not None:
            struct.fill(representative_system)

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
    def __init__(self, backend, logger):
        self.backend = backend
        self.logger = logger

    def atom_labels(self, material: Material, std_atoms: ase.Atoms) -> None:
        material.atom_labels = std_atoms.get_chemical_symbols()

    def atom_positions(self, material: Material, std_atoms: ase.Atoms) -> None:
        material.atom_positions = std_atoms.get_scaled_positions(wrap=False)

    @abstractmethod
    def cell_normalized(self, material: Material, std_atoms: ase.Atoms) -> None:
        pass

    @abstractmethod
    def cell_primitive(self, material: Material, std_atoms: ase.Atoms) -> None:
        pass

    def cell_volume(self, calculation: Calculation, std_atoms: ase.Atoms) -> None:
        calculation.cell_volume = float(std_atoms.get_volume() * 1e-10**3)

    def formula(self, material: Material, names: List[str], counts: List[int]) -> None:
        formula = structure.get_formula_string(names, counts)
        material.formula = formula

    def formula_reduced(self, material: Material, names: list, counts_reduced: list) -> None:
        formula = structure.get_formula_string(names, counts_reduced)
        material.formula_reduced = formula

    def material_hash(self, material: Material, symmetry_analyzer: SymmetryAnalyzer) -> None:
        wyckoff_sets = symmetry_analyzer.get_wyckoff_sets_conventional()
        space_group_number = symmetry_analyzer.get_space_group_number()

        # Create and store hash based on SHA512
        norm_hash_string = structure.get_symmetry_string(space_group_number, wyckoff_sets)
        material.material_hash = sha512(norm_hash_string.encode('utf-8')).hexdigest()

    def number_of_atoms(self, material: Material, std_atoms: ase.Atoms) -> None:
        material.number_of_atoms = len(std_atoms)


class StructureBulk(Structure):
    """Processes structure related metainfo for Encyclopedia bulk structures.
    """
    def atomic_density(self, calculation: Calculation, repr_system: ase.Atoms) -> None:
        orig_n_atoms = len(repr_system)
        orig_volume = repr_system.get_volume() * (1e-10)**3
        calculation.atomic_density = float(orig_n_atoms / orig_volume)

    def bravais_lattice(self, material: Material, section_symmetry: Dict) -> None:
        bravais_lattice = section_symmetry["bravais_lattice"]
        material.bravais_lattice = bravais_lattice

    def cell_normalized(self, material: Material, std_atoms: ase.Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        material.cell_normalized = cell_normalized

    def cell_primitive(self, material: Material, prim_atoms: ase.Atoms) -> None:
        cell_prim = prim_atoms.get_cell()
        cell_prim *= 1e-10
        material.cell_primitive = cell_prim

    def crystal_system(self, material: Material, section_symmetry: Dict) -> None:
        material.crystal_system = section_symmetry["crystal_system"]

    def has_free_wyckoff_parameters(self, material: Material, symmetry_analyzer: SymmetryAnalyzer) -> None:
        has_free_param = symmetry_analyzer.get_has_free_wyckoff_parameters()
        material.has_free_wyckoff_parameters = has_free_param

    def lattice_parameters(self, calculation: Calculation, std_atoms: ase.Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        calculation.lattice_parameters = structure.get_lattice_parameters(cell_normalized)

    def mass_density(self, calculation: Calculation, repr_system: ase.Atoms) -> None:
        mass = structure.get_summed_atomic_mass(repr_system.get_atomic_numbers())
        orig_volume = repr_system.get_volume() * (1e-10)**3
        calculation.mass_density = float(mass / orig_volume)

    def material_name(self, material: Material, symbols: list, numbers: list) -> None:
        # Systems with one element are named after it
        if len(symbols) == 1:
            number = ase.data.atomic_numbers[symbols[0]]
            name = ase.data.atomic_names[number]
            material.material_name = name

        # Binary systems have specific names
        if len(symbols) == 2:
            atomicnumbers = [ase.data.atomic_numbers[i] for i in symbols]
            names = [ase.data.atomic_names[i] for i in atomicnumbers]

            # Non-metal elements are anions in the binary compounds and receive the -ide suffix
            if names[1] == "Antimony":
                names[1] = names[1][:-1] + "ide"
            if names[1] == "Arsenic":
                names[1] = names[1][:-1] + "de"
            if names[1] == "Boron" or names[1] == "Carbon":
                names[1] = names[1][:-2] + "ide"
            if names[1] == "Chlorine" or names[1] == "Germanium" or names[1] == "Selenium" or names[1] == "Bromine" \
               or names[1] == "Tellurium" or names[1] == "Iodine" or names[1] == "Polonium" or names[1] == "Astatine" or \
               names[1] == "Fluorine":
                names[1] = names[1][:-2] + "de"
            if names[1] == "Silicon" or names[1] == "Sulfur":
                names[1] = names[1][:-2] + "ide"
            if names[1] == "Nitrogen" or names[1] == "Oxygen" or names[1] == "Hydrogen" or names[1] == "Phosphorus":
                names[1] = names[1][:-4] + "ide"

            name = names[0] + " " + names[1]

            if names[1] == "Fluoride" or names[1] == "Chloride" or names[1] == "Bromide" or \
               names[1] == "Iodide" or names[1] == "Hydride":

                # Non-metals with elements of variable valence, therefore we remove alkaline and
                # alkaline-earth elements, which have fixed valence
                # Only the most electronegative non-metals are supposed to make ionic compounds
                if names[0] != "Lithium" and names[0] != "Sodium" and names[0] != "Potassium" and \
                   names[0] != "Rubidium" and names[0] != "Cesium" and names[0] != "Francium" and \
                   names[0] != "Beryllium" and names[0] != "Magnesium" and names[0] != "Calcium" and \
                   names[0] != "Strontium" and names[0] != "Barium" and names[0] != "Radium" and \
                   names[0] != "Aluminum":

                    if numbers[1] == 2:
                        name = names[0] + "(II)" + " " + names[1]
                    elif numbers[1] == 3:
                        name = names[0] + "(III)" + " " + names[1]
                    elif numbers[1] == 4:
                        name = names[0] + "(IV)" + " " + names[1]
                    elif numbers[1] == 5:
                        name = names[0] + "(V)" + " " + names[1]
                    elif numbers[1] == 6:
                        name = names[0] + "(VI)" + " " + names[1]
                    elif numbers[1] == 7:
                        name = names[0] + "(VII)" + " " + names[1]

            if names[1] == "Oxide" or names[1] == "Sulfide" or names[1] == "Selenide":
                if names[0] != "Lithium" and names[0] != "Sodium" and names[0] != "Potassium" and \
                   names[0] != "Rubidium" and names[0] != "Cesium" and names[0] != "Francium" and \
                   names[0] != "Beryllium" and names[0] != "Magnesium" and names[0] != "Calcium" and \
                   names[0] != "Strontium" and names[0] != "Barium" and names[0] != "Radium" and \
                   names[0] != "Aluminum":

                    if numbers[0] == 1 and numbers[1] == 1:
                        name = names[0] + "(II)" + " " + names[1]
                    elif numbers[0] == 2 and numbers[1] == 1:
                        name = names[0] + "(I)" + " " + names[1]
                    elif numbers[0] == 1 and numbers[1] == 2:
                        name = names[0] + "(IV)" + " " + names[1]
                    elif numbers[0] == 2 and numbers[1] == 3:
                        name = names[0] + "(III)" + " " + names[1]
                    elif numbers[0] == 2 and numbers[1] == 5:
                        name = names[0] + "(V)" + " " + names[1]
                    elif numbers[0] == 1 and numbers[1] == 3:
                        name = names[0] + "(VI)" + " " + names[1]
                    elif numbers[0] == 2 and numbers[1] == 7:
                        name = names[0] + "(VII)" + " " + names[1]

            if names[1] == "Nitride" or names[1] == "Phosphide":
                if names[0] != "Lithium" and names[0] != "Sodium" and names[0] != "Potassium" and \
                   names[0] != "Rubidium" and names[0] != "Cesium" and names[0] != "Francium" and \
                   names[0] != "Beryllium" and names[0] != "Magnesium" and names[0] != "Calcium" and \
                   names[0] != "Strontium" and names[0] != "Barium" and names[0] != "Radium" and \
                   names[0] != "Aluminum":

                    if numbers[0] == 1 and numbers[1] == 1:
                        name = names[0] + "(III)" + " " + names[1]
                    if numbers[0] == 1 and numbers[1] == 2:
                        name = names[0] + "(VI)" + " " + names[1]
                    elif numbers[0] == 3 and numbers[1] == 2:
                        name = names[0] + "(II)" + " " + names[1]
                    elif numbers[0] == 3 and numbers[1] == 4:
                        name = names[0] + "(IV)" + " " + names[1]
                    elif numbers[0] == 3 and numbers[1] == 5:
                        name = names[0] + "(V)" + " " + names[1]
                    elif numbers[0] == 3 and numbers[1] == 7:
                        name = names[0] + "(VII)" + " " + names[1]

            if names[1] == "Carbide":
                if names[0] != "Lithium" and names[0] != "Sodium" and names[0] != "Potassium" and \
                   names[0] != "Rubidium" and names[0] != "Cesium" and names[0] != "Francium" and \
                   names[0] != "Beryllium" and names[0] != "Magnesium" and names[0] != "Calcium" and \
                   names[0] != "Strontium" and names[0] != "Barium" and names[0] != "Radium" and \
                   names[0] != "Aluminum":

                    if numbers[0] == 1 and numbers[1] == 1:
                        name = names[0] + "(IV)" + " " + names[1]
                    if numbers[0] == 2 and numbers[1] == 1:
                        name = names[0] + "(II)" + " " + names[1]
                    if numbers[0] == 4 and numbers[1] == 1:
                        name = names[0] + "(I)" + " " + names[1]
                    if numbers[0] == 4 and numbers[1] == 3:
                        name = names[0] + "(III)" + " " + names[1]
                    if numbers[0] == 4 and numbers[1] == 5:
                        name = names[0] + "(V)" + " " + names[1]
                    if numbers[0] == 2 and numbers[1] == 3:
                        name = names[0] + "(VI)" + " " + names[1]
                    if numbers[0] == 4 and numbers[1] == 7:
                        name = names[0] + "(VII)" + " " + names[1]

            material.material_name = name

    def periodicity(self, material: Material) -> None:
        material.periodicity = np.array([0, 1, 2], dtype=np.int8)

    def point_group(self, material: Material, section_symmetry: Dict) -> None:
        point_group = section_symmetry["point_group"]
        material.point_group = point_group

    def space_group_number(self, material: Material, symmetry_analyzer: SymmetryAnalyzer) -> None:
        spg_number = symmetry_analyzer.get_space_group_number()
        material.space_group_number = spg_number

    def space_group_international_short_symbol(self, material: Material, symmetry_analyzer: SymmetryAnalyzer) -> None:
        spg_int_symb = symmetry_analyzer.get_space_group_international_short()
        material.space_group_international_short_symbol = spg_int_symb

    def material_classification(self, material: Material, section_system) -> None:
        try:
            sec_springer = section_system["section_springer_material"][0]
        except Exception:
            return

        classes: Dict[str, List[str]] = {}
        try:
            classifications = sec_springer['springer_classification']
        except KeyError:
            pass
        else:
            classes["material_class_springer"] = classifications
        try:
            compound_classes = sec_springer['springer_compound_class']
        except KeyError:
            pass
        else:
            classes["compound_class_springer"] = compound_classes
        if classes:
            material.material_classification = json.dumps(classes)

    def structure_type(self, material: Material, section_system) -> None:
        try:
            sec_prototype = section_system["section_prototype"][0]
            notes = sec_prototype.tmp['prototype_notes']
        except Exception:
            return

        # Only relevant information hidden in "notes" is handed over TODO:
        # review and eventually add more ****ites which are commonly used
        # (see wurzite)
        if notes in {
           'perovskite',
           '4-member ring',
           'fct',
           'bct',
           'bct5',
           'wurtzite',
           'hcp',
           'half-Heusler',
           'zincblende',
           'cubic perovskite',
           'simple cubic',
           'clathrate',
           'cuprite',
           'Heusler',
           'rock salt',
           'fcc',
           'diamond',
           'bcc'}:
            material.structure_type = notes

    def structure_prototype(self, material: Material, section_system) -> None:
        try:
            sec_prototype = section_system["section_prototype"][0]
            name = sec_prototype.tmp['prototype_name']
        except Exception:
            return

        material.structure_prototype = name

    def strukturbericht_designation(self, material: Material, section_system) -> None:
        try:
            sec_prototype = section_system["section_prototype"][0]
            strukturbericht = sec_prototype.tmp["strukturbericht_designation"]
        except Exception:
            return

        material.strukturbericht_designation = strukturbericht

    def wyckoff_groups(self, material: Material, wyckoff_sets: Dict) -> None:
        wyckoff_list = []
        for group in wyckoff_sets:
            data = {
                "wyckoff_letter": group.wyckoff_letter,
                "element": group.element,
                "indices": group.indices,
                "variables": {
                    "x": group.x,
                    "y": group.y,
                    "z": group.z,
                },
            }
            wyckoff_list.append(data)
        material.wyckoff_groups = json.dumps(wyckoff_list, sort_keys=True)

    def fill(self, sec_system) -> None:
        # Fetch resources
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        material = sec_enc.material
        calculation = sec_enc.calculation
        sec_symmetry = sec_system["section_symmetry"][0]
        symmetry_analyzer = sec_system["section_symmetry"][0].tmp["symmetry_analyzer"]
        std_atoms = symmetry_analyzer.get_conventional_system()
        prim_atoms = symmetry_analyzer.get_primitive_system()
        repr_atoms = sec_system.tmp["representative_atoms"]  # Temporary value stored by SystemNormalizer
        wyckoff_sets = symmetry_analyzer.get_wyckoff_sets_conventional()
        names, counts = structure.get_hill_decomposition(prim_atoms.get_chemical_symbols(), reduced=False)
        greatest_common_divisor = reduce(gcd, counts)
        reduced_counts = np.array(counts) / greatest_common_divisor

        # Fill structural information
        self.mass_density(calculation, repr_atoms)
        self.material_hash(material, symmetry_analyzer)
        self.number_of_atoms(material, std_atoms)
        self.atom_labels(material, std_atoms)
        self.atom_positions(material, std_atoms)
        self.atomic_density(calculation, repr_atoms)
        self.bravais_lattice(material, sec_symmetry)
        self.cell_normalized(material, std_atoms)
        self.cell_volume(calculation, std_atoms)
        self.crystal_system(material, sec_symmetry)
        self.cell_primitive(material, prim_atoms)
        self.formula(material, names, counts)
        self.formula_reduced(material, names, reduced_counts)
        self.has_free_wyckoff_parameters(material, symmetry_analyzer)
        self.lattice_parameters(calculation, std_atoms)
        self.material_name(material, names, reduced_counts)
        self.material_classification(material, sec_system)
        self.periodicity(material)
        self.point_group(material, sec_symmetry)
        self.space_group_number(material, symmetry_analyzer)
        self.space_group_international_short_symbol(material, symmetry_analyzer)
        self.structure_type(material, sec_system)
        self.structure_prototype(material, sec_system)
        self.strukturbericht_designation(material, sec_system)
        self.wyckoff_groups(material, wyckoff_sets)


class Structure2D(Structure):
    """Processes structure related metainfo for Encyclopedia 2D structures.
    """
    def cell_normalized(self, material: Material, std_atoms: ase.Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        material.cell_normalized = cell_normalized

    def cell_primitive(self, material: Material, prim_atoms: ase.Atoms) -> None:
        cell_prim = prim_atoms.get_cell()
        cell_prim *= 1e-10
        material.cell_primitive = cell_prim

    def lattice_parameters(self, calculation: Calculation, std_atoms: Atoms, non_periodic_index_std: int) -> None:
        # Eliminate the parameters in the non-periodic dimension
        periodicity = np.array([True, True, True])
        periodicity[non_periodic_index_std] = False
        full_parameters = structure.get_lattice_parameters(std_atoms.get_cell())
        lengths = np.array(full_parameters[:3])
        angles = np.array(full_parameters[3:])

        a, b = lengths[periodicity] * 1e-10
        alpha = angles[non_periodic_index_std]
        calculation.lattice_parameters = np.array([a, b, 0.0, alpha, 0.0, 0.0])

    def periodicity(self, material: Material, non_periodic_index_std: int) -> None:
        periodic_indices = [0, 1, 2]
        del periodic_indices[non_periodic_index_std]
        material.periodicity = np.array(periodic_indices, dtype=np.int8)

    def get_symmetry_analyzer(self, original_system: Atoms) -> SymmetryAnalyzer:
        # Determine the periodicity by examining vacuum gaps
        vacuum_directions = structure.find_vacuum_directions(original_system, threshold=7.0)
        periodicity = np.invert(vacuum_directions)

        # If two axis are not periodic, return. This only happens if the vacuum
        # gap is not aligned with a cell vector.
        if sum(periodicity) != 2:
            self.logger.warn("Could not detect the periodic dimensions in a 2D system.")
            return False

        # Center the system in the non-periodic direction, also taking
        # periodicity into account. The get_center_of_mass()-function in MatID
        # takes into account periodicity and can produce the correct CM unlike
        # the similar function in ASE.
        pbc_cm = matid.geometry.get_center_of_mass(original_system)
        cell_center = 0.5 * np.sum(original_system.get_cell(), axis=0)
        translation = cell_center - pbc_cm
        translation[periodicity] = 0
        symm_system = original_system.copy()
        symm_system.translate(translation)
        symm_system.wrap()

        # Set the periodicity according to detected periodicity in order for
        # SymmetryAnalyzer to use the symmetry analysis designed for 2D
        # systems.
        symm_system.set_pbc(periodicity)
        symmetry_analyzer = SymmetryAnalyzer(
            symm_system,
            config.normalize.symmetry_tolerance,
            config.normalize.flat_dim_threshold
        )
        return symmetry_analyzer

    def get_non_periodic_index_std(self, symmetry_analyzer: SymmetryAnalyzer) -> int:
        # Get the periodicity as detected by classification
        pbc = symmetry_analyzer._original_system.get_pbc()
        non_periodic_index_orig = np.where(pbc == False)[0][0]  # noqa: E712

        # The index of the originally non-periodic dimension may not correspond
        # to the one in the normalized system, because the normalized system
        # may use a different coordinate system. We will have to get this
        # transformation_matrix = self.symmetry_dataset["transformation_matrix"]
        transformation_matrix = symmetry_analyzer._get_spglib_transformation_matrix()
        for i_axis, axis in enumerate(transformation_matrix):
            if axis[non_periodic_index_orig] != 0 and \
               axis[(non_periodic_index_orig + 1) % 3] == 0.0 and \
               axis[(non_periodic_index_orig + 2) % 3] == 0.0:
                non_periodic_index_std = i_axis
                break

        return non_periodic_index_std


    def fill(self, representative_system) -> None:
        # Fetch resources
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        material = sec_enc.material
        calculation = sec_enc.calculation
        repr_atoms = representative_system.tmp["representative_atoms"]  # Temporary value stored by SystemNormalizer
        symmetry_analyzer = self.get_symmetry_analyzer(repr_atoms)
        std_atoms = symmetry_analyzer.get_conventional_system()
        prim_atoms = symmetry_analyzer.get_primitive_system()
        names, counts = structure.get_hill_decomposition(prim_atoms.get_chemical_symbols(), reduced=False)
        greatest_common_divisor = reduce(gcd, counts)
        reduced_counts = np.array(counts) / greatest_common_divisor
        non_periodic_index_std = self.get_non_periodic_index_std(symmetry_analyzer)

        # Fill metainfo
        self.material_hash(material, symmetry_analyzer)
        self.number_of_atoms(material, std_atoms)
        self.atom_labels(material, std_atoms)
        self.atom_positions(material, std_atoms)
        self.cell_normalized(material, std_atoms)
        self.cell_primitive(material, prim_atoms)
        self.formula(material, names, counts)
        self.formula_reduced(material, names, reduced_counts)
        self.lattice_parameters(calculation, std_atoms, non_periodic_index_std)
        self.periodicity(material, non_periodic_index_std)


class Structure1D(Structure):
    """Processes structure related metainfo for Encyclopedia 1D structures.
    """
    def material_hash(self, material: Material, prim_atoms: Atoms) -> None:
        """Hash to be used as identifier for a material. Different 1D
        materials are defined by their Coulomb matrix eigenvalues and their
        Hill formulas.
        """
        fingerprint = self.get_structure_fingerprint(prim_atoms)
        formula = material.formula
        id_strings = []
        id_strings.append(formula)
        id_strings.append(fingerprint)
        hash_seed = ", ".join(id_strings)
        hash_val = sha512(hash_seed.encode('utf-8')).hexdigest()
        material.material_hash = hash_val

    def cell_normalized(self, material: Material, std_atoms: ase.Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        material.cell_normalized = cell_normalized

    def lattice_parameters(self, calculation: Calculation, std_atoms: Atoms, periodic_indices: np.array) -> None:
        # 1D systems only have one lattice parameter: length in periodic dimension
        cell = std_atoms.get_cell()
        a = np.linalg.norm(cell[periodic_indices[0], :]) * 1e-10
        calculation.lattice_parameters = np.array([a, 0.0, 0.0, 0.0, 0.0, 0.0])

    def periodicity(self, material: Material, prim_atoms: Atoms) -> None:
        # Determine the periodicity by examining vacuum gaps
        vacuum_directions = structure.find_vacuum_directions(
            prim_atoms,
            threshold=config.normalize.cluster_threshold
        )

        # If one axis is not periodic, return. This only happens if the vacuum
        # gap is not aligned with a cell vector.
        if sum(vacuum_directions) != 2:
            raise ValueError("Could not detect the periodic dimensions in a 1D system.")

        periodic_indices = np.where(vacuum_directions == False)[0]  # noqa: E712
        material.periodicity = np.array(periodic_indices, dtype=np.int8)

    def get_structure_fingerprint(self, prim_atoms: Atoms):
        """Calculates a numeric fingerprint that coarsely encodes the atomic
        positions and species.
        """
        # Calculate charge part
        q = prim_atoms.get_atomic_numbers()
        qiqj = np.sqrt(q[None, :] * q[:, None])

        # Calculate distance part. Notice that the minimum image convention
        # must be used. Without it, differently oriented atoms in the same cell
        # may be detected as the same material.
        pos = prim_atoms.get_positions()
        cell = prim_atoms.get_cell()
        cmat = 10 - matid.geometry.get_distance_matrix(pos, pos, cell, pbc=True, mic=True)
        cmat = np.clip(cmat, a_min=0, a_max=None)
        np.fill_diagonal(cmat, 0)
        cmat = qiqj * cmat

        # Calculate eigenvalues
        eigval, _ = np.linalg.eigh(cmat)

        # Sort eigenvalues
        eigval = np.array(sorted(eigval))

        # Go to smaller scale
        eigval /= 100

        # Create a slightly shifted versions of the values. The shift can go up or
        # down. If both of these shifted values gets rounded to the same number,
        # then the rounding is considered robust. Otherwise the rounding will be
        # done to a value that is between two consequent rounding points.
        padding = 0.1
        up_round = np.rint(eigval + padding)
        down_round = np.rint(eigval - padding)
        mask = (up_round == down_round)

        new_vals = np.array(up_round)
        for i, value in enumerate(mask):
            if not value:
                new_val = (up_round[i] + down_round[i]) / 2
                new_vals[i] = new_val

        # The values should be rounded to a rational number for the hashing to be
        # consistent
        strings = []
        for number in new_vals:
            num_str = str(number)
            strings.append(num_str)
        fingerprint = "; ".join(strings)

        return fingerprint

    def get_symmetry_analyzer(self, original_system: Atoms) -> SymmetryAnalyzer:
        """For 1D systems the symmetery is analyzed from the original system
        with enforced full periodicity.

        Args:
            original_system: The original simulation system.

        Returns:
            The SymmetryAnalyzer that is instantiated with the original system.
        """
        symm_system = original_system.copy()
        symm_system.set_pbc(True)
        symmetry_analyzer = SymmetryAnalyzer(
            symm_system,
            config.normalize.symmetry_tolerance,
            config.normalize.flat_dim_threshold
        )

        return symmetry_analyzer

    def get_std_atoms(self, periodicity: np.array, prim_atoms: Atoms) -> Atoms:
        """For 1D systems the standardized system is based on a primitive
        system. This primitive system is translated to the center of mass and
        the non-periodic dimensions are minimized so that the atoms just fit.

        Args:
            periodicity: List of periodic indices, in 1D case a list containing
                one index.
            prim_atoms: Primitive system

        Returns
            Standardized structure that represents this material and from which
            the material hash will be constructed from.
        """
        std_atoms = prim_atoms.copy()

        # Translate to center of mass
        pbc_cm = matid.geometry.get_center_of_mass(prim_atoms)
        cell_center = 0.5 * np.sum(std_atoms.get_cell(), axis=0)
        translation = cell_center - pbc_cm
        translation[periodicity] = 0
        std_atoms.translate(translation)
        std_atoms.wrap()

        # Reduce cell size to just fit the system in the non-periodic dimensions.
        indices = [0, 1, 2]
        for idx in periodicity:
            del indices[idx]
        pos = std_atoms.get_scaled_positions(wrap=False)
        cell = std_atoms.get_cell()
        new_cell = np.array(cell)
        translation = np.zeros(3)
        for index in indices:
            imin = np.min(pos[:, index])
            imax = np.max(pos[:, index])
            translation -= cell[index, :] * imin
            new_cell[index] = cell[index, :] * (imax - imin)
        std_atoms.translate(translation)
        std_atoms.set_cell(new_cell)

        return std_atoms

    def fill(self, representative_system) -> None:
        # Fetch resources
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        material = sec_enc.material
        calculation = sec_enc.calculation
        repr_atoms = representative_system.tmp["representative_atoms"]  # Temporary value stored by SystemNormalizer
        symmetry_analyzer = self.get_symmetry_analyzer(repr_atoms)
        prim_atoms = symmetry_analyzer.get_primitive_system()
        prim_atoms.set_pbc(True)
        names, counts = structure.get_hill_decomposition(prim_atoms.get_chemical_symbols(), reduced=False)
        greatest_common_divisor = reduce(gcd, counts)
        reduced_counts = np.array(counts) / greatest_common_divisor

        # Fill metainfo
        self.periodicity(material, prim_atoms)
        std_atoms = self.get_std_atoms(material.periodicity, prim_atoms)
        self.number_of_atoms(material, std_atoms)
        self.atom_labels(material, std_atoms)
        self.atom_positions(material, std_atoms)
        self.cell_normalized(material, std_atoms)
        self.formula(material, names, counts)
        self.formula_reduced(material, names, reduced_counts)
        self.material_hash(material, std_atoms)
        self.lattice_parameters(calculation, std_atoms, material.periodicity)
