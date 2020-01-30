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

from typing import Dict, List
from math import gcd as gcd
from functools import reduce
from abc import abstractmethod
from collections import OrderedDict
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
from nomad.utils import hash
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

    def method_type(self) -> str:
        pass

    def mainfile_uri(self, calculation: Calculation):
        entry_info = self._backend["section_entry_info"][0]
        upload_id = entry_info["upload_id"]
        mainfile_path = entry_info["mainfile"]
        uri = f"nmd://R{upload_id}/data/{mainfile_path}"
        calculation.mainfile_uri = uri

    # def similar_materials(self) -> None:
        # pass

    # def calculation_pid(self):
        # pass

    # def calculation(self) -> None:
        # pass

    # def contributor_first_name(self) -> None:
        # pass

    # def contributor_last_name(self) -> None:
        # pass

    # def contributor_type(self) -> None:
        # pass

    # def contributors(self) -> None:
        # pass

    # def number_of_calculations(self) -> None:
        # pass

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

        # Fill method related metainfo
        method = Method(self._backend, self.logger)
        method.fill()

    def normalize(self, logger=None) -> None:
        super().normalize(logger)
        system_enums = Material.system_type.type

        # Initialise metainfo structure
        sec_enc = Encyclopedia()
        material = sec_enc.m_create(Material)
        calculation = sec_enc.m_create(Calculation)

        # Get generic data
        self.mainfile_uri(calculation)

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

        # Process all present properties
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
        material.material_hash = hash(norm_hash_string)

    def number_of_atoms(self, material: Material, std_atoms: ase.Atoms) -> None:
        material.number_of_atoms = len(std_atoms)

    @abstractmethod
    def fill(self, sec_system) -> None:
        pass


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

    def lattice_parameters(self, calculation: Calculation, std_atoms: Atoms, periodic_indices: np.array) -> None:
        # 2D systems only have three lattice parameter: two length and angle between them
        cell = std_atoms.get_cell()
        a_vec = cell[periodic_indices[0], :] * 1e-10
        b_vec = cell[periodic_indices[1], :] * 1e-10
        a = np.linalg.norm(a_vec)
        b = np.linalg.norm(b_vec)
        alpha = np.clip(np.dot(a_vec, b_vec) / (a * b), -1.0, 1.0)
        alpha = np.arccos(alpha)
        calculation.lattice_parameters = np.array([a, b, 0.0, alpha, 0.0, 0.0])

    def periodicity(self, material: Material, std_atoms: Atoms) -> None:
        # Determine the periodicity by examining vacuum gaps
        vacuum_directions = structure.find_vacuum_directions(
            std_atoms,
            threshold=config.normalize.cluster_threshold
        )

        # If one axis is not periodic, return. This only happens if the vacuum
        # gap is not aligned with a cell vector.
        if sum(vacuum_directions) != 1:
            raise ValueError("Could not detect the periodic dimensions in a 2D system.")

        periodic_indices = np.where(vacuum_directions == False)[0]  # noqa: E712
        material.periodicity = np.sort(np.array(periodic_indices, dtype=np.int8))

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
        # non_periodic_index_std = self.get_non_periodic_index_std(symmetry_analyzer)

        # Fill metainfo
        self.periodicity(material, std_atoms)
        self.material_hash(material, symmetry_analyzer)
        self.number_of_atoms(material, std_atoms)
        self.atom_labels(material, std_atoms)
        self.atom_positions(material, std_atoms)
        self.cell_normalized(material, std_atoms)
        self.cell_primitive(material, prim_atoms)
        self.formula(material, names, counts)
        self.formula_reduced(material, names, reduced_counts)
        self.lattice_parameters(calculation, std_atoms, material.periodicity)


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
        hash_val = hash(hash_seed)
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
        material.periodicity = np.sort(np.array(periodic_indices, dtype=np.int8))

    def get_structure_fingerprint(self, prim_atoms: Atoms) -> str:
        """Calculates a numeric fingerprint that coarsely encodes the atomic
        positions and species.

        The fingerprint is based on calculating a discretized version of a
        sorted Coulomb matrix eigenspectrum (Grégoire Montavon, Katja Hansen,
        Siamac Fazli, Matthias Rupp, Franziska Biegler, Andreas Ziehe,
        Alexandre Tkatchenko, Anatole V. Lilienfeld, and Klaus-Robert Müller.
        Learning invariant representations of molecules for atomization energy
        prediction. In F. Pereira, C. J. C. Burges, L. Bottou, and K. Q.
        Weinberger, editors, Advances in Neural Information Processing Systems
        25, pages 440–448. Curran Associates, Inc., 2012.).

        The fingerprints are discretized in order to perform O(n) matching
        between structures (no need to compare fingerprints against each
        other). As regular discretization is susceptible to the "edge problem",
        a robust discretization is used instead (Birget, Jean-Camille & Hong,
        Dawei & Memon, Nasir. (2003). Robust discretization, with an
        application to graphical passwords. IACR Cryptology ePrint Archive.
        2003. 168.) Basically for the 1-dimensional domain two grids are
        created and the points are mapped to the first grid in which they are
        robust using a minimum tolerance parameter r, with the maximum
        tolerance being 5r.

        There are other robust discretization methods that can guarantee exact
        r-tolerance (e.g. Sonia Chiasson, Jayakumar Srinivasan, Robert Biddle,
        and P. C. van Oorschot. 2008. Centered discretization with application
        to graphical passwords. In Proceedings of the 1st Conference on
        Usability, Psychology, and Security (UPSEC’08). USENIX Association,
        USA, Article 6, 1–9.). This method however requires that a predefined
        "correct" structure exists against which the search is done.

        Args:
            prim_atoms: Primitive system.

        Returns:
            The numeric fingerprint for the system encoded as a string.
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

        # Perform robust discretization (see function docstring for details). r
        # = 0.5 ensures that all grids are integers which can be uniquely
        # mapped to strings. If finer grid is needed adjust the eigenvalue scale
        # instead.
        eigval /= 25  # Go to smaller scale where integer numbers are meaningful
        dimension = 1
        r = 0.5
        spacing = 2 * r * (dimension + 1)
        phi_k = 2 * r * np.array(range(dimension + 1))
        t = np.mod((eigval[None, :] + phi_k[:, None]), (2 * r * (dimension + 1)))
        grid_mask = (r <= t) & (t < r * (2 * dimension + 1))
        safe_grid_k = np.argmax(grid_mask == True, axis=0)   # noqa: E712
        discretization = spacing * np.floor((eigval + (2 * r * safe_grid_k)) / spacing)
        discretization[safe_grid_k == 1] += 2 * r

        # Form string
        strings = []
        for number in discretization:
            num_str = str(int(number))
            strings.append(num_str)
        fingerprint = ";".join(strings)

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


class Method():
    """A base class that is used for processing method related information
    in the Encylopedia.
    """
    def __init__(self, backend, logger):
        self.backend = backend
        self.logger = logger

    def code_name(self, calculation: Calculation) -> None:
        calculation.code_name = self.backend["program_name"]

    def code_version(self, calculation: Calculation) -> None:
        calculation.code_version = self.backend["program_version"]

    def method_hash(self, calculation: Calculation):
        # Create ordered dictionary with the values. Order is important for
        # consistent hashing. Not all the settings are set for this
        # calculation, in which case they will be simply set as None.
        hash_dict: OrderedDict = OrderedDict()
        hash_dict['program_name'] = calculation.code_name
        hash_dict['program_version'] = calculation.code_version

        # The subclasses may define their own method properties that are to be
        # included here.
        hash_dict.update(self.method_hash_dict())

        # Form a hash from the dictionary
        method_hash = self.group_dict_to_hash("method_hash", hash_dict)
        calculation.method_hash = method_hash

    def group_eos_hash(self, calculation: Calculation):
        # Create ordered dictionary with the values. Order is important for
        # consistent hashing.
        hash_dict: OrderedDict = OrderedDict()
        hash_dict['upload_id'] = self.backend["section_entry_info"][0]["upload_id"]  # Only calculations from the same upload are grouped

        # The subclasses may define their own method properties that are to be
        # included here.
        hash_dict.update(self.group_eos_hash_dict())

        # Form a hash from the dictionary
        group_eos_hash = self.group_dict_to_hash('group_eos_hash', hash_dict)
        calculation.group_eos_hash = group_eos_hash

    def group_parametervariation_hash(self, calculation: Calculation):
        # Create ordered dictionary with the values. Order is important for
        # consistent hashing.
        hash_dict: OrderedDict = OrderedDict()
        hash_dict['upload_id'] = self.backend["section_entry_info"][0]["upload_id"]  # Only calculations from the same upload are grouped

        # The subclasses may define their own method properties that are to be
        # included here.
        hash_dict.update(self.group_parametervariation_hash_dict())

        # Form a hash from the dictionary
        group_eos_hash = self.group_dict_to_hash('group_eos_hash', hash_dict)
        calculation.group_parametervariation_hash = group_eos_hash

    def group_e_min(self) -> None:
        pass

    def group_type(self) -> None:
        pass

    def k_point_grid_description(self) -> None:
        pass

    def basis_set_short_name(self) -> None:
        pass

    def basis_set_type(self) -> None:
        pass

    def core_electron_treatment(self) -> None:
        pass

    def functional_long_name(self) -> None:
        pass

    def functional_type(self) -> None:
        pass

    def gw_starting_point(self) -> None:
        pass

    def gw_type(self) -> None:
        pass

    def pseudopotential_type(self) -> None:
        pass

    def scf_threshold(self) -> None:
        pass

    def smearing_kind(self) -> None:
        pass

    def smearing_parameters(self) -> None:
        pass

    @abstractmethod
    def method_hash_dict(self):
        return OrderedDict()

    @abstractmethod
    def group_eos_hash_dict(self):
        return OrderedDict()

    @abstractmethod
    def group_parametervariation_hash_dict(self):
        return OrderedDict()

    def group_dict_to_hash(self, name, src_dict: OrderedDict):
        """Given a dictionary of computational settings, this function forms a
        hash code from them.
        """
        nones = self.find_nones(src_dict)
        if len(nones) > 0:
            self.logger.warning(
                '%s: missing data for hash: %s',
                name, ', '.join(sorted(nones))
            )
        hash_str = json.dumps(src_dict)

        return hash(hash_str)

    def find_nones(self, data, parent='') -> List[str]:
        """Recursively finds values that are set as None. Returns a list of
        identifiers for the None-values.
        """
        result: List[str] = []
        if data is None:
            return [parent]
        elif isinstance(data, (bytes, str)):
            return []
        if getattr(data, 'items', None) is not None:
            for (k, v) in data.items():
                result += self.find_nones(v, "%s['%s']" % (parent, k))
        elif getattr(data, '__len__', None) is not None:
            for i in range(len(data)):
                result += self.find_nones(data[i], "%s[%d]" % (parent, i))
        return result

    def fill(self) -> None:
        # Fetch resources
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        calculation = sec_enc.calculation

        # Fill metainfo
        self.code_name(calculation)
        self.code_version(calculation)
        self.method_hash(calculation)


class Properties():
    """A base class that is used for processing information that is specific to
    a type of calculation.
    """
    def __init__(self, backend, logger):
        self.backend = backend
        self.logger = logger

    def band_gap(self) -> None:
        pass

    def band_gap_position(self) -> None:
        pass

    def band_gap_type(self) -> None:
        pass

    def band_structure(self) -> None:
        pass

    def brillouin_zone(self) -> None:
        pass

    def brillouin_zone_viewer(self) -> None:
        pass

    def dos(self) -> None:
        pass

    def elastic_constants_matrix(self) -> None:
        pass

    def elastic_deformation_energies(self) -> None:
        pass

    def elastic_fitting_parameters(self) -> None:
        pass

    def elastic_moduli(self) -> None:
        pass

    def elastic_properties(self) -> None:
        pass

    def fermi_surface(self) -> None:
        pass

    def has_bs(self) -> None:
        pass

    def has_dos(self) -> None:
        pass

    def has_fermi_surface(self) -> None:
        pass

    def has_thermal_properties(self) -> None:
        pass

    def phonon_dispersion(self) -> None:
        pass

    def phonon_dos(self) -> None:
        pass

    def specific_heat_cv(self) -> None:
        pass

    def helmholtz_free_energy(self) -> None:
        pass

    def energies(self) -> None:
        pass

    @abstractmethod
    def fill(self, sec_system) -> None:
        pass
