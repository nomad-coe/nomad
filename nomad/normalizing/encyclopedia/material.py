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

from typing import Dict, List
from nptyping import NDArray
from math import gcd, isnan
from functools import reduce
from abc import abstractmethod
import re
import json
import ase
import ase.data
from ase import Atoms
import numpy as np
from matid import SymmetryAnalyzer
import matid.geometry

from nomad.datamodel.encyclopedia import (
    Material,
    Properties,
    Bulk,
    IdealizedStructure,
    WyckoffSet,
    WyckoffVariables,
    LatticeParameters,
)
from nomad.normalizing.encyclopedia.context import Context
from nomad.metainfo import Section
from nomad import atomutils
from nomad.utils import hash
from nomad import config

J_to_Ry = 4.587425e+17


class MaterialNormalizer():
    """A base class that is used for processing material-related information
    in the Encylopedia.
    """
    def __init__(self, entry_archive, logger):
        self.logger = logger
        self.entry_archive = entry_archive

    def atom_labels(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        ideal.atom_labels = std_atoms.get_chemical_symbols()

    def atom_positions(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        ideal.atom_positions = std_atoms.get_scaled_positions(wrap=False)

    @abstractmethod
    def lattice_vectors(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        pass

    def cell_volume(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        ideal.cell_volume = float(std_atoms.get_volume() * 1e-10**3)

    def formula(self, material: Material, names: List[str], counts: List[int]) -> None:
        formula = atomutils.get_formula_string(names, counts)
        material.formula = formula

    def formula_reduced(self, material: Material, names: list, counts_reduced: list) -> None:
        formula = atomutils.get_formula_string(names, counts_reduced)
        material.formula_reduced = formula

    def species_and_counts(self, material: Material, names: List[str], reduced_counts: List[int]) -> None:
        parts = []
        for name, count in zip(names, reduced_counts):
            if count == 1:
                parts.append(name)
            else:
                parts.append("{}{}".format(name, int(count)))
        material.species_and_counts = " ".join(parts)

    def species(self, material: Material, names: List[str]) -> None:
        material.species = " ".join(names)

    def material_id(self, material: Material, spg_number: int, wyckoff_sets: List[WyckoffSet]) -> None:
        # Create and store hash based on SHA512
        norm_hash_string = atomutils.get_symmetry_string(spg_number, wyckoff_sets)
        material.material_id = hash(norm_hash_string)

    def number_of_atoms(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        ideal.number_of_atoms = len(std_atoms)

    @abstractmethod
    def normalize(self, context: Context) -> None:
        pass


class MaterialBulkNormalizer(MaterialNormalizer):
    """Processes structure related metainfo for Encyclopedia bulk structures.
    """
    def atomic_density(self, properties: Properties, repr_system: Atoms) -> None:
        orig_n_atoms = len(repr_system)
        orig_volume = repr_system.get_volume() * (1e-10)**3
        properties.atomic_density = float(orig_n_atoms / orig_volume)

    def bravais_lattice(self, bulk: Bulk, section_symmetry: Section) -> None:
        bravais_lattice = section_symmetry["bravais_lattice"]
        bulk.bravais_lattice = bravais_lattice

    def lattice_vectors(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        ideal.lattice_vectors = cell_normalized

    def lattice_vectors_primitive(self, ideal: IdealizedStructure, prim_atoms: Atoms) -> None:
        cell_prim = prim_atoms.get_cell()
        cell_prim *= 1e-10
        ideal.lattice_vectors_primitive = cell_prim

    def crystal_system(self, bulk: Bulk, section_symmetry: Section) -> None:
        bulk.crystal_system = section_symmetry["crystal_system"]

    def has_free_wyckoff_parameters(self, bulk: Bulk, symmetry_analyzer: SymmetryAnalyzer) -> None:
        has_free_param = symmetry_analyzer.get_has_free_wyckoff_parameters()
        bulk.has_free_wyckoff_parameters = has_free_param

    def lattice_parameters(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        cell_normalized = std_atoms.get_cell() * 1E-10
        param_values = atomutils.cell_to_cellpar(cell_normalized)
        param_section = ideal.m_create(LatticeParameters)
        param_section.a = float(param_values[0])
        param_section.b = float(param_values[1])
        param_section.c = float(param_values[2])
        param_section.alpha = float(param_values[3])
        param_section.beta = float(param_values[4])
        param_section.gamma = float(param_values[5])

    def mass_density(self, properties: Properties, repr_system: Atoms) -> None:
        mass = atomutils.get_summed_atomic_mass(repr_system.get_atomic_numbers())
        orig_volume = repr_system.get_volume() * (1e-10)**3
        mass_density = float(mass / orig_volume)
        if isnan(mass_density):
            properties.mass_density = 0
        else:
            properties.mass_density = mass_density

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

    def periodicity(self, ideal: IdealizedStructure) -> None:
        ideal.periodicity = np.array([True, True, True], dtype=np.bool_)

    def point_group(self, bulk: Bulk, section_symmetry: Section) -> None:
        point_group = section_symmetry["point_group"]
        bulk.point_group = point_group

    def space_group_number(self, bulk: Bulk, spg_number: int) -> None:
        bulk.space_group_number = spg_number

    def space_group_international_short_symbol(self, bulk: Bulk, symmetry_analyzer: SymmetryAnalyzer) -> None:
        spg_int_symb = symmetry_analyzer.get_space_group_international_short()
        bulk.space_group_international_short_symbol = spg_int_symb

    def material_classification(self, material: Material, section_system: Section) -> None:
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

    def structure_type(self, bulk: Bulk, section_system: Section) -> None:
        try:
            sec_prototype = section_system["section_prototype"][0]
            notes = sec_prototype.m_cache['prototype_notes']
        except Exception:
            return

        # Only relevant information hidden in "notes" is handed over TODO:
        # review and eventually add more ****ites which are commonly used
        # (see wurzite)
        note_map = {
            "CaTiO<sub>3</sub> Pnma Perovskite Structure": "perovskite",
            "Hypothetical Tetrahedrally Bonded Carbon with 4&ndash;Member Rings": "4-member ring",
            "In (A6) Structure": "fct",
            "$\\alpha$&ndash;Pa (A<sub>a</sub>) Structure": "bct",
            "Hypothetical BCT5 Si Structure": "bct5",
            "Wurtzite (ZnS, B4) Structure": "wurtzite",
            "Hexagonal Close Packed (Mg, A3) Structure": "hcp",
            "Half&ndash;Heusler (C1<sub>b</sub>) Structure": "half-Heusler",
            "Zincblende (ZnS, B3) Structure": "zincblende",
            "Cubic Perovskite (CaTiO<sub>3</sub>, E2<sub>1</sub>) Structure": "cubic perovskite",
            "$\\alpha$&ndash;Po (A<sub>h</sub>) Structure": "simple cubic",
            "Si<sub>46</sub> Clathrate Structure": "clathrate",
            "Cuprite (Cu<sub>2</sub>O, C3) Structure": "cuprite",
            "Heusler (L2<sub>1</sub>) Structure": "Heusler",
            "Rock Salt (NaCl, B1) Structure": "rock salt",
            "Face&ndash;Centered Cubic (Cu, A1) Structure": "fcc",
            "Diamond (A4) Structure": "diamond",
            "Body&ndash;Centered Cubic (W, A2) Structure": "bcc",
        }
        enc_note = note_map.get(notes, None)
        if enc_note is not None:
            bulk.structure_type = enc_note

    def structure_prototype(self, bulk: Bulk, section_system: Section) -> None:
        try:
            sec_prototype = section_system["section_prototype"][0]
            name = sec_prototype.m_cache['prototype_name']
        except Exception:
            return

        bulk.structure_prototype = name

    def strukturbericht_designation(self, bulk: Bulk, section_system: Section) -> None:
        try:
            sec_prototype = section_system["section_prototype"][0]
            strukturbericht = sec_prototype.m_cache["strukturbericht_designation"]
        except Exception:
            return

        # In the current GUI we replace LaTeX with plain text
        strukturbericht = re.sub('[$_{}]', '', strukturbericht)
        bulk.strukturbericht_designation = strukturbericht

    def wyckoff_sets(self, ideal: IdealizedStructure, wyckoff_sets: Dict) -> None:
        for group in wyckoff_sets:
            wset = ideal.m_create(WyckoffSet)
            if group.x is not None or group.y is not None or group.z is not None:
                variables = wset.m_create(WyckoffVariables)
                if group.x is not None:
                    variables.x = float(group.x)
                if group.y is not None:
                    variables.y = float(group.y)
                if group.z is not None:
                    variables.z = float(group.z)
            wset.indices = group.indices
            wset.element = group.element
            wset.wyckoff_letter = group.wyckoff_letter

    def normalize(self, context: Context) -> None:
        # Fetch resources
        sec_system = context.representative_system
        sec_enc = self.entry_archive.section_metadata.encyclopedia
        material = sec_enc.material
        properties = sec_enc.properties
        sec_symmetry = sec_system["section_symmetry"][0]
        symmetry_analyzer = sec_system["section_symmetry"][0].m_cache["symmetry_analyzer"]
        spg_number = symmetry_analyzer.get_space_group_number()
        std_atoms = symmetry_analyzer.get_conventional_system()
        prim_atoms = symmetry_analyzer.get_primitive_system()
        repr_atoms = sec_system.m_cache["representative_atoms"]  # Temporary value stored by SystemNormalizer
        try:
            wyckoff_sets = symmetry_analyzer.get_wyckoff_sets_conventional(return_parameters=True)
        except Exception:
            self.logger.error('Error resolving Wyckoff sets.')
            wyckoff_sets = []
        names, counts = atomutils.get_hill_decomposition(prim_atoms.get_chemical_symbols(), reduced=False)
        greatest_common_divisor = reduce(gcd, counts)
        context.greatest_common_divisor = greatest_common_divisor
        reduced_counts = np.array(counts) / greatest_common_divisor

        # Fill structural information
        bulk = material.m_create(Bulk)
        ideal = material.m_create(IdealizedStructure)
        self.mass_density(properties, repr_atoms)
        self.material_id(material, spg_number, wyckoff_sets)
        self.number_of_atoms(ideal, std_atoms)
        self.atom_labels(ideal, std_atoms)
        self.atom_positions(ideal, std_atoms)
        self.atomic_density(properties, repr_atoms)
        self.bravais_lattice(bulk, sec_symmetry)
        self.lattice_vectors(ideal, std_atoms)
        self.cell_volume(ideal, std_atoms)
        self.crystal_system(bulk, sec_symmetry)
        self.lattice_vectors_primitive(ideal, prim_atoms)
        self.formula(material, names, counts)
        self.formula_reduced(material, names, reduced_counts)  # type: ignore
        self.species(material, names)
        self.species_and_counts(material, names, reduced_counts)  # type: ignore
        self.has_free_wyckoff_parameters(bulk, symmetry_analyzer)
        self.lattice_parameters(ideal, std_atoms)
        self.material_name(material, names, reduced_counts)  # type: ignore
        self.material_classification(material, sec_system)
        self.periodicity(ideal)
        self.point_group(bulk, sec_symmetry)
        self.space_group_number(bulk, spg_number)
        self.space_group_international_short_symbol(bulk, symmetry_analyzer)
        self.structure_type(bulk, sec_system)
        self.structure_prototype(bulk, sec_system)
        self.strukturbericht_designation(bulk, sec_system)
        self.wyckoff_sets(ideal, wyckoff_sets)


class Material2DNormalizer(MaterialNormalizer):
    """Processes structure related metainfo for Encyclopedia 2D structures.
    """
    def material_id(self, material: Material, spg_number: int, wyckoff_sets: List[WyckoffSet]) -> None:
        # The hash is based on the symmetry analysis of the structure when it
        # is treated as a 3D structure. Due to this the hash may overlap with
        # real 3D structures unless we include a distinguishing label for 2D
        # structures in the hash seed.
        norm_hash_string = atomutils.get_symmetry_string(spg_number, wyckoff_sets, is_2d=True)
        material.material_id = hash(norm_hash_string)

    def lattice_vectors(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        ideal.lattice_vectors = cell_normalized

    def lattice_vectors_primitive(self, ideal: IdealizedStructure, prim_atoms: Atoms) -> None:
        cell_prim = prim_atoms.get_cell()
        cell_prim *= 1e-10
        ideal.lattice_vectors_primitive = cell_prim

    def lattice_parameters(self, ideal: IdealizedStructure, std_atoms: Atoms, periodicity: NDArray) -> None:
        # 2D systems only have three lattice parameter: two lengths and angle between them
        periodic_indices = np.where(np.array(periodicity) == True)[0]  # noqa: E712
        cell = std_atoms.get_cell()
        a_vec = cell[periodic_indices[0], :] * 1e-10
        b_vec = cell[periodic_indices[1], :] * 1e-10
        a = np.linalg.norm(a_vec)
        b = np.linalg.norm(b_vec)
        gamma = np.clip(np.dot(a_vec, b_vec) / (a * b), -1.0, 1.0)
        gamma = np.arccos(gamma)
        param_section = ideal.m_create(LatticeParameters)
        param_section.a = float(a)
        param_section.b = float(b)
        param_section.gamma = float(gamma)

    def periodicity(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        # MatID already provides the correct periodicity
        ideal.periodicity = std_atoms.get_pbc()

    def get_symmetry_analyzer(self, original_system: Atoms) -> SymmetryAnalyzer:
        # Get dimension of system by also taking into account the covalent radii
        dimensions = matid.geometry.get_dimensions(original_system, [True, True, True])
        basis_dimensions = np.linalg.norm(original_system.get_cell(), axis=1)
        gaps = basis_dimensions - dimensions
        periodicity = gaps <= config.normalize.cluster_threshold

        # If two axis are not periodic, return. This only happens if the vacuum
        # gap is not aligned with a cell vector or if the linear gap search is
        # unsufficient (the structure is "wavy" making also the gap highly
        # nonlinear).
        if sum(periodicity) != 2:
            raise ValueError("Could not detect the periodic dimensions in a 2D system.")

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

    def normalize(self, context: Context) -> None:
        # Fetch resources
        sec_enc = self.entry_archive.section_metadata.encyclopedia
        material = sec_enc.material
        repr_atoms = context.representative_system.m_cache["representative_atoms"]  # Temporary value stored by SystemNormalizer
        try:
            symmetry_analyzer = self.get_symmetry_analyzer(repr_atoms)
        except Exception:
            self.logger.error('Error setting up symmetry analyzer.')
            return

        spg_number = symmetry_analyzer.get_space_group_number()
        wyckoff_sets = symmetry_analyzer.get_wyckoff_sets_conventional(return_parameters=False)
        std_atoms = symmetry_analyzer.get_conventional_system()
        prim_atoms = symmetry_analyzer.get_primitive_system()
        names, counts = atomutils.get_hill_decomposition(prim_atoms.get_chemical_symbols(), reduced=False)
        greatest_common_divisor = reduce(gcd, counts)
        context.greatest_common_divisor = greatest_common_divisor
        reduced_counts = np.array(counts) / greatest_common_divisor

        # Fill metainfo
        ideal = material.m_create(IdealizedStructure)
        self.periodicity(ideal, std_atoms)
        self.material_id(material, spg_number, wyckoff_sets)
        self.number_of_atoms(ideal, std_atoms)
        self.atom_labels(ideal, std_atoms)
        self.atom_positions(ideal, std_atoms)
        self.lattice_vectors(ideal, std_atoms)
        self.lattice_vectors_primitive(ideal, prim_atoms)
        self.formula(material, names, counts)
        self.formula_reduced(material, names, reduced_counts)  # type: ignore
        self.species(material, names)
        self.species_and_counts(material, names, reduced_counts)  # type: ignore
        self.lattice_parameters(ideal, std_atoms, ideal.periodicity)


class Material1DNormalizer(MaterialNormalizer):
    """Processes structure related metainfo for Encyclopedia 1D structures.
    """
    def material_id_1d(self, material: Material, prim_atoms: Atoms) -> None:
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
        material.material_id = hash_val

    def lattice_vectors(self, ideal: IdealizedStructure, std_atoms: Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        ideal.lattice_vectors = cell_normalized

    def lattice_parameters(self, ideal: IdealizedStructure, std_atoms: Atoms, periodicity: NDArray) -> None:
        # 1D systems only have one lattice parameter: length in periodic dimension
        periodic_indices = np.where(np.array(periodicity) == True)[0]  # noqa: E712
        if len(periodic_indices) == 0:
            return
        cell = std_atoms.get_cell()
        a = np.linalg.norm(cell[periodic_indices[0], :]) * 1e-10
        params = ideal.m_create(LatticeParameters)
        params.a = float(a)

    def periodicity(self, ideal: IdealizedStructure, prim_atoms: Atoms) -> None:
        # Get dimension of system by also taking into account the covalent radii
        dimensions = matid.geometry.get_dimensions(prim_atoms, [True, True, True])
        basis_dimensions = np.linalg.norm(prim_atoms.get_cell(), axis=1)
        gaps = basis_dimensions - dimensions
        periodicity = gaps <= config.normalize.cluster_threshold

        # If one axis is not periodic, return. This only happens if the vacuum
        # gap is not aligned with a cell vector.
        if sum(periodicity) != 1:
            self.logger.error("Could not detect the periodic dimensions in a 1D system.")

        ideal.periodicity = periodicity

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
        """For 1D systems the symmetry is analyzed from the original system
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

    def get_std_atoms(self, periodicity: NDArray, prim_atoms: Atoms) -> Atoms:
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
        pos = std_atoms.get_scaled_positions(wrap=False)
        cell = std_atoms.get_cell()
        new_cell = np.array(cell)
        translation = np.zeros(3)
        for index, periodic in enumerate(periodicity):
            if not periodic:
                imin = np.min(pos[:, index])
                imax = np.max(pos[:, index])
                translation -= cell[index, :] * imin
                new_cell[index] = cell[index, :] * (imax - imin)
        std_atoms.translate(translation)
        std_atoms.set_cell(new_cell)

        return std_atoms

    def normalize(self, context: Context) -> None:
        # Fetch resources
        sec_system = context.representative_system
        sec_enc = self.entry_archive.section_metadata.encyclopedia
        material = sec_enc.material
        repr_atoms = sec_system.m_cache["representative_atoms"]  # Temporary value stored by SystemNormalizer
        symmetry_analyzer = self.get_symmetry_analyzer(repr_atoms)
        prim_atoms = symmetry_analyzer.get_primitive_system()
        prim_atoms.set_pbc(True)
        names, counts = atomutils.get_hill_decomposition(prim_atoms.get_chemical_symbols(), reduced=False)
        greatest_common_divisor = reduce(gcd, counts)
        context.greatest_common_divisor = greatest_common_divisor
        reduced_counts = np.array(counts) / greatest_common_divisor

        # Fill metainfo
        ideal = material.m_create(IdealizedStructure)
        self.periodicity(ideal, prim_atoms)
        std_atoms = self.get_std_atoms(ideal.periodicity, prim_atoms)
        self.number_of_atoms(ideal, std_atoms)
        self.atom_labels(ideal, std_atoms)
        self.atom_positions(ideal, std_atoms)
        self.lattice_vectors(ideal, std_atoms)
        self.formula(material, names, counts)
        self.formula_reduced(material, names, reduced_counts)  # type: ignore
        self.species(material, names)
        self.species_and_counts(material, names, reduced_counts)  # type: ignore
        self.material_id_1d(material, std_atoms)
        self.lattice_parameters(ideal, std_atoms, ideal.periodicity)
