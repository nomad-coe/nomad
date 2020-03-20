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

from typing import Dict, List, Any
from math import gcd as gcd
from functools import reduce
from abc import abstractmethod
from collections import OrderedDict
import re
import json
import ase
import ase.data
from ase import Atoms
import numpy as np
from matid import SymmetryAnalyzer
import matid.geometry

from nomad.normalizing.normalizer import (
    Normalizer,
    s_run,
    s_scc,
    s_system,
    s_method,
    s_frame_sequence,
    r_frame_sequence_to_sampling,
    s_sampling_method,
    r_frame_sequence_local_frames,
)
from nomad.metainfo.encyclopedia import (
    Encyclopedia,
    Material,
    Method,
    Properties,
    RunType,
    WyckoffSet,
    WyckoffVariables,
    ElectronicBandStructure,
    BandGap,
)
from nomad.parsing.backend import Section, LocalBackend
from nomad.normalizing.settingsbasisset import get_basis_set_settings
from nomad.normalizing import structure
from nomad.utils import hash, RestrictedDict
from nomad import config

J_to_Ry = 4.587425e+17


class Context():
    """A simple class for holding the context related to an Encylopedia entry.
    """
    def __init__(
        self,
        material_type: str,
        method_type: str,
        run_type: str,
        representative_system,
        representative_method,
        representative_scc,
        representative_scc_idx,
    ):
        self.material_type = material_type
        self.method_type = method_type
        self.run_type = run_type
        self.representative_system = representative_system
        self.representative_method = representative_method
        self.representative_scc = representative_scc
        self.representative_scc_idx = representative_scc_idx
        self.greatest_common_divisor: int = None


class EncyclopediaNormalizer(Normalizer):
    """
    This normalizer emulates the functionality of the old Encyclopedia backend.
    The data used by the encyclopedia have been assigned under new metainfo
    within a new section called "Encyclopedia". In the future these separate
    metainfos could be absorbed into the existing metainfo hiearchy.
    """
    def __init__(self, backend: LocalBackend):
        super().__init__(backend)
        self.backend: LocalBackend = backend

    def run_type(self, run_type_sec: RunType) -> str:
        """Decides what type of calculation this is: single_point, md,
        geometry_optimization, etc.
        """
        run_enums = RunType.run_type.type
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

        # No sequences, only a few calculations
        if n_scc <= 3 and n_frame_seq == 0:
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
                    "value for frame_sequence_to_sampling_ref."
                )
                return run_type

            # See if local frames are present
            try:
                frames = frame_seq[r_frame_sequence_local_frames]
            except KeyError:
                self.logger.info(
                    "section_frame_sequence_local_frames not found although a "
                    "frame_sequence exists."
                )
                return run_type
            if len(frames) == 0:
                self.logger.info("No frames referenced in section_frame_sequence_local_frames.")
                return run_type

            section_sampling_method = self._backend[s_sampling_method][i_sampling_method]
            sampling_method = section_sampling_method["sampling_method"]

            if sampling_method == "molecular_dynamics":
                run_type = run_enums.molecular_dynamics
            if sampling_method == "geometry_optimization":
                run_type = run_enums.geometry_optimization
            if sampling_method == "taylor_expansion":
                run_type = run_enums.phonon_calculation

        run_type_sec.run_type = run_type
        return run_type

    def material_type(self, material: Material) -> tuple:
        # Try to fetch representative system
        system = None
        material_type = config.services.unavailable_value
        material_enums = Material.material_type.type
        system_idx = self._backend["section_run"][0].tmp["representative_system_idx"]
        if system_idx is not None:
            # Try to find system type information from backend for the selected system.
            try:
                system = self._backend[s_system][system_idx]
                stype = system["system_type"]
            except KeyError:
                pass
            else:
                if stype == material_enums.one_d or stype == material_enums.two_d:
                    material_type = stype
                # For bulk systems we also ensure that the symmetry information is available
                if stype == material_enums.bulk:
                    try:
                        system["section_symmetry"][0]
                    except (KeyError, IndexError):
                        self.logger.info("Symmetry information is not available for a bulk system. No Encylopedia entry created.")
                    else:
                        material_type = stype

        material.material_type = material_type
        return system, material_type

    def method_type(self, method: Method) -> tuple:
        repr_method = None
        method_id = config.services.unavailable_value
        methods = self._backend[s_method]
        n_methods = len(methods)

        if n_methods == 1:
            repr_method = methods[0]
            method_id = repr_method.get("electronic_structure_method", config.services.unavailable_value)
        elif n_methods > 1:
            for sec_method in self._backend[s_method]:
                # GW
                electronic_structure_method = sec_method.get("electronic_structure_method", None)
                if electronic_structure_method in {"G0W0", "scGW"}:
                    repr_method = sec_method
                    method_id = "GW"
                    break

                # Methods linked to each other through references. Get all
                # linked methods, try to get electronic_structure_method from
                # each.
                try:
                    refs = sec_method["section_method_to_method_refs"]
                except KeyError:
                    pass
                else:
                    linked_methods = [sec_method]
                    for ref in refs:
                        method_to_method_kind = ref["method_to_method_kind"]
                        method_to_method_ref = ref["method_to_method_ref"]
                        if method_to_method_kind == "core_settings":
                            linked_methods.append(methods[method_to_method_ref])

                    for i_method in linked_methods:
                        try:
                            electronic_structure_method = i_method["electronic_structure_method"]
                        except KeyError:
                            pass
                        else:
                            repr_method = sec_method
                            method_id = electronic_structure_method

        method.method_type = method_id
        return repr_method, method_id

    def mainfile_uri(self, encyclopedia: Encyclopedia):
        entry_info = self._backend["section_entry_info"][0]
        upload_id = entry_info["upload_id"]
        mainfile_path = entry_info["mainfile"]
        uri = f"nmd://R{upload_id}/data/{mainfile_path}"
        encyclopedia.mainfile_uri = uri

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

    def fill(self, ctx: Context):
        # Fill structure related metainfo
        struct: Any = None
        if ctx.material_type == Material.material_type.type.bulk:
            struct = MaterialBulkNormalizer(self.backend, self.logger)
        elif ctx.material_type == Material.material_type.type.two_d:
            struct = Material2DNormalizer(self.backend, self.logger)
        elif ctx.material_type == Material.material_type.type.one_d:
            struct = Material1DNormalizer(self.backend, self.logger)
        if struct is not None:
            struct.normalize(ctx)

        # Fill method related metainfo
        method = None
        if ctx.method_type == Method.method_type.type.DFT or ctx.method_type == Method.method_type.type.DFTU:
            method = MethodDFTNormalizer(self._backend, self.logger)
        elif ctx.method_type == Method.method_type.type.GW:
            method = MethodGWNormalizer(self._backend, self.logger)
        if method is not None:
            method.normalize(ctx)

        # Fill properties related metainfo
        properties = PropertiesNormalizer(self.backend, self.logger)
        properties.normalize(ctx)

    def normalize(self, logger=None) -> None:
        """The caller will automatically log if the normalizer succeeds or ends
        up with an exception.
        """
        try:
            super().normalize(logger)

            # Initialise metainfo structure
            sec_enc = Encyclopedia()
            material = sec_enc.m_create(Material)
            method = sec_enc.m_create(Method)
            sec_enc.m_create(Properties)
            run_type = sec_enc.m_create(RunType)

            # Get generic data
            self.mainfile_uri(sec_enc)

            # Determine run type, stop if unknown
            run_type_name = self.run_type(run_type)
            if run_type_name == config.services.unavailable_value:
                self.logger.info(
                    "Unsupported run type for encyclopedia, encyclopedia metainfo not created.",
                    enc_status="unsupported_run_type",
                )
                return

            # Get the system type, stop if unknown
            material_enums = Material.material_type.type
            representative_system, material_type = self.material_type(material)
            if material_type != material_enums.bulk and material_type != material_enums.two_d and material_type != material_enums.one_d:
                self.logger.info(
                    "Unsupported system type for encyclopedia, encyclopedia metainfo not created.",
                    enc_status="unsupported_material_type",
                )
                return

            # Get the method type, stop if unknown
            representative_method, method_type = self.method_type(method)
            if method_type == config.services.unavailable_value:
                self.logger.info(
                    "Unsupported method type for encyclopedia, encyclopedia metainfo not created.",
                    enc_status="unsupported_method_type",
                )
                return

            # Get representative scc
            try:
                representative_scc_idx = self._backend[s_run][0].tmp["representative_scc_idx"]
                representative_scc = self._backend[s_scc][representative_scc_idx]
            except (KeyError, IndexError):
                representative_scc = None
                representative_scc_idx = None

            # Create one context that holds all details
            context = Context(
                material_type=material_type,
                method_type=method_type,
                run_type=run_type_name,
                representative_system=representative_system,
                representative_method=representative_method,
                representative_scc=representative_scc,
                representative_scc_idx=representative_scc_idx,
            )

            # Put the encyclopedia section into backend
            self._backend.add_mi2_section(sec_enc)
            self.fill(context)
        except Exception:
            self.logger.error(
                "Failed to create an Encyclopedia entry due to an unhandlable exception.",
                enc_status="failure",
            )
            raise  # Reraise for the caller to log the exception as well
        else:
            self.logger.info(
                "Successfully created metainfo for Encyclopedia.",
                enc_status="success",
            )


class MaterialNormalizer():
    """A base class that is used for processing material-related information
    in the Encylopedia.
    """
    def __init__(self, backend: LocalBackend, logger):
        self.backend = backend
        self.logger = logger

    def atom_labels(self, material: Material, std_atoms: Atoms) -> None:
        material.atom_labels = std_atoms.get_chemical_symbols()

    def atom_positions(self, material: Material, std_atoms: Atoms) -> None:
        material.atom_positions = std_atoms.get_scaled_positions(wrap=False)

    @abstractmethod
    def cell_normalized(self, material: Material, std_atoms: Atoms) -> None:
        pass

    def cell_volume(self, material: Material, std_atoms: Atoms) -> None:
        material.cell_volume = float(std_atoms.get_volume() * 1e-10**3)

    def formula(self, material: Material, names: List[str], counts: List[int]) -> None:
        formula = structure.get_formula_string(names, counts)
        material.formula = formula

    def formula_reduced(self, material: Material, names: list, counts_reduced: list) -> None:
        formula = structure.get_formula_string(names, counts_reduced)
        material.formula_reduced = formula

    def material_hash(self, material: Material, spg_number: int, wyckoff_sets: List[WyckoffSet]) -> None:
        # Create and store hash based on SHA512
        norm_hash_string = structure.get_symmetry_string(spg_number, wyckoff_sets)
        material.material_hash = hash(norm_hash_string)

    def number_of_atoms(self, material: Material, std_atoms: Atoms) -> None:
        material.number_of_atoms = len(std_atoms)

    @abstractmethod
    def normalize(self, ctx: Context) -> None:
        pass


class MaterialBulkNormalizer(MaterialNormalizer):
    """Processes structure related metainfo for Encyclopedia bulk structures.
    """
    def atomic_density(self, properties: Properties, repr_system: Atoms) -> None:
        orig_n_atoms = len(repr_system)
        orig_volume = repr_system.get_volume() * (1e-10)**3
        properties.atomic_density = float(orig_n_atoms / orig_volume)

    def bravais_lattice(self, material: Material, section_symmetry: Section) -> None:
        bravais_lattice = section_symmetry["bravais_lattice"]
        material.bravais_lattice = bravais_lattice

    def cell_normalized(self, material: Material, std_atoms: Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        material.cell_normalized = cell_normalized

    def cell_primitive(self, material: Material, prim_atoms: Atoms) -> None:
        cell_prim = prim_atoms.get_cell()
        cell_prim *= 1e-10
        material.cell_primitive = cell_prim

    def crystal_system(self, material: Material, section_symmetry: Section) -> None:
        material.crystal_system = section_symmetry["crystal_system"]

    def has_free_wyckoff_parameters(self, material: Material, symmetry_analyzer: SymmetryAnalyzer) -> None:
        has_free_param = symmetry_analyzer.get_has_free_wyckoff_parameters()
        material.has_free_wyckoff_parameters = has_free_param

    def lattice_parameters(self, material: Material, std_atoms: Atoms) -> None:
        cell_normalized = std_atoms.get_cell() * 1E-10
        material.lattice_parameters = structure.get_lattice_parameters(cell_normalized)

    def mass_density(self, properties: Properties, repr_system: Atoms) -> None:
        mass = structure.get_summed_atomic_mass(repr_system.get_atomic_numbers())
        orig_volume = repr_system.get_volume() * (1e-10)**3
        properties.mass_density = float(mass / orig_volume)

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
        material.periodicity = np.array([True, True, True], dtype=np.bool)

    def point_group(self, material: Material, section_symmetry: Section) -> None:
        point_group = section_symmetry["point_group"]
        material.point_group = point_group

    def space_group_number(self, material: Material, spg_number: int) -> None:
        material.space_group_number = spg_number

    def space_group_international_short_symbol(self, material: Material, symmetry_analyzer: SymmetryAnalyzer) -> None:
        spg_int_symb = symmetry_analyzer.get_space_group_international_short()
        material.space_group_international_short_symbol = spg_int_symb

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

    def structure_type(self, material: Material, section_system: Section) -> None:
        try:
            sec_prototype = section_system["section_prototype"][0]
            notes = sec_prototype.tmp['prototype_notes']
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
            material.structure_type = enc_note

    def structure_prototype(self, material: Material, section_system: Section) -> None:
        try:
            sec_prototype = section_system["section_prototype"][0]
            name = sec_prototype.tmp['prototype_name']
        except Exception:
            return

        material.structure_prototype = name

    def strukturbericht_designation(self, material: Material, section_system: Section) -> None:
        try:
            sec_prototype = section_system["section_prototype"][0]
            strukturbericht = sec_prototype.tmp["strukturbericht_designation"]
        except Exception:
            return

        # In the current GUI we replace LaTeX with plain text
        strukturbericht = re.sub('[$_{}]', '', strukturbericht)
        material.strukturbericht_designation = strukturbericht

    def wyckoff_sets(self, material: Material, wyckoff_sets: Dict) -> None:
        for group in wyckoff_sets:
            wset = material.m_create(WyckoffSet)
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

    def normalize(self, ctx: Context) -> None:
        # Fetch resources
        sec_system = ctx.representative_system
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        material = sec_enc.material
        properties = sec_enc.properties
        sec_symmetry = sec_system["section_symmetry"][0]
        symmetry_analyzer = sec_system["section_symmetry"][0].tmp["symmetry_analyzer"]
        spg_number = symmetry_analyzer.get_space_group_number()
        std_atoms = symmetry_analyzer.get_conventional_system()
        prim_atoms = symmetry_analyzer.get_primitive_system()
        repr_atoms = sec_system.tmp["representative_atoms"]  # Temporary value stored by SystemNormalizer
        wyckoff_sets = symmetry_analyzer.get_wyckoff_sets_conventional(return_parameters=True)
        names, counts = structure.get_hill_decomposition(prim_atoms.get_chemical_symbols(), reduced=False)
        greatest_common_divisor = reduce(gcd, counts)
        ctx.greatest_common_divisor = greatest_common_divisor
        reduced_counts = np.array(counts) / greatest_common_divisor

        # Fill structural information
        self.mass_density(properties, repr_atoms)
        self.material_hash(material, spg_number, wyckoff_sets)
        self.number_of_atoms(material, std_atoms)
        self.atom_labels(material, std_atoms)
        self.atom_positions(material, std_atoms)
        self.atomic_density(properties, repr_atoms)
        self.bravais_lattice(material, sec_symmetry)
        self.cell_normalized(material, std_atoms)
        self.cell_volume(material, std_atoms)
        self.crystal_system(material, sec_symmetry)
        self.cell_primitive(material, prim_atoms)
        self.formula(material, names, counts)
        self.formula_reduced(material, names, reduced_counts)
        self.has_free_wyckoff_parameters(material, symmetry_analyzer)
        self.lattice_parameters(material, std_atoms)
        self.material_name(material, names, reduced_counts)
        self.material_classification(material, sec_system)
        self.periodicity(material)
        self.point_group(material, sec_symmetry)
        self.space_group_number(material, spg_number)
        self.space_group_international_short_symbol(material, symmetry_analyzer)
        self.structure_type(material, sec_system)
        self.structure_prototype(material, sec_system)
        self.strukturbericht_designation(material, sec_system)
        self.wyckoff_sets(material, wyckoff_sets)


class Material2DNormalizer(MaterialNormalizer):
    """Processes structure related metainfo for Encyclopedia 2D structures.
    """
    def cell_normalized(self, material: Material, std_atoms: Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        material.cell_normalized = cell_normalized

    def cell_primitive(self, material: Material, prim_atoms: Atoms) -> None:
        cell_prim = prim_atoms.get_cell()
        cell_prim *= 1e-10
        material.cell_primitive = cell_prim

    def lattice_parameters(self, material: Material, std_atoms: Atoms, periodicity: np.array) -> None:
        # 2D systems only have three lattice parameter: two length and angle between them
        periodic_indices = np.where(np.array(periodicity) == True)[0]  # noqa: E712
        cell = std_atoms.get_cell()
        a_vec = cell[periodic_indices[0], :] * 1e-10
        b_vec = cell[periodic_indices[1], :] * 1e-10
        a = np.linalg.norm(a_vec)
        b = np.linalg.norm(b_vec)
        alpha = np.clip(np.dot(a_vec, b_vec) / (a * b), -1.0, 1.0)
        alpha = np.arccos(alpha)
        material.lattice_parameters = np.array([a, b, 0.0, alpha, 0.0, 0.0])

    def periodicity(self, material: Material, std_atoms: Atoms) -> None:
        # MatID already provides the correct periodicity
        material.periodicity = std_atoms.get_pbc()

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

    def normalize(self, ctx: Context) -> None:
        # Fetch resources
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        material = sec_enc.material
        repr_atoms = ctx.representative_system.tmp["representative_atoms"]  # Temporary value stored by SystemNormalizer
        symmetry_analyzer = self.get_symmetry_analyzer(repr_atoms)
        spg_number = symmetry_analyzer.get_space_group_number()
        wyckoff_sets = symmetry_analyzer.get_wyckoff_sets_conventional(return_parameters=False)
        std_atoms = symmetry_analyzer.get_conventional_system()
        prim_atoms = symmetry_analyzer.get_primitive_system()
        names, counts = structure.get_hill_decomposition(prim_atoms.get_chemical_symbols(), reduced=False)
        greatest_common_divisor = reduce(gcd, counts)
        ctx.greatest_common_divisor = greatest_common_divisor
        reduced_counts = np.array(counts) / greatest_common_divisor

        # Fill metainfo
        self.periodicity(material, std_atoms)
        self.material_hash(material, spg_number, wyckoff_sets)
        self.number_of_atoms(material, std_atoms)
        self.atom_labels(material, std_atoms)
        self.atom_positions(material, std_atoms)
        self.cell_normalized(material, std_atoms)
        self.cell_primitive(material, prim_atoms)
        self.formula(material, names, counts)
        self.formula_reduced(material, names, reduced_counts)
        self.lattice_parameters(material, std_atoms, material.periodicity)


class Material1DNormalizer(MaterialNormalizer):
    """Processes structure related metainfo for Encyclopedia 1D structures.
    """
    def material_hash_1d(self, material: Material, prim_atoms: Atoms) -> None:
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

    def cell_normalized(self, material: Material, std_atoms: Atoms) -> None:
        cell_normalized = std_atoms.get_cell()
        cell_normalized *= 1e-10
        material.cell_normalized = cell_normalized

    def lattice_parameters(self, material: Material, std_atoms: Atoms, periodicity: np.array) -> None:
        # 1D systems only have one lattice parameter: length in periodic dimension
        periodic_indices = np.where(np.array(periodicity) == True)[0]  # noqa: E712
        cell = std_atoms.get_cell()
        a = np.linalg.norm(cell[periodic_indices[0], :]) * 1e-10
        material.lattice_parameters = np.array([a, 0.0, 0.0, 0.0, 0.0, 0.0])

    def periodicity(self, material: Material, prim_atoms: Atoms) -> None:
        # Get dimension of system by also taking into account the covalent radii
        dimensions = matid.geometry.get_dimensions(prim_atoms, [True, True, True])
        basis_dimensions = np.linalg.norm(prim_atoms.get_cell(), axis=1)
        gaps = basis_dimensions - dimensions
        periodicity = gaps <= config.normalize.cluster_threshold

        # If one axis is not periodic, return. This only happens if the vacuum
        # gap is not aligned with a cell vector.
        if sum(periodicity) != 1:
            raise ValueError("Could not detect the periodic dimensions in a 1D system.")

        material.periodicity = periodicity

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

    def normalize(self, ctx: Context) -> None:
        # Fetch resources
        sec_system = ctx.representative_system
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        material = sec_enc.material
        repr_atoms = sec_system.tmp["representative_atoms"]  # Temporary value stored by SystemNormalizer
        symmetry_analyzer = self.get_symmetry_analyzer(repr_atoms)
        prim_atoms = symmetry_analyzer.get_primitive_system()
        prim_atoms.set_pbc(True)
        names, counts = structure.get_hill_decomposition(prim_atoms.get_chemical_symbols(), reduced=False)
        greatest_common_divisor = reduce(gcd, counts)
        ctx.greatest_common_divisor = greatest_common_divisor
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
        self.material_hash_1d(material, std_atoms)
        self.lattice_parameters(material, std_atoms, material.periodicity)


class MethodNormalizer():
    """A base class that is used for processing method related information
    in the Encylopedia.
    """
    def __init__(self, backend, logger):
        self.backend = backend
        self.logger = logger

    def code_name(self, method: Method) -> None:
        method.code_name = self.backend["program_name"]

    def code_version(self, method: Method) -> None:
        method.code_version = self.backend["program_version"]

    def method_hash(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section):
        method_dict = RestrictedDict(
            mandatory_keys=[
                "program_name",
                "subsettings",
            ],
            forbidden_values=[None]
        )
        method_dict['program_name'] = method.code_name

        # The subclasses may define their own method properties that are to be
        # included here.
        subsettings = self.method_hash_dict(method, settings_basis_set, repr_method)
        method_dict["subsettings"] = subsettings

        # If all required information is present, safe the hash
        try:
            method_dict.check(recursive=True)
        except (KeyError, ValueError) as e:
            self.logger.info("Could not create method hash, missing required information: {}".format(str(e)))
        else:
            method.method_hash = method_dict.hash()

    @abstractmethod
    def method_hash_dict(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section) -> RestrictedDict:
        pass

    def group_eos_hash(self, method: Method, material: Material, repr_method: Section):
        eos_dict = RestrictedDict(
            mandatory_keys=(
                "upload_id",
                "method_hash",
                "formula",
            ),
            forbidden_values=[None]
        )

        # Only calculations from the same upload are grouped
        eos_dict['upload_id'] = self.backend["section_entry_info"][0]["upload_id"]

        # Method
        eos_dict["method_hash"] = method.method_hash

        # The formula should be same for EoS (maybe even symmetries)
        eos_dict["formula"] = material.formula

        # Form a hash from the dictionary
        try:
            eos_dict.check(recursive=True)
        except (KeyError, ValueError) as e:
            self.logger.info("Could not create EOS hash, missing required information: {}".format(str(e)))
        else:
            method.group_eos_hash = eos_dict.hash()

    def group_parametervariation_hash(self, method: Method, settings_basis_set: RestrictedDict, repr_system: Section, repr_method: Section):
        # Create ordered dictionary with the values. Order is important for
        param_dict = RestrictedDict(
            mandatory_keys=[
                "upload_id",
                "program_name",
                "program_version",
                "settings_geometry",
                "subsettings",
            ],
            forbidden_values=[None]
        )

        # Only calculations from the same upload are grouped
        param_dict['upload_id'] = self.backend["section_entry_info"][0]["upload_id"]

        # The same code and functional type is required
        param_dict['program_name'] = method.code_name
        param_dict['program_version'] = method.code_version

        # Get a string representation of the geometry. It is included as the
        # geometry should remain the same during parameter variation. By simply
        # using the atom labels and positions we assume that their
        # order/translation/rotation does not change.
        geom_dict: OrderedDict = OrderedDict()
        sec_sys = repr_system
        atom_labels = sec_sys['atom_labels']
        geom_dict['atom_labels'] = ', '.join(atom_labels)
        atom_positions = sec_sys['atom_positions']
        geom_dict['atom_positions'] = np.array2string(
            atom_positions * 1e10,  # convert to Angstrom
            formatter={'float_kind': lambda x: "%.6f" % x},
        ).replace('\n', '')
        cell = sec_sys['lattice_vectors']
        geom_dict['simulation_cell'] = np.array2string(
            cell * 1e10,  # convert to Angstrom
            formatter={'float_kind': lambda x: "%.6f" % x},
        ).replace('\n', '')
        param_dict['settings_geometry'] = geom_dict

        # The subclasses may define their own method properties that are to be
        # included here.
        subsettings = self.group_parametervariation_hash_dict(method, settings_basis_set, repr_method)
        param_dict["subsettings"] = subsettings

        # Form a hash from the dictionary
        try:
            param_dict.check(recursive=True)
        except (KeyError, ValueError) as e:
            self.logger.info("Could not create parameter variation hash, missing required information: {}".format(str(e)))
        else:
            method.group_parametervariation_hash = param_dict.hash()

    @abstractmethod
    def group_parametervariation_hash_dict(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section) -> RestrictedDict:
        pass

    def group_e_min(self) -> None:
        pass

    def group_type(self) -> None:
        pass

    @abstractmethod
    def normalize(self, ctx: Context) -> None:
        pass


class MethodDFTNormalizer(MethodNormalizer):
    """A base class that is used for processing method related information
    in the Encylopedia.
    """
    def basis_set_type(self, method: Method) -> None:
        """Type of basis set used by the code"""
        basis_set_type = config.services.unavailable_value
        archive_basis_set = self.backend["program_basis_set_type"]
        # TODO: this translation should not be necessary if parsers did their
        # work correctly (using the metainfo doc as reference) until then we
        # keep this mapping.
        basis_set_type_ambiguity = {
            "numeric AOs": "Numeric AOs",
            "gaussians": "Gaussians",
            "plane waves": "Plane waves",
            "plane_waves": "Plane waves",
            "real-space grid": "Real-space grid"
        }
        if archive_basis_set in basis_set_type_ambiguity:
            self.logger.info(
                "Basis set type '{}' does not correspond to valid options in "
                "metainfo documentation and was corrected."
                .format(basis_set_type)
            )
            archive_basis_set = basis_set_type_ambiguity[archive_basis_set]

        if archive_basis_set is not None:
            basis_set_type = archive_basis_set

        method.basis_set_type = basis_set_type

    def core_electron_treatment(self, method: Method) -> None:
        treatment = config.services.unavailable_value
        code_name = method.code_name
        if code_name is not None:
            core_electron_treatments = {
                'VASP': 'pseudopotential',
                'FHI-aims': 'full all electron',
                'exciting': 'full all electron',
                'quantum espresso': 'pseudopotential'
            }
            treatment = core_electron_treatments.get(code_name, config.services.unavailable_value)
        method.core_electron_treatment = treatment

    def functional_long_name(self, method: Method, repr_method: Section) -> None:
        """'Long' form of exchange-correlation functional, list of components
        and parameters as a string: see
        https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional
        """
        xc_functional = MethodDFTNormalizer.functional_long_name_from_method(repr_method, self.backend[s_method])
        if xc_functional is config.services.unavailable_value:
            self.logger.warning(
                "Metainfo for 'XC_functional' not found, and could not "
                "compose name from 'section_XC_functionals'."
            )
        method.functional_long_name = xc_functional

    @staticmethod
    def functional_long_name_from_method(repr_method: Section, methods: List[Section]):
        """'Long' form of exchange-correlation functional, list of components
        and parameters as a string: see
        https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info/wikis/metainfo/XC-functional
        """
        linked_methods = [repr_method]
        try:
            refs = repr_method["section_method_to_method_refs"]
        except KeyError:
            pass
        else:
            for ref in refs:
                method_to_method_kind = ref["method_to_method_kind"]
                method_to_method_ref = ref["method_to_method_ref"]
                if method_to_method_kind == "core_settings":
                    linked_methods.append(methods[method_to_method_ref])

        xc_functional = config.services.unavailable_value
        for method in linked_methods:
            try:
                section_xc_functionals = method["section_XC_functionals"]
            except KeyError:
                pass
            else:
                components = {}
                for component in section_xc_functionals:
                    try:
                        cname = component['XC_functional_name']
                    except KeyError:
                        pass
                    else:
                        this_component = ''
                        if 'XC_functional_weight' in component:
                            this_component = str(component['XC_functional_weight']) + '*'
                        this_component += cname
                        components[cname] = this_component
                result_array = []
                for name in sorted(components):
                    result_array.append(components[name])
                if len(result_array) >= 1:
                    xc_functional = '+'.join(result_array)

        return xc_functional

    def functional_type(self, method: Method) -> None:
        long_name = method.functional_long_name
        if long_name is not None:
            short_name = self.create_xc_functional_shortname(long_name)
            method.functional_type = short_name

    def smearing_kind(self, method: Method, representative_method: Section) -> None:
        try:
            smearing_kind = representative_method['smearing_kind']
        except KeyError:
            pass
        else:
            method.smearing_kind = smearing_kind

    def smearing_parameter(self, method: Method, representative_method) -> None:
        try:
            smearing_width = representative_method['smearing_width']
        except KeyError:
            pass
        else:
            method.smearing_parameter = smearing_width

    def method_hash_dict(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section) -> RestrictedDict:
        # Extend by DFT settings.
        hash_dict = RestrictedDict(
            mandatory_keys=(
                "functional_long_name",
                "settings_basis_set",
                "scf_threshold_energy_change",
            ),
            optional_keys=(
                "smearing_kind",
                "smearing_parameter",
                "number_of_eigenvalues_kpoints",
            ),
            forbidden_values=[None]
        )
        # Functional settings
        hash_dict['functional_long_name'] = method.functional_long_name

        # Basis set settings
        hash_dict['settings_basis_set'] = settings_basis_set

        # k-point sampling settings if present. Add number of kpoints as
        # detected from eigenvalues. TODO: we would like to have info on the
        # _reducible_ k-point-mesh:
        #    - grid dimensions (e.g. [ 4, 4, 8 ])
        #    - or list of reducible k-points
        if method.smearing_kind is not None:
            hash_dict['smearing_kind'] = method.smearing_kind
        if method.smearing_parameter is not None:
            smearing_parameter = '%.4f' % (method.smearing_parameter * J_to_Ry)
            hash_dict['smearing_parameter'] = smearing_parameter
        try:
            scc = self.backend[s_scc][-1]
            eigenvalues = scc['eigenvalues']
            kpt = eigenvalues[-1]['eigenvalues_kpoints']
        except (KeyError, IndexError):
            pass
        else:
            hash_dict['number_of_eigenvalues_kpoints'] = str(len(kpt))

        # SCF convergence settings
        conv_thr = repr_method.get('scf_threshold_energy_change', None)
        if conv_thr is not None:
            conv_thr = '%.13f' % (conv_thr * J_to_Ry)
            hash_dict['scf_threshold_energy_change'] = conv_thr

        return hash_dict

    def group_parametervariation_hash_dict(self, method: Method, settings_basis_set: RestrictedDict, repr_method: Section):
        """Dictionary containing the parameters used for convergence test
        grouping
        This is the source for generating the related hash."""
        param_dict = RestrictedDict(
            mandatory_keys=(
                "functional_long_name",
                "scf_threshold_energy_change",
            ),
            optional_keys=(
                "atoms_pseudopotentials",
            ),
            forbidden_values=[None]
        )

        # TODO: Add other DFT-specific properties
        # considered variations:
        #   - smearing kind/width
        #   - k point grids
        #   - basis set parameters
        # convergence threshold should be kept constant during convtest
        param_dict['functional_long_name'] = method.functional_long_name
        conv_thr = repr_method.get('scf_threshold_energy_change', None)
        if conv_thr is not None:
            conv_thr = '%.13f' % (conv_thr * J_to_Ry)
        param_dict['scf_threshold_energy_change'] = conv_thr

        # Pseudopotentials are kept constant, if applicable
        if settings_basis_set is not None:
            pseudos = settings_basis_set.get('atoms_pseudopotentials', None)
            if pseudos is not None:
                param_dict['atoms_pseudopotentials'] = pseudos

        return param_dict

    def create_xc_functional_shortname(self, xc_longname):
        """Use lookup table to transform xc functional long- into shortname.
        """
        # Loof for "special" functional names listed in table
        """Easily editable table of 'short' XC functional names"""
        xc_functional_shortname = {
            'HF_X': 'HF',
            'HYB_GGA_XC_B3LYP5': 'hybrid-GGA',
            'HYB_GGA_XC_HSE06': 'hybrid-GGA',
            'BEEF-vdW': 'vdW-DF'
        }
        shortname = xc_functional_shortname.get(xc_longname, None)

        # If not, look into other options:
        if shortname is None:
            xc_functional_starts = {
                "LDA": "LDA",
                "GGA": "GGA",
                "HYB_GGA": "hybrid-GGA",
                "MGGA": "meta-GGA",
                "HYB_MGGA": "hybrid-meta-GGA",
                "HF": "HF"
            }
            sections = xc_longname.split("+")
            # decompose long name, this could be done more consistent with the
            # composition of the long name
            funcnames = []
            for section in sections:
                funcname = section.split('*')[-1]
                for func_start in xc_functional_starts:
                    if funcname.startswith(func_start):
                        funcnames.append(func_start)
                        break
            funcnames = set(funcnames)

            # Only one functional is defined
            # (usually for correlation and exchange)
            if len(funcnames) == 1:
                shortname = xc_functional_starts[func_start]
            # Two functionals that give a hybrid-GGA functional
            elif "GGA" in funcnames and "HF" in funcnames:
                shortname = "hybrid-GGA"

        if shortname is None:
            self.logger.info(
                "Could not find a functional shortname for xc_functional {}."
                .format(xc_longname)
            )

        return shortname

    def normalize(self, ctx: Context) -> None:
        # Fetch resources
        repr_method = ctx.representative_method
        repr_system = ctx.representative_system
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        method = sec_enc.method
        material = sec_enc.material
        settings_basis_set = get_basis_set_settings(ctx, self.backend, self.logger)

        # Fill metainfo
        self.basis_set_type(method)
        self.code_name(method)
        self.code_version(method)
        self.core_electron_treatment(method)
        self.functional_long_name(method, repr_method)
        self.functional_type(method)
        self.smearing_kind(method, repr_method)
        self.smearing_parameter(method, repr_method)
        self.method_hash(method, settings_basis_set, repr_method)
        self.group_eos_hash(method, material, repr_method)
        self.group_parametervariation_hash(method, settings_basis_set, repr_system, repr_method)


class MethodGWNormalizer(MethodDFTNormalizer):
    """A base class that is used for processing GW calculations.
    """
    def gw_starting_point(self, method: Method, repr_method: Section) -> None:
        try:
            ref = repr_method["section_method_to_method_refs"][0]
            method_to_method_kind = ref["method_to_method_kind"]
            method_to_method_ref = ref["method_to_method_ref"]
        except KeyError:
            pass
        else:
            if method_to_method_kind == "starting_point":
                methods = self.backend[s_method]
                start_method = methods[method_to_method_ref]
                xc_functional = MethodDFTNormalizer.functional_long_name_from_method(start_method, methods)
                method.gw_starting_point = xc_functional

    def functional_type(self, method: Method) -> None:
        method.functional_type = "GW"

    def gw_type(self, method: Method, repr_method: Section) -> None:
        method.gw_type = repr_method["electronic_structure_method"]

    def normalize(self, ctx: Context) -> None:
        # Fetch resources
        repr_method = ctx.representative_method
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        method = sec_enc.method

        # Fill metainfo
        self.code_name(method)
        self.code_version(method)
        self.functional_type(method)
        self.gw_type(method, ctx.representative_method)
        self.gw_starting_point(method, repr_method)


class PropertiesNormalizer():
    """A base class that is used for processing calculated quantities that
    should be extracted to Encyclopedia.
    """
    def __init__(self, backend: LocalBackend, logger):
        self.backend = backend
        self.logger = logger

    def get_reciprocal_cell(self, band_structure: ElectronicBandStructure, prim_atoms: Atoms, orig_atoms: Atoms):
        """A reciprocal cell for this calculation. If the original unit cell is
        not a primitive one, then we will use the one given by spglib.

        If the used code is calculating it's own primitive cell, then the
        reciprocal cell used in e.g. band structure calculation might differ
        from the one given by spglib. To overcome this we would need to get the
        exact reciprocal cell used by the calculation or then test for
        different primitive cells
        (https://atztogo.github.io/spglib/definition.html#transformation-to-the-primitive-cell)
        whether the k-point path stays inside the first Brillouin zone.

        Args:
            prim_atoms: The primitive system as detected by MatID.
            orig_atoms: The original simulation cell.

        Returns:
            Reciprocal cell as 3x3 numpy array. Units are in SI (1/m).
        """
        primitive_cell = prim_atoms.get_cell()
        source_cell = orig_atoms.get_cell()

        volume_primitive = primitive_cell.volume
        volume_source = source_cell.volume
        volume_diff = abs(volume_primitive - volume_source)

        if volume_diff > (0.001)**3:
            recip_cell = primitive_cell.reciprocal() * 1e10
        else:
            recip_cell = source_cell.reciprocal() * 1e10

        band_structure.reciprocal_cell = recip_cell

    def get_band_gaps(self, band_structure: ElectronicBandStructure, energies: np.array, path: np.array) -> None:
        """Given the band structure and fermi level, calculates the band gap
        for spin channels and also reports the total band gap as the minum gap
        found.
        """
        # Handle spin channels separately to find gaps for spin up and down
        reciprocal_cell = band_structure.reciprocal_cell
        fermi_level = band_structure.fermi_level
        n_channels = energies.shape[0]

        gaps: List[BandGap] = [None, None]
        for channel in range(n_channels):
            channel_energies = energies[channel, :, :]
            lower_defined = False
            upper_defined = False
            num_bands = channel_energies.shape[0]
            band_indices = np.arange(num_bands)
            band_minima_idx = channel_energies.argmin(axis=1)
            band_maxima_idx = channel_energies.argmax(axis=1)
            band_minima = channel_energies[band_indices, band_minima_idx]
            band_maxima = channel_energies[band_indices, band_maxima_idx]

            # Add a tolerance to minima and maxima
            band_minima_tol = band_minima + config.normalize.fermi_level_precision
            band_maxima_tol = band_maxima - config.normalize.fermi_level_precision

            for band_idx in range(num_bands):
                band_min = band_minima[band_idx]
                band_max = band_maxima[band_idx]
                band_min_tol = band_minima_tol[band_idx]
                band_max_tol = band_maxima_tol[band_idx]

                # If any of the bands band crosses the Fermi level, there is no
                # band gap
                if band_min_tol <= fermi_level and band_max_tol >= fermi_level:
                    break
                # Whole band below Fermi level, save the current highest
                # occupied band point
                elif band_min_tol <= fermi_level and band_max_tol <= fermi_level:
                    gap_lower_energy = band_max
                    gap_lower_idx = band_maxima_idx[band_idx]
                    lower_defined = True
                # Whole band above Fermi level, save the current lowest
                # unoccupied band point
                elif band_min_tol >= fermi_level:
                    gap_upper_energy = band_min
                    gap_upper_idx = band_minima_idx[band_idx]
                    upper_defined = True
                    break

            # If a highest point of the valence band and a lowest point of the
            # conduction band are found, and the difference between them is
            # positive, save the information location and value.
            if lower_defined and upper_defined and gap_upper_energy - gap_lower_energy >= 0:

                # See if the gap is direct or indirect by comparing the k-point
                # locations with some tolerance
                k_point_lower = path[gap_lower_idx]
                k_point_upper = path[gap_upper_idx]
                k_point_distance = self.get_k_space_distance(reciprocal_cell, k_point_lower, k_point_upper)
                is_direct_gap = k_point_distance <= config.normalize.k_space_precision

                gap = BandGap()
                gap.type = "direct" if is_direct_gap else "indirect"
                gap.value = float(gap_upper_energy - gap_lower_energy)
                gap.conduction_band_min_k_point = k_point_upper
                gap.conduction_band_min_energy = float(gap_upper_energy)
                gap.valence_band_max_k_point = k_point_lower
                gap.valence_band_max_energy = float(gap_lower_energy)
                gaps[channel] = gap

        # For unpolarized calculations we simply report the gap if it is found.
        if n_channels == 1:
            if gaps[0] is not None:
                band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[0])
        # For polarized calculations we report the gap separately for both
        # channels. Also we report the smaller gap as the the total gap for the
        # calculation.
        elif n_channels == 2:
            if gaps[0] is not None:
                band_structure.m_add_sub_section(ElectronicBandStructure.band_gap_spin_up, gaps[0])
            if gaps[1] is not None:
                band_structure.m_add_sub_section(ElectronicBandStructure.band_gap_spin_down, gaps[1])
            if gaps[0] is not None and gaps[1] is not None:
                if gaps[0].value <= gaps[1].value:
                    band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[0])
                else:
                    band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[1])
            else:
                if gaps[0] is not None:
                    band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[0])
                elif gaps[1] is not None:
                    band_structure.m_add_sub_section(ElectronicBandStructure.band_gap, gaps[1])

    def get_k_space_distance(self, reciprocal_cell: np.array, point1: np.array, point2: np.array) -> float:
        """Used to calculate the Euclidean distance of two points in k-space,
        given relative positions in the reciprocal cell.

        Args:
            reciprocal_cell: Reciprocal cell.
            point1: The first position in k-space.
            point2: The second position in k-space.

        Returns:
            float: Euclidian distance of the two points in k-space in SI units.
        """
        k_point_displacement = np.dot(reciprocal_cell, point1 - point2)
        k_point_distance = np.linalg.norm(k_point_displacement)

        return k_point_distance

    def get_brillouin_zone(self, band_structure: ElectronicBandStructure) -> None:
        """Returns a dictionary containing the information needed to display
        the Brillouin zone for this material. This functionality could be put
        into the GUI directly, with the Brillouin zone construction performed
        from the reciprocal cell.

        The Brillouin Zone is a Wigner-Seitz cell, and is thus uniquely
        defined. It's shape does not depend on the used primitive cell.
        """
        recip_cell = band_structure.reciprocal_cell.magnitude
        brillouin_zone = structure.get_brillouin_zone(recip_cell)
        bz_json = json.dumps(brillouin_zone)
        band_structure.brillouin_zone = bz_json

    def band_structure(self, properties: Properties, run_type: str, material_type: str, representative_scc: Section, sec_system: Section) -> None:
        """Band structure data following arbitrary path.

        Currently this function is only taking into account the normalized band
        structures, and if they have a k-point that is labeled as '?', the
        whole band strcture will be ignored.
        """
        # Band structure data is extracted only from bulk calculations. This is
        # because for now we need the reciprocal cell of the primitive cell
        # that is only available from the symmetry analysis. Once the
        # reciprocal cell is directly reported with the band structure this
        # restriction can go away.
        if run_type != RunType.run_type.type.single_point or material_type != Material.material_type.type.bulk:
            return

        orig_atoms = sec_system.tmp["representative_atoms"]
        symmetry_analyzer = sec_system["section_symmetry"][0].tmp["symmetry_analyzer"]
        prim_atoms = symmetry_analyzer.get_primitive_system()

        # Try to find an SCC with band structure data, give priority to
        # normalized data
        for src_name in ["section_k_band_normalized", "section_k_band"]:
            bands = representative_scc.get(src_name)
            if bands is None:
                continue
            norm = "_normalized" if src_name == "section_k_band_normalized" else ""

            # Loop over bands
            for band_data in bands:
                band_structure = ElectronicBandStructure()
                kpoints = []
                energies = []
                try:
                    segments = band_data['section_k_band_segment' + norm]
                except KeyError:
                    return

                # Loop over segments
                for segment_src in segments:
                    try:
                        seg_k_points = segment_src["band_k_points" + norm]
                        seg_energies = segment_src["band_energies" + norm]
                        seg_labels = segment_src['band_segm_labels' + norm]
                    except KeyError:
                        return
                    if "?" in seg_labels:
                        return

                    seg_energies = np.swapaxes(seg_energies, 1, 2)
                    kpoints.append(seg_k_points)
                    energies.append(seg_energies)

                # Continue to calculate band gaps and other properties.
                kpoints = np.concatenate(kpoints, axis=0)
                energies = np.concatenate(energies, axis=2)
                self.get_reciprocal_cell(band_structure, prim_atoms, orig_atoms)
                self.get_brillouin_zone(band_structure)

                # If we are using a normalized band structure (or the Fermi level
                # is defined somehow else), we can add band gap information
                fermi_level = None
                if src_name == "section_k_band_normalized":
                    fermi_level = 0.0
                if fermi_level is not None:
                    band_structure.fermi_level = fermi_level
                    self.get_band_gaps(band_structure, energies, kpoints)

                properties.m_add_sub_section(Properties.electronic_band_structure, band_structure)

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

    def energies(self, properties: Properties, gcd: int, representative_scc: Section) -> None:
        energy_dict = {}
        if representative_scc is not None:
            energies_entries = {
                "energy_total": "Total E",
                "energy_total_T0": "Total E projected to T=0",
                "energy_free": "Free E",
            }
            for energy_name, label in energies_entries.items():
                result = representative_scc.get(energy_name)
                if result is not None:
                    energy_dict[label] = result / gcd

            if len(energy_dict) == 0:
                energy_dict = None
        energies = json.dumps(energy_dict)
        properties.energies = energies

    def normalize(self, ctx: Context) -> None:
        # There needs to be a valid SCC in order to extract any properties
        representative_scc = ctx.representative_scc
        if representative_scc is None:
            return

        # Fetch resources
        sec_enc = self.backend.get_mi2_section(Encyclopedia.m_def)
        properties = sec_enc.properties
        properties.scc_index = int(ctx.representative_scc_idx)
        run_type = ctx.run_type
        material_type = ctx.material_type
        sec_system = ctx.representative_system
        gcd = ctx.greatest_common_divisor

        # Save metainfo
        self.band_structure(properties, run_type, material_type, representative_scc, sec_system)
        self.energies(properties, gcd, representative_scc)
