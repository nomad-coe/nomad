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

from nomad.normalizing.normalizer import Normalizer, s_scc, s_system, s_frame_sequence, r_frame_sequence_to_sampling, s_sampling_method, r_frame_sequence_local_frames
from nomad.metainfo.encyclopedia import Encyclopedia, Material, Calculation
from nomad.normalizing import symmetry
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

    # NOTE: Enc specific visualization
    def get_atom_labels(self) -> None:
        pass

    # NOTE: Enc specific visualization
    def get_atom_positions(self) -> None:
        pass

    # NOTE: System normalizer
    def get_atomic_density(self) -> None:
        pass

    # NOTE: Enc specific visualization
    def get_atomistic_structure(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def get_band_gap(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def get_band_gap_position(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def get_band_gap_type(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def get_band_structure(self) -> None:
        pass

    # NOTE: Method normalizer
    def get_basis_set_short_name(self) -> None:
        pass

    # NOTE: Method normalizer
    def get_basis_set_type(self) -> None:
        pass

    # NOTE: System normalizer
    def get_bravais_lattice(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def get_brillouin_zone(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def get_brillouin_zone_viewer(self) -> None:
        pass

    def get_calculation(self) -> None:
        pass

    def get_calculation_pid(self) -> None:
        pass

    # NOTE: Enc specific visualization
    def get_cell_angles(self) -> None:
        pass

    # NOTE: System normalizer
    def get_cell_normalized(self) -> None:
        pass

    # NOTE: System normalizer
    def get_cell_primitive(self) -> None:
        pass

    # NOTE: System normalizer
    def get_cell_volume(self) -> None:
        pass

    # NOTE: Parser
    def get_code_name(self) -> None:
        pass

    # NOTE: Parser
    def get_code_version(self) -> None:
        pass

    # NOTE: Repo
    def get_contributor_first_name(self) -> None:
        pass

    # NOTE: Repo
    def get_contributor_last_name(self) -> None:
        pass

    # NOTE: Repo
    def get_contributor_type(self) -> None:
        pass

    # NOTE: Repo
    def get_contributors(self) -> None:
        pass

    # NOTE: Method normalizer
    def get_core_electron_treatment(self) -> None:
        pass

    # NOTE: System normalizer
    def get_crystal_system(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def get_dos(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def get_elastic_constants_matrix(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def get_elastic_deformation_energies(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def get_elastic_fitting_parameters(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def get_elastic_moduli(self) -> None:
        pass

    # NOTE: Elastic properties normalizer
    def get_elastic_properties(self) -> None:
        pass

    def get_energies(self) -> None:
        pass

    # NOTE: Band structure normalizer
    def get_fermi_surface(self) -> None:
        pass

    # NOTE: System normalizer
    def get_formula(self) -> None:
        pass

    # NOTE: System normalizer
    def get_formula_cell(self) -> None:
        pass

    # NOTE: System normalizer
    def get_formula_reduced(self) -> None:
        pass

    # NOTE: System normalizer
    def get_free_wyckoff_parameters(self) -> None:
        pass

    # NOTE: Method normalizer
    def get_functional_long_name(self) -> None:
        pass

    # NOTE: Method normalizer
    def get_functional_type(self) -> None:
        pass

    # TODO: ??
    def get_group_e_min(self) -> None:
        pass

    # TODO: ??
    def get_group_type(self) -> None:
        pass

    # TODO: Method normalizer
    def get_gw_starting_point(self) -> None:
        pass

    # TODO: Method normalizer
    def get_gw_type(self) -> None:
        pass

    # NOTE: Enc specific
    def get_has_bs(self) -> None:
        pass

    # NOTE: Enc specific
    def get_has_dos(self) -> None:
        pass

    # NOTE: Enc specific
    def get_has_fermi_surface(self) -> None:
        pass

    # NOTE: Enc specific
    def get_has_thermal_properties(self) -> None:
        pass

    def get_helmholtz_free_energy(self) -> None:
        pass

    def get_k_point_grid_description(self) -> None:
        pass

    def get_lattice_parameters(self) -> None:
        pass

    def get_mainfile_uri(self) -> None:
        pass

    # NOTE: System normalizer
    def get_mass_density(self) -> None:
        pass

    def get_material_hash(self, material: Material, section_system: Dict) -> None:
        # Get symmetry information from the section
        section_symmetry = section_system["section_symmetry"][0]
        space_group_number = section_symmetry["space_group_number"]
        section_std_system = section_symmetry["section_std_system"][0]
        wyckoff_sets = section_std_system.tmp["wyckoff_sets"]

        # Create and store hash based on SHA512
        norm_hash_string = symmetry.create_symmetry_string(space_group_number, wyckoff_sets)
        material.material_hash = sha512(norm_hash_string.encode('utf-8')).hexdigest()

    def get_material_name(self) -> None:
        pass

    # NOTE: System normalizer
    def get_number_of_atoms(self) -> None:
        pass

    def get_number_of_calculation(self) -> None:
        pass

    # NOTE: System normalizer
    def get_periodic_dimensions(self) -> None:
        pass

    # NOTE: Phonon normalizer
    def get_phonon_dispersion(self) -> None:
        pass

    # NOTE: Phonon normalizer
    def get_phonon_dos(self) -> None:
        pass

    # NOTE: System normalizer
    def get_point_group(self) -> None:
        pass

    # NOTE: Method normalizer
    def get_pseudopotential_type(self) -> None:
        pass

    # NOTE: Repo
    def get_repository_dowload_uri(self) -> None:
        pass

    # NOTE: Repo
    def get_repository_upload_comment(self) -> None:
        pass

    # NOTE: Repo
    def get_repository_uri(self) -> None:
        pass

    # NOTE: Enc specific
    def get_run_type(self, calculation) -> str:
        """Decides what type of calculation this is: single_point, md,
        geometry_optimization, etc.
        """
        run_type = config.services.unavailable_value
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
                run_type = "elastic constants"
            else:
                run_type = "single point"
        # One sequence. Currently calculations with multiple sequences are
        # unsupported.
        elif n_frame_seq == 1:
            frame_seq = frame_sequences[0]
            i_sampling_method = frame_seq[r_frame_sequence_to_sampling]
            section_sampling_method = self._backend[s_sampling_method][i_sampling_method]
            sampling_method = section_sampling_method["sampling_method"]

            if sampling_method == "molecular_dynamics":
                run_type = "molecular dynamics"
            if sampling_method == "geometry_optimization":
                run_type = "geometry optimization"
            if sampling_method == "taylor_expansion":
                run_type = "phonon calculation"

        calculation.run_type = run_type
        return run_type

    def get_scf_threshold(self) -> None:
        pass

    # NOTE: Enc specific
    def get_similar_materials(self) -> None:
        pass

    # NOTE: Method normalizer
    def get_smearing_kind(self) -> None:
        pass

    # NOTE: Method normalizer
    def get_smearing_parameters(self) -> None:
        pass

    # NOTE: System normalizer
    def get_space_group(self) -> None:
        pass

    # NOTE: System normalizer
    def get_space_group_international_short_symbol(self) -> None:
        pass

    # NOTE: System normalizer
    def get_space_group_number(self) -> None:
        pass

    # NOTE: Phonon normalizer
    def get_specific_heat_cv(self) -> None:
        pass

    # NOTE: System normalizer
    def get_springer_classification(self) -> None:
        pass

    # NOTE: System normalizer
    def get_springer_compound_class(self) -> None:
        pass

    # NOTE: System normalizer
    def get_springer_prototype(self) -> None:
        pass

    # NOTE: System normalizer
    def get_structure_prototype(self) -> None:
        pass

    # NOTE: System normalizer
    def get_structure_type(self) -> None:
        pass

    # NOTE: System normalizer
    def get_strukturbericht_designation(self) -> None:
        pass

    # NOTE: System normalizer
    def get_system_type(self) -> None:
        pass

    def get_template(self) -> None:
        pass

    # NOTE: System normalizer
    def get_wyckoff_groups(self) -> None:
        pass

    def normalize(self, logger=None) -> None:
        super().normalize(logger)

        # Initialise metainfo structure
        sec_enc = Encyclopedia()
        material = sec_enc.m_create(Material)
        calculation = sec_enc.m_create(Calculation)

        # Determine run type, stop if unknown
        run_type = self.get_run_type(calculation)
        if run_type == config.services.unavailable_value:
            self.logger.info("unknown run type for encyclopedia")
            return

        # Get the system type, stop if unknown
        if run_type in {"geometry optimization", "molecular dynamics"}:
            frame_seqs = self._backend[s_frame_sequence]
            frame_seq = frame_seqs[0]
            frames = frame_seq[r_frame_sequence_local_frames]
            systems = self._backend[s_system]
            if run_type == "geometry optimization":
                system = systems[frames[-1]]
            elif run_type == "moleculardynamics":
                system = systems[frames[0]]
        elif run_type == "single point":
            system = self._backend[s_system][0]
        system_type = system["system_type"]
        if system_type != "bulk" and system_type != "surface" and system_type != "2D":
            self.logger.info("unknown system type for encyclopedia")
            return

        # Proceed to fill rest of data
        self.get_material_hash(material, system)

        # Put the encyclopedia section into backend
        self._backend.add_mi2_section(sec_enc)
