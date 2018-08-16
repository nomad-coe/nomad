import numpy as np
import os.path
import logging
import sys
import setup_paths
from ase import Atoms
from systax import Class3DAnalyzer
from systax.exceptions import SystaxError
import systax.geometry

from nomad.normalizing import Normalizer

class SymmetryNormalizer(Normalizer):
    """
    This is basically a copy of the legace NOMAD-coe symmetry normalizer.
    """

    def normalize(self) -> None:
        for g_index in self._backend.get_sections('section_system'):
            section_system = dict_stream.readNextDict()
            if section_system is None:
                break

            try:
                # TODO maybe write a reasonable backend first?
                atom_labels = self._backend.get_value('atom_labels', g_index)
                atom_pos = self._backend.get_value('atom_positions', g_index)

                # Try to first read the cell information from the renamed metainfo
                try:
                    cell = section_system["lattice_vectors"]
                except KeyError:
                    cell = section_system["simulation_cell"]

                pbc = section_system["configuration_periodic_dimensions"]
                pbc = pbc[0]

            # If these values could not be read, then skip this system
            except KeyError:
                logging.exception(
                    "The necessary information is not available for the system: {}"
                    .format(section_system)
                )
                continue

            # Make the data into numpy arrays
            atom_labels = np.array(atom_labels["flatData"]).reshape(atom_labels["shape"])
            atom_pos = np.array(atom_pos["flatData"]).reshape(atom_pos["shape"])
            cell = np.array(cell["flatData"]).reshape(cell["shape"])
            pbc = np.array(pbc["flatData"]).reshape(pbc["shape"])

            # The pbc should be defined as a single-dimensional list
            if len(pbc.shape) == 2:
                pbc = pbc[0, :]

            # If all dimensions are not defined to be periodic, skip this system as
            # it cannot represent a crystal with symmetries
            if not pbc.all():
                continue
            try:
                atoms = Atoms(
                    positions=1e10*atom_pos,
                    symbols=atom_labels,
                    cell=1e10*cell,
                    pbc=pbc
                )
            # If there is an issue in creating the Atoms object, e.g. because the
            # labels are invalid, then skip this system
            except Exception:
                logging.exception(
                    "Could not create an ASE.Atoms object for the system: {}. Could"
                    " be that the calculation is using customized atomic "
                    "labels.".format(section_system)
                )
                continue

            # Figure out the dimensionality of the system by using the
            # dimensionality detection included in the systax package. If the
            # system is not 3D, then it is skipped.
            try:
                dim, _ = systax.geometry.get_dimensionality(
                    atoms,
                    cluster_threshold=3.0)
            # If this exception is thrown, the dimensionality could not be detected
            # because there are multiple energetically isolated components in the
            # system.
            except SystaxError:
                continue
            if dim != 3:
                continue
            try:
                analyzer = Class3DAnalyzer(atoms)

                space_group_number = analyzer.get_space_group_number()
                hall_number = analyzer.get_hall_number()
                hall_symbol = analyzer.get_hall_symbol()
                international_short = analyzer.get_space_group_international_short()
                point_group = analyzer.get_point_group()
                crystal_system = analyzer.get_crystal_system()
                bravais_lattice = analyzer.get_bravais_lattice()

                conv_sys = analyzer._get_spglib_conventional_system()
                conv_pos = conv_sys.get_scaled_positions()
                conv_cell = conv_sys.get_cell()
                conv_num = conv_sys.get_atomic_numbers()
                conv_wyckoff = analyzer._get_spglib_wyckoff_letters_conventional()
                conv_equivalent_atoms = analyzer._get_spglib_equivalent_atoms_conventional()

                prim_sys = analyzer._get_spglib_primitive_system()
                prim_pos = prim_sys.get_scaled_positions()
                prim_cell = prim_sys.get_cell()
                prim_num = prim_sys.get_atomic_numbers()
                prim_wyckoff = analyzer._get_spglib_wyckoff_letters_primitive()
                prim_equivalent_atoms = analyzer._get_spglib_equivalent_atoms_primitive()

                orig_wyckoff = analyzer.get_wyckoff_letters_original()
                orig_equivalent_atoms = analyzer.get_equivalent_atoms_original()
                transform = analyzer._get_spglib_transformation_matrix()
                origin_shift = analyzer._get_spglib_origin_shift()
            except:
                # If there is an issue getting the symmetry data (happens e.g. when
                # atoms overlap), then skip this system
                logging.exception(
                    "Error in getting the symmetry information for system: {}. This"
                    " can be e.g. caused by overlapping atoms."
                    .format(section_system)
                )
                continue

            # Push the derived values to the backend
            # print(context, file=sys.stderr)
            context = section_system["uri"]
            backend.openContext(context)
            symGid = backend.openSection("section_symmetry")

            backend.addValue("symmetry_method", "spg_normalized")

            backend.addValue("space_group_number", space_group_number)
            backend.addValue("hall_number", hall_number)
            backend.addValue("hall_symbol", hall_symbol)
            backend.addValue("international_short_symbol", international_short)
            backend.addValue("point_group", point_group)
            backend.addValue("crystal_system", crystal_system)
            backend.addValue("bravais_lattice", bravais_lattice)
            backend.addArrayValues("origin_shift", origin_shift)
            backend.addArrayValues("transformation_matrix", transform)

            stdGid = backend.openSection("section_std_system")
            backend.addArrayValues("lattice_vectors_std", conv_cell)
            backend.addArrayValues("atom_positions_std", conv_pos)
            backend.addArrayValues("atomic_numbers_std", conv_num)
            backend.addArrayValues("wyckoff_letters_std", conv_wyckoff)
            backend.addArrayValues("equivalent_atoms_std", conv_equivalent_atoms)
            backend.closeSection("section_std_system", stdGid)

            primGid = backend.openSection("section_primitive_system")
            backend.addArrayValues("lattice_vectors_primitive", prim_cell)
            backend.addArrayValues("atom_positions_primitive", prim_pos)
            backend.addArrayValues("atomic_numbers_primitive", prim_num)
            backend.addArrayValues("wyckoff_letters_primitive", prim_wyckoff)
            backend.addArrayValues("equivalent_atoms_primitive", prim_equivalent_atoms)
            backend.closeSection("section_primitive_system", primGid)

            origGid = backend.openSection("section_original_system")
            backend.addArrayValues("wyckoff_letters_original", orig_wyckoff)
            backend.addArrayValues("equivalent_atoms_original", orig_equivalent_atoms)
            backend.closeSection("section_original_system", origGid)

            backend.closeSection("section_symmetry", symGid)
            backend.closeContext(context)
            sys.stdout.flush()

        backend.finishedParsingSession("ParseSuccess", None)
        sys.stdout.flush()
