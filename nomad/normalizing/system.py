# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import ase
import numpy
import spglib

from nomadcore.json_support import addShasOfJson
from statsnormalizer.stats import crystalSystem
from statsnormalizer.classify_structure import ClassifyStructure

from nomad.normalizing.normalizer import SystemBasedNormalizer

# TODO: check what is wrong, the commented meta names seem not to exist
#       in the current meta info


class SystemNormalizer(SystemBasedNormalizer):
    """
    This normalizer performs all system (atoms, cells, etc.) related normalizations
    of the legacy NOMAD-coe *stats* normalizer.
    """

    @staticmethod
    def atom_label_to_num(atom_label):
        atom_label = atom_label[:3].title()

        for symbol_length in reversed(range(1, 4)):
            symbol = atom_label[:symbol_length]
            if symbol in ase.data.chemical_symbols:
                return ase.data.chemical_symbols.index(symbol)

        return 0

    def normalize_system(self, section_system) -> None:
        results = dict()

        atom_labels = section_system['atom_labels']
        atom_species = section_system['atom_atom_numbers']
        if atom_labels is not None and atom_species is None:
            atom_label_to_num = SystemNormalizer.atom_label_to_num
            atom_species = [atom_label_to_num(atom_label) for atom_label in atom_labels]

        periodic_dirs = section_system['configuration_periodic_dimensions']
        formula = None
        if atom_species:
            results['atom_species'] = atom_species
            atom_symbols = [ase.data.chemical_symbols[atom_number] for atom_number in atom_species]
            formula = ase.Atoms(atom_symbols).get_chemical_formula(mode='all')
            formula_reduced = ase.Atoms(atom_symbols).get_chemical_formula(mode='reduce')
            # if periodic_dirs is not None and any(periodic_dirs):
            #     formula_bulk = formula_reduced
            # else:
            #     formula_bulk = formula

        cell = section_system['simulation_cell']
        if cell is not None:
            results['lattice_vectors'] = cell

        positions = section_system['atom_positions']
        if positions is not None:
            results['atom_positions'] = positions
            if not formula:
                formula = 'X%d' % len(positions) if len(positions) != 1 else 'X'

        if periodic_dirs is not None:
            results['configuration_periodic_dimensions'] = periodic_dirs.tolist()

        configuration_id = 's' + addShasOfJson(results).b64digests()[0][0:28]
        if cell is not None and atom_labels is not None:
            if cell is not None:
                results['simulation_cell'] = cell
            if atom_labels is not None:
                results['atom_labels'] = atom_labels

            results['gIndex'] = section_system['gIndex']
            results['name'] = 'section_system'
            structure = ClassifyStructure(None, jsonValue={
                "sections": [{
                    "name": "section_run",
                    "gIndex": 1,
                    "sections": [results]
                }]
            })
            classification = structure.classify()

            if classification.get('classificationStatus', None) == 'ClassificationSuccess':
                classType = classification['sections'][0]['sections'][0]['structure_kind']
            else:
                classType = 'NoClassification'

            if classType == 'Bulk' and positions is not None and atom_species is not None and cell is not None:
                acell = numpy.asarray(cell) * 1.0e10
                cellInv = numpy.linalg.inv(cell)
                symm = spglib.get_symmetry_dataset(
                    (acell, numpy.dot(positions, cellInv), atom_species),
                    0.002, -1)  # use m instead of Angstrom?
                if symm:
                    symm['configuration_raw_gid'] = configuration_id

        self._backend.addValue("configuration_raw_gid", configuration_id)
        self._backend.addValue("atom_species", atom_species)

        if symm:
            # for quantity in ["number", "international", "hall", "choice", "pointgroup"]:
            #     v = symm.get(quantity)
            #     if v is not None:
            #         self._backend.addValue("spacegroup_3D_" + quantity, v)

            # for quantity in ["transformation_matrix"]:
            #     v = symm.get(quantity)
            #     if v is not None:
            #         self._backend.addArrayValues(
            #             "spacegroup_3D_" + quantity, numpy.asarray(v))

            n = symm.get("number")
            if n:
                self._backend.openNonOverlappingSection('section_symmetry')
                self._backend.addValue("bravais_lattice", crystalSystem(n))
                self._backend.closeNonOverlappingSection('section_symmetry')

            # self._backend.addValue("chemical_composition", formula)
            # self._backend.addValue("chemical_composition_reduced", formula_reduced)
            # self._backend.addValue("chemical_composition_bulk_reduced", formula_bulk)

            # for quantity in ["origin_shift", "std_lattice"]:
            #     v = symm.get(quantity)
            #     if v is not None:
            #         backend.addArrayValues(
            #             "spacegroup_3D_" + quantity, 1.0e-10 * numpy.asarray(v, dtype=float))

            # for (r, t) in zip(symm.get("rotations", []), symm.get("translations", [])):
            #     self._backend.openNonOverlappingSection("section_spacegroup_3D_operation")
            #     self._backend.addArrayValues("spacegroup_3D_rotation", numpy.asarray(r))
            #     self._backend.addArrayValues(
            #         "spacegroup_3D_translation", 1.0e-10 * numpy.asarray(t, dtype=float))
            #     self._backend.closeNonOverlappingSection("section_spacegroup_3D_operation")

            # v = symm.get("wyckoffs")
            # if v is not None:
            #     for w in v:
            #         self._backend.addValue("spacegroup_3D_wyckoff", w)
