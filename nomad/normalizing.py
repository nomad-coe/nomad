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

from abc import ABCMeta, abstractmethod
import ase
import numpy
import spglib

from nomadcore.json_support import addShasOfJson
from statsnormalizer.stats import crystalSystem
from statsnormalizer.classify_structure import ClassifyStructure

from nomad.parsing import AbstractParserBackend

"""
After parsing calculations have to be normalized with a set of *normalizers*.
In NOMAD-coe those were programmed in python (we'll reuse) and scala (we'll rewrite).
"""


class Normalizer(metaclass=ABCMeta):
    """
    A base class for normalizers. Normalizers work on a :class:`AbstractParserBackend` instance
    for read and write.
    """
    def __init__(self, backend: AbstractParserBackend) -> None:
        self._backend = backend

    @abstractmethod
    def normalize(self) -> None:
        pass


class SystemNomalizer(Normalizer):
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

    def _normalize_section_system(self, g_index: int) -> None:
        backend = self._backend
        results = dict()

        atom_labels = backend.get_value('atom_labels', g_index)
        atom_species = backend.get_value('atom_atom_numbers', g_index)
        if atom_labels is not None and atom_species is None:
            atom_label_to_num = SystemNomalizer.atom_label_to_num
            atom_species = [atom_label_to_num(atom_label) for atom_label in atom_labels]

        periodic_dirs = backend.get_value('configuration_periodic_dimensions', g_index)
        formula = None
        if atom_species:
            results['atom_species'] = atom_species
            atom_symbols = [ase.data.chemical_symbols[atom_number] for atom_number in atom_species]
            formula = ase.Atoms(atom_symbols).get_chemical_formula(mode='all')
            formula_reduced = ase.Atoms(atom_symbols).get_chemical_formula(mode='reduce')
            if periodic_dirs is not None and any(periodic_dirs):
                formula_bulk = formula_reduced
            else:
                formula_bulk = formula

        cell = backend.get_value('simulation_cell', g_index)
        if cell is not None:
            results['lattice_vectors'] = cell

        positions = backend.get_value('atom_positions', g_index)
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

            results['gIndex'] = g_index
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

        backend.addValue("configuration_raw_gid", configuration_id, g_index)

        backend.openNonOverlappingSection('section_topology')
        backend.addValue("atom_species", atom_species)
        backend.closeNonOverlappingSection('section_topology')

        if symm:
            # TODO: check what is wrong, the commented meta names seem not to exist
            #       in the current meta info
            # for quantity in ["number", "international", "hall", "choice", "pointgroup"]:
            #     v = symm.get(quantity)
            #     if v is not None:
            #         backend.addValue("spacegroup_3D_" + quantity, v, g_index)

            # for quantity in ["transformation_matrix"]:
            #     v = symm.get(quantity)
            #     if v is not None:
            #         backend.addArrayValues(
            #             "spacegroup_3D_" + quantity, numpy.asarray(v), g_index)

            n = symm.get("number")
            if n:
                backend.openNonOverlappingSection('section_symmetry')
                backend.addValue("bravais_lattice", crystalSystem(n), g_index)
                backend.closeNonOverlappingSection('section_symmetry')

            # backend.addValue("chemical_composition", formula, g_index)
            # backend.addValue("chemical_composition_reduced", formula_reduced, g_index)
            # backend.addValue("chemical_composition_bulk_reduced", formula_bulk, g_index)

            # for quantity in ["origin_shift", "std_lattice"]:
            #     v = symm.get(quantity)
            #     if v is not None:
            #         backend.addArrayValues(
            #             "spacegroup_3D_" + quantity, 1.0e-10 * numpy.asarray(v, dtype=float),
            #             g_index)

            # for (r, t) in zip(symm.get("rotations", []), symm.get("translations", [])):
            #     backend.openNonOverlappingSection("section_spacegroup_3D_operation")
            #     backend.addArrayValues("spacegroup_3D_rotation", numpy.asarray(r), g_index)
            #     backend.addArrayValues(
            #         "spacegroup_3D_translation", 1.0e-10 * numpy.asarray(t, dtype=float),
            #         g_index)
            #     backend.closeNonOverlappingSection("section_spacegroup_3D_operation")

            # v = symm.get("wyckoffs")
            # if v is not None:
            #     for w in v:
            #         backend.addValue("spacegroup_3D_wyckoff", w, g_index)

    def normalize(self) -> None:
        for g_index in self._backend.get_sections('section_system'):
            self._normalize_section_system(g_index)

normalizers = [SystemNomalizer]
