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

from nomadcore.json_support import addShasOfJson
from nomad.normalizing.normalizer import SystemBasedNormalizer


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
            if periodic_dirs is not None and any(periodic_dirs):
                formula_bulk = formula_reduced
            else:
                formula_bulk = formula

        cell = section_system.get('simulation_cell', None)
        if cell is not None:
            results['lattice_vectors'] = cell

        positions = section_system['atom_positions']
        if positions is not None:
            results['atom_positions'] = positions
            if not formula:
                formula = 'X%d' % len(positions) if len(positions) != 1 else 'X'

        if periodic_dirs is not None:
            results['configuration_periodic_dimensions'] = periodic_dirs.tolist()
        # TODO: @dts, might be good to clean this up so it is more readable in the
        # future.
        configuration_id = 's' + addShasOfJson(results).b64digests()[0][0:28]

        self._backend.addValue("configuration_raw_gid", configuration_id)
        self._backend.addValue("atom_species", atom_species)
        self._backend.addValue("chemical_composition", formula)
        self._backend.addValue("chemical_composition_reduced", formula_reduced)
        self._backend.addValue("chemical_composition_bulk_reduced", formula_bulk)

    # def map_matid_to_nomad_system_types(self, system_type):
    #     """ We map the system type classification from matid to Nomad values.

    #     Args:
    #         system_type: This is
    #     Returns:

