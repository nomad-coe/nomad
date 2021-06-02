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


from typing import Dict, Any, List
import os
import logging
import numpy as np
import re
from ase.data import chemical_symbols

from nomad.datamodel import EntryArchive
from nomad.parsing import FairdiParser
from nomad.parsing.file_parser import TextParser, Quantity
from nomad.datamodel.metainfo.common_dft import Run, Method, System, SingleConfigurationCalculation


class BasicParser(FairdiParser):
    '''
    Defines a fairdi parser that parse basic quantities for sections method, system and
    single_configuration_calculation.

    Arguments:
        specifications: dictionary that will be passed on to  FairdiParser
        units_mapping: dictionary of nomad units for basic quantities such as length
        auxiliary_files: re pattern to match auxilliary files from mainfile. If no files
            are found will match files in working directory.
        kwargs: metainfo_key: re pattern pairs used to parse quantity
    '''
    def __init__(self, specifications: Dict[str, Any], **kwargs):
        super().__init__(**specifications)
        self.specifications = specifications
        self.units_mapping = kwargs.get('units_mapping', {})
        self.auxilliary_files = kwargs.get('auxilliary_files', '')
        self.mainfile_parser = TextParser()
        for key, pattern in kwargs.items():
            if isinstance(pattern, str):
                self.mainfile_parser._quantities.append(
                    Quantity(key, pattern, repeats=True, flatten=False))
            elif isinstance(pattern, tuple) and isinstance(pattern[0], str):
                self.mainfile_parser._quantities.append(
                    Quantity(key, pattern[0], str_operation=pattern[1], repeats=True))
        self._re_float = r'\-*\d+\.\d+E*e*\-*\+*\d*'
        self.auxilliary_parsers: List[TextParser] = []

    def init_parser(self):
        '''
        Initializes the mainfile and auxiliary parsers.
        '''
        self.mainfile_parser.mainfile = self.mainfile
        self.mainfile_parser.logger = self.logger

        auxilliary_files = self.mainfile_parser.get('auxilliary_files', os.listdir(self.maindir))
        # remove duplicates, maintain order
        auxilliary_files = [f for n, f in enumerate(auxilliary_files) if f not in auxilliary_files[:n]]
        self.auxilliary_parsers = []
        for filename in auxilliary_files:
            filename = os.path.basename(filename)
            if self.mainfile_parser.get('auxilliary_files') is None:
                if not self.auxilliary_files or not re.match(self.auxilliary_files, filename):
                    continue
            filename = os.path.join(self.maindir, filename)
            if not os.path.isfile(filename):
                continue
            parser = self.mainfile_parser.copy()
            parser.mainfile = filename
            parser.logger = self.logger
            self.auxilliary_parsers.append(parser)

    def parse(self, mainfile: str, archive: EntryArchive, logger=None) -> None:
        '''
        Triggers parsing of mainfile and writing parsed quantities to archive.
        '''
        self.mainfile = os.path.abspath(mainfile)
        self.maindir = os.path.dirname(self.mainfile)
        self.archive = archive
        self.logger = logger if logger is not None else logging

        self.init_parser()

        def set_value(section, key, value, unit=None, shape=None, dtype=None):
            dtype = dtype if dtype is not None else type(value)
            if value is None:
                return
            try:
                if not hasattr(value, 'units'):
                    value = np.reshape(np.array(
                        value, dtype=np.dtype(dtype)), shape) if shape is not None else dtype(value)
                    value = value * unit if unit is not None else value
                setattr(section, key, value)
            except Exception:
                pass

        def get_value(source, pattern, key=None):
            if isinstance(source, str):
                val = re.findall(pattern, source)
                return val[0] if len(val) == 1 else val
            elif isinstance(source, list):
                return [get_value(s, pattern) for s in source]
            elif isinstance(source, dict):
                return source.get(key)
            else:
                return source

        def remove_empty_section(sections, definition):
            for n in range(len(sections) - 1, -1, -1):
                empty = True
                for _ in sections[n].m_traverse():
                    empty = False
                    break
                if empty:
                    sec_run.m_remove_sub_section(definition, n)

        sec_run = self.archive.m_create(Run)
        sec_run.program_name = self.specifications.get('code_name', '')

        energy_unit = self.units_mapping.get('energy', 1.0)
        length_unit = self.units_mapping.get('length', 1.0)
        mass_unit = self.units_mapping.get('mass')
        time_unit = self.units_mapping.get('time')

        re_f = r'\-*\d+\.\d+E*e*\-*\+*\d*'

        for key, values in self.mainfile_parser.items():
            if values is None:
                # get if from auxiliary files
                values = []
                for parser in self.auxilliary_parsers:
                    values.extend(parser.get(key, []))
            if values is None or len(values) == 0:
                continue
            # set header quantities
            set_value(sec_run, key, values[0])
            for n, value in enumerate(values):
                # method related quantities
                if len(sec_run.section_method) <= n:
                    sec_run.m_create(Method)
                sec_method = sec_run.section_method[n]
                set_value(sec_method, key, value)

                # system related quantities
                if len(sec_run.section_system) <= n:
                    sec_run.m_create(System)
                sec_system = sec_run.section_system[n]
                set_value(sec_system, key, value)

                # calculation related quantities
                if len(sec_run.section_single_configuration_calculation) <= n:
                    sec_run.m_create(SingleConfigurationCalculation)
                sec_scc = sec_run.section_single_configuration_calculation[n]
                set_value(sec_scc, key, value)

                # specific quantities that need formatting
                if 'energy' in key:
                    shape = np.shape(value) if key == 'energy_reference_fermi' else None
                    set_value(sec_scc, key, value, energy_unit, shape, np.float64)

                if 'atom_forces' in key:
                    val = get_value(value, rf'.*({re_f}) +({re_f}) +({re_f}).*', 'atom_forces')
                    if mass_unit is not None and time_unit is not None:
                        unit = mass_unit * length_unit / time_unit ** 2
                    else:
                        unit = energy_unit / length_unit
                    set_value(sec_scc, 'atom_forces', val, unit, (np.size(val) // 3, 3), np.float64)

                if 'lattice_vectors' in key:
                    val = get_value(value, rf'({re_f}) +({re_f}) +({re_f}).*', 'lattice_vectors')
                    set_value(sec_system, 'lattice_vectors', val, length_unit, (3, 3), np.float64)
                    if val is not None:
                        sec_system.configuration_periodic_dimensions = [True, True, True]

                if 'atom_positions' in key:
                    sub_key = 'atom_positions_scaled' if 'atom_positions_scaled' in key else 'atom_positions'
                    val = get_value(value, rf'({re_f}) +({re_f}) +({re_f}).*', sub_key)
                    unit = length_unit
                    if sub_key == 'atom_positions_scaled':
                        try:
                            val = np.dot(np.array(val, dtype=np.dtype(np.float64)), sec_system.lattice_vectors.magnitude)
                            unit = 1.0
                        except Exception:
                            pass
                    set_value(sec_system, 'atom_positions', val, unit, (np.size(val) // 3, 3), np.float64)

                if 'atom_velocities' in key:
                    val = get_value(value, rf'({re_f}) +({re_f}) +({re_f}).*', 'atom_velocities')
                    set_value(sec_system, 'atom_velocities', val, length_unit / time_unit, (np.size(val) // 3, 3), np.float64)

                if 'atom_labels' in key:
                    val = get_value(value, r'([A-Z][a-z]*)\s', 'atom_labels')
                    val = [val] if isinstance(val, str) else val
                    set_value(sec_system, 'atom_labels', val, shape=(len(val)), dtype=str)

                if 'atom_atom_number' in key:
                    val = get_value(value, r'(\d+)\s', 'atom_atom_number')
                    val = [val] if isinstance(val, str) else val
                    set_value(sec_system, 'atom_atom_number', val, shape=(len(val)), dtype=np.int32)
                    set_value(sec_system, 'atom_labels', [chemical_symbols[int(n)] for n in sec_system.atom_atom_number], shape=(len(val)))

        # remove unfilled sections
        remove_empty_section(sec_run.section_method, Run.section_method)
        remove_empty_section(sec_run.section_system, Run.section_system)
        remove_empty_section(sec_run.section_single_configuration_calculation, Run.section_single_configuration_calculation)
