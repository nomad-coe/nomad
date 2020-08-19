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

'''
Parser for creating artificial test, brenchmark, and demonstration data.
'''

import json
import os.path
import numpy as np
import random
from ase.data import chemical_symbols
import numpy
import sys
import time
import os
import signal

from nomad import metainfo
from nomad.datamodel.metainfo import m_env as general_nomad_metainfo_env

from .legacy import Backend
from .parser import Parser, MatchingParser


class ArtificalParser(Parser):
    ''' Base class for artifical parsers based on VASP metainfo. '''
    def __init__(self):
        super().__init__()
        self.backend = None

    def init_backend(self):
        self.backend = Backend(metainfo='vasp')


class EmptyParser(MatchingParser):
    '''
    Implementation that produces an empty code_run
    '''
    name = "parsers/empty"

    def run(self, mainfile: str, logger=None) -> Backend:
        backend = Backend(metainfo=general_nomad_metainfo_env, domain=self.domain, logger=logger)
        backend.openSection('section_run')
        backend.addValue('program_name', self.code_name)
        backend.closeSection('section_run', 0)
        return backend


class TemplateParser(ArtificalParser):
    '''
    A parser that generates data based on a template given via the
    mainfile. The template is basically some archive json. Only
    '''
    name = 'parsers/template'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        from nomad.datamodel.metainfo import m_env as metainfo_env
        self._metainfo_env = metainfo_env

    def is_mainfile(
            self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
            compression: str = None) -> bool:
        return filename.endswith('template.json')

    def transform_value(self, name, value):
        ''' allow subclasses to modify values '''
        return value

    def transform_section(self, name, section):
        ''' allow subclasses to modify sections '''
        return section

    def add_section(self, section):
        name = section['_name']
        index = self.backend.openSection(name)

        for key, value in section.items():
            if key.startswith('x_') or key.startswith('_'):
                continue

            if key.startswith('section_'):
                values = value if isinstance(value, list) else [value]
                for value in values:
                    self.add_section(self.transform_section(key, value))
            else:
                value = self.transform_value(key, value)
                if isinstance(value, list):
                    quantity_def = self.backend.env.resolve_definition(key, metainfo.Quantity)
                    if quantity_def.is_scalar:
                        for single_value in value:
                            self.backend.addValue(key, single_value, index)
                    else:
                        self.backend.addArrayValues(key, np.asarray(value), index)
                else:
                    self.backend.addValue(key, value, index)

        self.backend.closeSection(name, index)

    def run(self, mainfile: str, logger=None) -> Backend:
        # tell tests about received logger
        if logger is not None:
            logger.debug('received logger')

        self.init_backend()

        if 'warning' in mainfile:
            self.backend.pwarn('A test warning.')

        template_json = json.load(open(mainfile, 'r'))
        self.add_section(template_json['section_run'][0])
        if 'section_workflow' in template_json:
            self.add_section(template_json['section_workflow'])
        self.backend.finishedParsingSession('ParseSuccess', [])
        logger.debug('a test log entry')
        return self.backend


class ChaosParser(ArtificalParser):
    '''
    Parser that emulates typical error situations. Files can contain a json string (or
    object with key `chaos`) with one of the following string values:
    - exit
    - deadlock
    - consume_ram
    - exception
    - segfault
    - random
    '''
    name = 'parsers/chaos'

    def is_mainfile(
            self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
            compression: str = None) -> bool:
        return filename.endswith('chaos.json')

    def run(self, mainfile: str, logger=None) -> Backend:
        self.init_backend()

        chaos_json = json.load(open(mainfile, 'r'))
        if isinstance(chaos_json, str):
            chaos = chaos_json
        elif isinstance(chaos_json, dict):
            chaos = chaos_json.get('chaos', None)
        else:
            chaos = None

        if chaos == 'random':
            chaos = random.choice(['exit', 'deadlock', 'consume_ram', 'exception', 'segfault'])

        if chaos == 'exit':
            sys.exit(1)
        elif chaos == 'deadlock':
            while True:
                time.sleep(1)
        elif chaos == 'consume_ram':
            data = []
            i = 0
            while True:
                data.append('a' * 10**6)
                i += 1
                logger.info('ate %d mb' % i)
        elif chaos == 'exception':
            raise Exception('Some chaos happened, muhuha...')
        elif chaos == 'segfault':
            os.kill(os.getpid(), signal.SIGSEGV)

        raise Exception('Unknown chaos')


class GenerateRandomParser(TemplateParser):
    name = 'parsers/random'

    basis_set_types = [
        'Numeric AOs', 'Gaussians', '(L)APW+lo', 'FLAPW', 'Plane waves',
        'Real-space grid', 'Local-orbital minimum-basis']

    electronic_structure_methods = [
        'DFT', 'DFT+U', 'full-CI', 'CIS', 'CISD'
        'CCS', 'CCS(D)', 'CCSD', 'CCSD(T)', 'CCSDT(Q)', 'MP2', 'MP3', 'MP4', 'MP5', 'MP6',
        'G0W0', 'scGW', 'LDA', 'hybrid', 'CASPT2', 'MRCIS', 'MRCISD', 'RAS-CI']

    XC_functional_names = [
        'LDA_X', 'LDA_C', 'GGA_X', 'GGA_C', 'HYB_GGA_XC', 'MGGA_X', 'MGGA_C', 'HYB_MGGA_XC']

    low_numbers = [1, 1, 2, 2, 2, 2, 2, 3, 3, 4]

    def __init__(self):
        super(GenerateRandomParser, self).__init__()
        file_dir = os.path.dirname(os.path.abspath(__file__))
        relative_template_file = "random_template.json"
        template_file = os.path.normpath(os.path.join(file_dir, relative_template_file))
        self.template = json.load(open(template_file, 'r'))
        self.random = None

    def is_mainfile(
            self, filename: str, mime: str, buffer: bytes, decoded_buffer: str,
            compression: str = None) -> bool:
        return os.path.basename(filename).startswith('random_')

    def transform_section(self, name, section):
        ''' allow subclasses to modify sections '''
        if name == 'section_system':
            atoms = []
            atom_positions = []
            # different atoms
            for _ in range(0, random.choice(GenerateRandomParser.low_numbers)):
                # ase data starts with a place holder value X, we don't want that
                atom = random.choice(chemical_symbols[1:])
                # atoms of the same number
                for _ in range(0, random.choice(GenerateRandomParser.low_numbers)):
                    atoms.append(atom)
                    atom_positions.append([.0, .0, .0])

            section['atom_labels'] = atoms
            section['atom_positions'] = atom_positions
            return section
        else:
            return section

    def transform_value(self, name, value):
        if name == 'program_basis_set_type':
            return random.choice(GenerateRandomParser.basis_set_types)
        elif name == 'electronic_structure_method':
            return random.choice(GenerateRandomParser.electronic_structure_methods)
        elif name == 'XC_functional_name':
            return random.choice(GenerateRandomParser.XC_functional_names)
        elif name == 'atom_positions':
            return [numpy.random.uniform(0e-10, 1e-9, 3) for _ in value]
        else:
            return value

    def run(self, mainfile: str, logger=None) -> Backend:
        # tell tests about received logger
        if logger is not None:
            logger.debug('received logger')

        self.init_backend()
        seed = int(os.path.basename(mainfile).split('_')[1])
        random.seed(seed)
        numpy.random.seed(seed)
        section = self.template['section_run'][0]
        self.add_section(section)
        self.backend.finishedParsingSession('ParseSuccess', [])
        return self.backend
