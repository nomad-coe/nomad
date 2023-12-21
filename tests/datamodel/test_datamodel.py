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

"""
A generator for random test calculations.
"""
import random
from essential_generators import DocumentGenerator

from nomad.datamodel import EntryArchive, EntryMetadata
from nomad.metainfo import MSection, Quantity, SubSection
from nomad.parsing.parsers import parser_dict

number_of = 20

random.seed(0)
gen = DocumentGenerator()

users = [
    '20bb9766-d338-4314-be43-7906042a5086',
    'a03af8b6-3aa7-428a-b3b1-4a6317e576b6',
    '54cb1f64-f84e-4815-9ade-440ce0b5430f',
]
basis_sets = ['Numeric AOs', 'Gaussians', '(L)APW+lo', 'Plane waves']
xc_functionals = ['LDA', 'GGA', 'hybrid', 'meta-GGA', 'GW', 'unknown']
crystal_systems = [
    'triclinic',
    'monoclinic',
    'orthorombic',
    'tetragonal',
    'hexagonal',
    'cubic',
]
systems = ['atom', 'molecule/cluster', '2D/surface', 'bulk']
comments = [gen.sentence() for _ in range(0, number_of)]
references = [(i + 1, gen.url()) for i in range(0, number_of)]
datasets = [(i + 1, gen.slug()) for i in range(0, number_of)]
codes = list(
    set(
        [
            parser.code_name
            for parser in parser_dict.values()
            if hasattr(parser, 'code_name')
        ]
    )
)  # type: ignore
filepaths = ['/'.join(gen.url().split('/')[3:]) for _ in range(0, number_of)]

low_numbers_for_atoms = [1, 1, 2, 2, 2, 2, 2, 3, 3, 4]
low_numbers_for_files = [1, 2, 2, 3, 3, 3, 3, 3, 4, 4]
low_numbers_for_refs_and_datasets = [0, 0, 0, 0, 1, 1, 1, 2]
low_numbers_for_total_energies = [1, 2, 2, 2, 3, 4, 5, 6, 10, 100]
low_numbers_for_geometries = [1, 2, 2, 3, 3, 4, 4]


def test_common_metainfo():
    from nomad.datamodel.metainfo.simulation.run import Run
    from nomad.datamodel.metainfo.simulation.system import System, Atoms

    run = Run()
    system = run.m_create(System)
    system.atoms = Atoms(labels=['H', 'H', 'O'])

    assert run.system[0].atoms.labels == ['H', 'H', 'O']


def test_vasp_metainfo():
    from nomad.datamodel.metainfo.simulation.run import Run
    from electronicparsers.vasp.metainfo import m_env  # pylint: disable=unused-import

    run = Run()
    assert 'vasp_src_date' in run.m_def.all_quantities
