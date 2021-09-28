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

'''
A generator for random test calculations.
'''

import random
from essential_generators import DocumentGenerator

from nomad import datamodel, utils
from nomad.parsing.parsers import parser_dict

number_of = 20

random.seed(0)
gen = DocumentGenerator()

users = ['20bb9766-d338-4314-be43-7906042a5086', 'a03af8b6-3aa7-428a-b3b1-4a6317e576b6', '54cb1f64-f84e-4815-9ade-440ce0b5430f']
basis_sets = ['Numeric AOs', 'Gaussians', '(L)APW+lo', 'Plane waves']
xc_functionals = ['LDA', 'GGA', 'hybrid', 'meta-GGA', 'GW', 'unknown']
crystal_systems = ['triclinic', 'monoclinic', 'orthorombic', 'tetragonal', 'hexagonal', 'cubic']
systems = ['atom', 'molecule/cluster', '2D/surface', 'bulk']
comments = [gen.sentence() for _ in range(0, number_of)]
references = [(i + 1, gen.url()) for i in range(0, number_of)]
datasets = [(i + 1, gen.slug()) for i in range(0, number_of)]
codes = list(set([parser.code_name for parser in parser_dict.values() if hasattr(parser, 'code_name')]))  # type: ignore
filepaths = ['/'.join(gen.url().split('/')[3:]) for _ in range(0, number_of)]

low_numbers_for_atoms = [1, 1, 2, 2, 2, 2, 2, 3, 3, 4]
low_numbers_for_files = [1, 2, 2, 3, 3, 3, 3, 3, 4, 4]
low_numbers_for_refs_and_datasets = [0, 0, 0, 0, 1, 1, 1, 2]
low_numbers_for_total_energies = [1, 2, 2, 2, 3, 4, 5, 6, 10, 100]
low_numbers_for_geometries = [1, 2, 2, 3, 3, 4, 4]


def _gen_user():
    return random.choice(users)


def _gen_dataset():
    id, dataset_name = random.choice(datasets)
    id_str = str(id)
    if datamodel.Dataset.m_def.a_mongo.objects(dataset_id=id_str).first() is None:
        datamodel.Dataset(
            user_id=random.choice(users), dataset_id=id_str, dataset_name=dataset_name,
            doi=_gen_ref().value).a_mongo.create()
    return id_str


def _gen_ref():
    return random.choice(references)


def generate_calc(pid: int = 0, calc_id: str = None, upload_id: str = None, with_embargo=None) -> datamodel.EntryMetadata:
    random.seed(pid)

    entry = datamodel.EntryMetadata()

    entry.upload_id = upload_id if upload_id is not None else utils.create_uuid()
    entry.calc_id = calc_id if calc_id is not None else utils.create_uuid()

    entry.calc_hash = utils.create_uuid()
    entry.pid = str(pid)
    entry.mainfile = random.choice(filepaths)
    entry.files = list([entry.mainfile] + random.choices(filepaths, k=random.choice(low_numbers_for_files)))
    entry.uploader = _gen_user()

    entry.with_embargo = with_embargo if with_embargo is not None else random.choice([True, False])
    entry.published = True
    entry.coauthors = list(_gen_user() for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))
    entry.shared_with = list(_gen_user() for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))
    entry.comment = random.choice(comments)
    entry.references = list(_gen_ref() for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))
    entry.datasets = list(
        _gen_dataset()
        for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))

    return entry


def test_common_metainfo():
    from nomad.datamodel.metainfo.simulation.run import Run
    from nomad.datamodel.metainfo.simulation.system import System, Atoms

    run = Run()
    system = run.m_create(System)
    system.atoms = Atoms(labels=['H', 'H', 'O'])

    assert run.system[0].atoms.labels == ['H', 'H', 'O']


def test_vasp_metainfo():
    from nomad.datamodel.metainfo.simulation.run import Run
    from vaspparser.metainfo import m_env  # pylint: disable=unused-import
    run = Run()
    assert 'vasp_src_date' in run.m_def.all_quantities
