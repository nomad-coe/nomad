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
import datetime
from ase.data import chemical_symbols
from ase.spacegroup import Spacegroup

from nomad import datamodel, utils, files
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
    id, name = random.choice(datasets)
    id_str = str(id)
    if datamodel.Dataset.m_def.a_mongo.objects(dataset_id=id_str).first() is None:
        datamodel.Dataset(
            user_id=random.choice(users), dataset_id=id_str, name=name,
            doi=_gen_ref().value).a_mongo.create()
    return id_str


def _gen_ref():
    return random.choice(references)


def generate_calc(pid: int = 0, calc_id: str = None, upload_id: str = None) -> datamodel.EntryMetadata:
    random.seed(pid)

    entry = datamodel.EntryMetadata()

    entry.upload_id = upload_id if upload_id is not None else utils.create_uuid()
    entry.calc_id = calc_id if calc_id is not None else utils.create_uuid()

    entry.upload_time = datetime.datetime.utcnow()
    entry.calc_hash = utils.create_uuid()
    entry.pid = str(pid)
    entry.mainfile = random.choice(filepaths)
    entry.files = list([entry.mainfile] + random.choices(filepaths, k=random.choice(low_numbers_for_files)))
    entry.uploader = _gen_user()

    entry.with_embargo = random.choice([True, False])
    entry.published = True
    entry.coauthors = list(_gen_user() for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))
    entry.shared_with = list(_gen_user() for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))
    entry.comment = random.choice(comments)
    entry.references = list(_gen_ref() for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))
    entry.datasets = list(
        _gen_dataset()
        for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))

    entry.atoms = list(random.choices(chemical_symbols[1:], k=random.choice(low_numbers_for_atoms)))
    entry.formula = ''.join('%s%d' % (atom, random.choice(low_numbers_for_atoms)) for atom in entry.atoms)
    entry.formula = entry.formula.replace('1', '')

    dft_metadata = entry.m_create(datamodel.DFTMetadata)
    dft_metadata.basis_set = random.choice(basis_sets)
    dft_metadata.xc_functional = random.choice(xc_functionals)
    dft_metadata.system = random.choice(systems)
    dft_metadata.crystal_system = random.choice(crystal_systems)
    spacegroup = random.randint(1, 225)
    dft_metadata.spacegroup = str(spacegroup)
    dft_metadata.spacegroup_symbol = Spacegroup(spacegroup).symbol
    dft_metadata.code_name = random.choice(codes)
    dft_metadata.code_version = '1.0.0'

    dft_metadata.n_total_energies = random.choice(range(0, 5))
    dft_metadata.geometries = ['%d' % random.randint(1, 500), '%d' % random.randint(1, 500)]

    return entry


def test_common_metainfo():
    from nomad.datamodel.metainfo.run.run import Run
    from nomad.datamodel.metainfo.run.system import System, Atoms

    run = Run()
    system = run.m_create(System)
    system.atoms = Atoms(labels=['H', 'H', 'O'])

    assert run.system[0].atoms.labels == ['H', 'H', 'O']


def test_vasp_metainfo():
    from nomad.datamodel.metainfo.run import Run
    from vaspparser.metainfo import m_env  # pylint: disable=unused-import
    run = Run()
    assert 'vasp_src_date' in run.m_def.all_quantities


if __name__ == '__main__':
    import sys
    from elasticsearch.helpers import bulk

    from nomad import infrastructure

    print('Generate test data and add it to search and files')
    print('  first arg is number of calcs (code runs)')
    print('  second arg is number uploads to spread calcs over')

    infrastructure.setup_mongo()
    infrastructure.setup_elastic()

    n_calcs, n_uploads = int(sys.argv[1]), int(sys.argv[2])
    pid = 1

    for calcs_per_upload in utils.chunks(range(0, n_calcs), int(n_calcs / n_uploads)):
        upload_id = utils.create_uuid()
        upload_files = files.StagingUploadFiles(
            upload_id=upload_id, create=True, is_authorized=lambda: True)

        search_entries = []
        calcs = []
        for _ in calcs_per_upload:
            calc = generate_calc(pid, upload_id=upload_id)
            assert calc.upload_id == upload_id
            calc.published = True

            for filepath in calc.files:
                if len(filepath) > 0:
                    with upload_files.raw_file(filepath, 'wt') as f:
                        f.write('this is a generated test file')

            upload_files.write_archive(calc.calc_id, {
                'section_run': [{'test': 'this is a generated test files'}],
                'processing_logs': [{'event': 'this is a generated test file'}]
            })

            search_entry = calc.a_elastic.create_index_entry()
            search_entry.n_total_energies = random.choice(low_numbers_for_total_energies)
            search_entry.n_geometries = low_numbers_for_geometries
            for _ in range(0, random.choice(search_entry.n_geometries)):
                search_entry.geometries.append(utils.create_uuid())
            search_entries.append(search_entry)

            pid += 1
            calcs.append(calc)

        bulk(
            infrastructure.elastic_client,
            [entry.to_dict(include_meta=True) for entry in search_entries])

        upload_files.pack(calcs)
        upload_files.delete()
