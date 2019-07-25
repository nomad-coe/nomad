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

"""
A generator for random test calculations.
"""

import names
import random
from essential_generators import DocumentGenerator
import datetime
from ase.data import chemical_symbols
from ase.spacegroup import Spacegroup

from nomad import datamodel, parsing, utils, files

number_of = 20

random.seed(0)
gen = DocumentGenerator()

users = [(i + 1, names.get_first_name(), names.get_last_name(), gen.email()) for i in range(0, number_of)]
basis_sets = ['Numeric AOs', 'Gaussians', '(L)APW+lo', 'Plane waves']
xc_functionals = ['LDA', 'GGA', 'hybrid', 'meta-GGA', 'GW', 'unknown']
crystal_systems = ['triclinic', 'monoclinic', 'orthorombic', 'tetragonal', 'hexagonal', 'cubic']
systems = ['atom', 'molecule/cluster', '2D/surface', 'bulk']
comments = [gen.sentence() for _ in range(0, number_of)]
references = [(i + 1, gen.url()) for i in range(0, number_of)]
datasets = [(i + 1, gen.slug()) for i in range(0, number_of)]
codes = [parser[8:] for parser in parsing.parser_dict.keys()]
filepaths = ['/'.join(gen.url().split('/')[3:]) for _ in range(0, number_of)]

low_numbers_for_atoms = [1, 1, 2, 2, 2, 2, 2, 3, 3, 4]
low_numbers_for_files = [1, 2, 2, 3, 3, 3, 3, 3, 4, 4]
low_numbers_for_refs_and_datasets = [0, 0, 0, 0, 1, 1, 1, 2]
low_numbers_for_total_energies = [1, 2, 2, 2, 3, 4, 5, 6, 10, 100]
low_numbers_for_geometries = [1, 2, 2, 3, 3, 4, 4]


def _gen_user():
    id, first, last, email = random.choice(users)
    return utils.POPO(id=id, first_name=first, last_name=last, email=email)


def _gen_dataset():
    id, name = random.choice(datasets)
    return utils.POPO(id=id, name=name, doi=_gen_ref())


def _gen_ref():
    id, value = random.choice(references)
    return utils.POPO(id=id, value=value)


def generate_calc(pid: int = 0, calc_id: str = None, upload_id: str = None) -> datamodel.CalcWithMetadata:
    random.seed(pid)

    self = datamodel.DFTCalcWithMetadata()

    self.upload_id = upload_id if upload_id is not None else utils.create_uuid()
    self.calc_id = calc_id if calc_id is not None else utils.create_uuid()

    self.upload_time = datetime.datetime.utcnow()
    self.calc_hash = utils.create_uuid()
    self.pid = pid
    self.mainfile = random.choice(filepaths)
    self.files = list([self.mainfile] + random.choices(filepaths, k=random.choice(low_numbers_for_files)))
    self.uploader = _gen_user()

    self.with_embargo = random.choice([True, False])
    self.published = True
    self.coauthors = list(_gen_user() for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))
    self.shared_with = list(_gen_user() for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))
    self.comment = random.choice(comments)
    self.references = list(_gen_ref() for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))
    self.datasets = list(
        _gen_dataset()
        for _ in range(0, random.choice(low_numbers_for_refs_and_datasets)))

    self.atoms = list(random.choices(chemical_symbols[1:], k=random.choice(low_numbers_for_atoms)))
    self.formula = ''.join('%s%d' % (atom, random.choice(low_numbers_for_atoms)) for atom in self.atoms)
    self.formula = self.formula.replace('1', '')

    self.basis_set = random.choice(basis_sets)
    self.xc_functional = random.choice(xc_functionals)
    self.system = random.choice(systems)
    self.crystal_system = random.choice(crystal_systems)
    spacegroup = random.randint(1, 225)
    self.spacegroup = str(spacegroup)
    self.spacegroup_symbol = Spacegroup(spacegroup).symbol
    self.code_name = random.choice(codes)
    self.code_version = '1.0.0'

    return self


if __name__ == '__main__':
    import sys
    import json
    from elasticsearch.helpers import bulk

    from nomad import infrastructure, search

    print('Generate test data and add it to search and files')
    print('  first arg is number of calcs (code runs)')
    print('  second arg is number uploads to spread calcs over')

    infrastructure.setup_logging()
    infrastructure.setup_elastic()

    n_calcs, n_uploads = int(sys.argv[1]), int(sys.argv[2])
    pid = 1

    for calcs_per_upload in utils.chunks(range(0, n_calcs), int(n_calcs / n_uploads)):
        upload_id = utils.create_uuid()
        upload = datamodel.UploadWithMetadata(upload_id=upload_id)
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
            with upload_files.archive_file(calc.calc_id, 'wt') as f:
                f.write(json.dumps({'section_run': [{'test': 'this is a generated test files'}]}))
            with upload_files.archive_log_file(calc.calc_id, 'wt') as f:
                f.write('this is a generated test file')

            search_entry = search.Entry.from_calc_with_metadata(calc)
            search_entry.n_total_energies = random.choice(low_numbers_for_total_energies)
            search_entry.n_geometries = low_numbers_for_geometries
            for _ in range(0, random.choice(search_entry.n_geometries)):
                search_entry.geometries.append(utils.create_uuid())
            search_entries.append(search_entry)

            pid += 1
            calcs.append(calc)

        upload.calcs = calcs

        bulk(
            infrastructure.elastic_client,
            [entry.to_dict(include_meta=True) for entry in search_entries])

        upload_files.pack(upload)
        upload_files.delete()
