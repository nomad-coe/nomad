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

import names
import random
from essential_generators import DocumentGenerator
import datetime
from ase.data import chemical_symbols

from nomad import datamodel, parsing, utils

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
files = ['/'.join(gen.url().split('/')[3:]) for _ in range(0, number_of)]

low_numbers_for_atoms = [1, 1, 2, 2, 2, 2, 2, 3, 3, 4]
low_numbers_for_refs_and_datasets = [0, 0, 0, 0, 1, 1, 1, 2]


def _gen_user():
    id, first, last, email = random.choice(users)
    return utils.POPO(id=id, first_name=first, last_name=last, email=email)


def _gen_dataset():
    id, name = random.choice(datasets)
    return utils.POPO(id=id, name=name, doi=_gen_ref())


def _gen_ref():
    id, value = random.choice(references)
    return utils.POPO(id=id, value=value)


def generate_calc(pid: int = 0) -> datamodel.CalcWithMetadata:
    random.seed(pid)

    self = datamodel.CalcWithMetadata()

    self.upload_id = utils.create_uuid()
    self.calc_id = utils.create_uuid()

    self.upload_time = datetime.datetime.now()
    self.calc_hash = utils.create_uuid()
    self.pid = pid
    self.mainfile = random.choice(files)
    self.files = list([self.mainfile] + random.choices(files, k=random.choice(low_numbers_for_atoms)))
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
    self.spacegroup = '1'
    self.code_name = random.choice(codes)
    self.code_version = '1.0.0'

    return self


if __name__ == '__main__':
    import time
    n = 2
    start = time.time()
    for pid in range(0, n):
        calc = generate_calc(pid)
        print(calc.to_dict())

    print('%f' % ((time.time() - start) / n))
