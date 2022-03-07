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

from nomad import infrastructure
from nomad.utils import create_uuid
from nomad.utils.exampledata import ExampleData


def search():
    '''
    Used to construct an API state that is suitable for several different kinds of
    search tests. Constructs an upload with multiple entries containing a wide
    variety of data.
    '''
    infrastructure.setup()
    main_author = infrastructure.keycloak.get_user(username='test')
    data = ExampleData(main_author=main_author)
    upload_id = create_uuid()
    data.create_upload(upload_id=upload_id, published=True, embargo_length=0)

    material_h2o = {
        'structural_type': 'molecule / cluster',
        'elements': ['C', 'H'],
        'chemical_formula_hill': 'CH3',
        'chemical_formula_anonymous': 'AB2',
        'chemical_formula_descriptive': 'CH3',
        'chemical_formula_reduced': 'CH3',
    }
    material_high_entropy_alloy = {
        'structural_type': 'bulk',
        'elements': ['Hf', 'Nb', 'Ta', 'Ti', 'Zr'],
        'chemical_formula_hill': 'HfNbTaTiZr',
        'chemical_formula_anonymous': 'ABCDE',
        'chemical_formula_descriptive': 'HfNbTaTiZr',
        'chemical_formula_reduced': 'HfNbTaTiZr',
    }
    material_graphene = {
        'structural_type': '2D',
        'elements': ['C'],
        'chemical_formula_hill': 'C',
        'chemical_formula_anonymous': 'A',
        'chemical_formula_descriptive': 'C',
        'chemical_formula_reduced': 'C',
    }
    method_dft_vasp = {
        'simulation': {
            'program_name': 'VASP',
            'dft': {
                'xc_functional_names': ['GGA_X_PBE_SOL', 'GGA_C_PBE_SOL']
            }
        }
    }
    method_dft_exciting = {
        'simulation': {
            'program_name': 'exciting',
            'dft': {
                'xc_functional_names': ['LDA_X_PZ', 'LDA_C_PZ']
            }
        }
    }

    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': material_h2o,
            'method': method_dft_vasp,
            'properties': {}
        }
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': material_graphene,
            'method': method_dft_exciting,
            'properties': {}
        }
    )
    data.create_entry(
        upload_id=upload_id,
        entry_id=create_uuid(),
        mainfile='upload/archive.json',
        results={
            'material': material_high_entropy_alloy,
            'method': method_dft_vasp,
            'properties': {}
        }
    )

    data.save()
