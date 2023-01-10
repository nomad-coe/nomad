'''
Demonstrates how to use requests for a simple query and archive access.
'''

import requests
from ase import Atoms

base_url = 'http://nomad-lab.eu/prod/v1/api/v1'
response = requests.post(
    f'{base_url}/entries/archive/query',
    json={
        'pagination': {
            'page_size': 1
        },
        'required': {
            'run': {
                'system[-1]': {
                    'atoms': '*'
                }
            }
        }
    })
response_json = response.json()
nomad_atoms = response_json['data'][0]['archive']['run'][0]['system'][-1]['atoms']
atoms = Atoms(
    symbols=nomad_atoms['labels'],
    positions=nomad_atoms['positions'],
    cell=nomad_atoms['lattice_vectors'],
    pbc=nomad_atoms['periodic']
)

print(atoms)


from nomad.client import ArchiveQuery

query = ArchiveQuery(
    required={
        'run': {
            'system[-1]': {
                'atoms': '*'
            }
        }
    })

result = query.download(1)[0]
atoms = result.run[0].system[-1].atoms.to_ase()

print(atoms.get_chemical_formula(mode='reduce'))
