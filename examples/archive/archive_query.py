'''
A simple example used in the NOMAD webinar API tutorial
'''

from nomad.client import ArchiveQuery
from nomad.metainfo import units

query = ArchiveQuery(
    query={
        'results.method.simulation.program_name': 'VASP',
        'results.material.elements': ['Ti', 'O'],
        'results.properties.geometry_optimization': {
            'final_energy_difference:lt': 1e-22,
        }
    },
    required={
        'workflow': {
            'calculation_result_ref': {
                'energy': '*',
                'system_ref': {
                    'chemical_composition_reduced': '*'
                }
            }
        }
    })

for result in query.download(100):
    calc = result.workflow[0].calculation_result_ref
    formula = calc.system_ref.chemical_composition_reduced
    if calc.energy.total:
        total_energy = calc.energy.total.value.to(units.eV)
    else:
        total_energy = 'N/A'
    print(f'{formula}: {total_energy}')
