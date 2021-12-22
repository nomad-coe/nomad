'''
A simple example used in the NOMAD webinar API tutorial
'''

from nomad.client import ArchiveQuery
from nomad.metainfo import units

query = ArchiveQuery(
    query={
        'results.method.simulation.program_name': 'VASP',
        'results.material.elements': ['Ti', 'O'],
        'results.method.simulation.geometry_optimization': {
            'convergence_tolerance_energy_difference:lt': 1e-22,
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
    },
    parallel=10,
    max=100)

print(query)

for result in query:
    calc = result.workflow[0].calculation_result_ref
    formula = calc.system_ref.chemical_composition_reduced
    total_energy = calc.energy.total.value.to(units.eV)
    print(f'{formula}: {total_energy}')
