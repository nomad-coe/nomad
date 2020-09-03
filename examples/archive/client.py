'''
A simple example that uses the NOMAD client library to access the archive.
'''

from nomad.client import ArchiveQuery
from nomad.metainfo import units


query = ArchiveQuery(
    # url='http://nomad-lab.eu/prod/rae/beta/api',
    query={
        '$and': [
            {'dft.code_name': 'VASP'},
            {'$gte': {'n_atoms': 3}},
            {'$lte': {'dft.workflow.section_relaxation.final_energy_difference': 1e-24}}
        ]
    },
    required={
        'section_workflow': {
            'section_relaxation': {
                'final_calculation_ref': {
                    'energy_total': '*',
                    'single_configuration_calculation_to_system_ref': {
                        'chemical_composition_reduced': '*'
                    }
                }
            }
        }
    },
    parallel=10,
    per_page=1000,
    max=10000)

for i, result in enumerate(query):
    if i < 10:
        calc = result.section_workflow.section_relaxation.final_calculation_ref
        energy = calc.energy_total
        formula = calc.single_configuration_calculation_to_system_ref.chemical_composition_reduced
        print('%s: energy %s' % (formula, energy.to(units.hartree)))

print(query)
