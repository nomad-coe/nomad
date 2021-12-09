'''
A simple example used in the NOMAD webinar API tutorial
'''

from nomad.client import ArchiveQuery

query = ArchiveQuery(
    url='http://nomad-lab.eu/prod/rae/api',
    query={
        'dft.code_name': 'VASP',
        'atoms': ['Ti', 'O']
    },
    required={
        'section_run': {
            'section_single_configuration_calculation[-1]': {
                'energy_total': '*',
                'section_dos': '*'
            }
        }
    },
    parallel=1,
    max=10)

print(query)
result = query[0]
print(result.section_run[0].section_single_configuration_calculation[-1].section_dos[0].dos_energies)
