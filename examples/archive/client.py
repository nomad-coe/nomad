'''
A simple example that uses the NOMAD client library to access the archive.
'''

from nomad.client import ArchiveQuery
from nomad.metainfo import units


query = ArchiveQuery(
    # url='http://nomad-lab.eu/prod/rae/beta/api',
    query={
        '$and': {
            'dft.code_name': 'VASP',
            '$not': {
                'atoms': ["Ti", "O"]
            }
        }
    },
    required={
        'section_run': {
            'section_single_configuration_calculation': {
                'energy_total': '*'
            },
            'section_system': '*'
        }
    },
    per_page=10,
    max=None)


print(query)

# for i, result in enumerate(query):
#     if i < 10:
#         calc = result.section_run[0].section_single_configuration_calculation[-1]
#         energy = calc.energy_total
#         formula = calc.single_configuration_calculation_to_system_ref.chemical_composition_reduced
#         print('%s: energy %s' % (formula, energy.to(units.hartree)))
