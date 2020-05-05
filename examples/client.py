from nomad import config
from nomad.client import query_archive
from nomad.metainfo import units

# this will not be necessary, once this is the official NOMAD version
config.client.url = 'https://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/testing-major/api'

query = ArchiveQuery(
    query={
        'dft.compound_type': 'binary',
        'dft.crystal_system': 'cubic',
        'dft.code_name': 'FHI-aims',
        'atoms': ['O']
    },
    required={
        'section_run': {
            'section_single_configuration_calculation[0]': {
                'energy_total': '*'
            }
        }
    },
    per_page=10,
    max=1000)

print(query)

for result in query[0:10]:
    energy = result.section_run[0].section_single_configuration_calculation[0].energy_total
    print('Energy %s' % energy.to(units.hartree))
