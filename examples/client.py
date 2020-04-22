from nomad import config
from nomad.client import query_archive
from nomad.metainfo import units

# this will not be necessary, once this is the official NOMAD version
config.client.url = 'https://labdev-nomad.esc.rzg.mpg.de/dev/nomad/v0-8-0/api'


aq = query_archive(
    query={
        'upload_id': ['6LUBCju3T3KK3D_fRCJ4qw']
    },
    required={
        'section_run': {
            'section_single_configuration_calculation[0]': {
                'energy_total': '*'
            }
        }
    },
    per_page=100, max=1000)

print('total', aq.total)

for i, e in enumerate(aq):
    if i % 200 == 0:
        print(e.section_run[0].section_single_configuration_calculation[0].energy_total.to(units.hartree))

print(aq)
