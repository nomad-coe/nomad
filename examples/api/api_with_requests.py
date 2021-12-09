# type: ignore

'''
A simple example used in the NOMAD webinar API tutorial
'''

import requests
import json

base_url = 'http://nomad-lab.eu/prod/rae/api'

# response = requests.get(base_url + '/repo?datasets.name=NOMAD%20webinar')
response = requests.get(
    base_url + '/repo',
    params={'datasets.name': 'NOMAD webinar', 'per_page': 1})

data = response.json()
upload_id = data['results'][0]['upload_id']
calc_id = data['results'][0]['calc_id']

response = requests.get(
    base_url + '/archive/%s/%s' % (upload_id, calc_id))

print(json.dumps(response.json(), indent=2))
print(
    response.json()['section_run'][0]['section_single_configuration_calculation'][-1]['section_dos'][0]['dos_energies'])
