'''
This is a brief example on how to use the public API.
'''
import requests

from nomad import config


nomad_url = config.client.url

# perform the search request to print number of public entries
response = requests.post(
    f'{nomad_url}/v1/entries',
    json={
        'query': {
            'results.material.elements:any': ['Si', 'O']
        }
    })
response_data = response.json()

# print the total ammount of search results
print(response_data['pagination']['total'])

# print the data of the first result
print(response_data['data'][0])
