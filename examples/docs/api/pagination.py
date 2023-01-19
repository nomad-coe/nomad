import requests

base_url = 'http://nomad-lab.eu/prod/v1/api/v1'
json_body = {
    'query': {
        'results.material.elements': {
            'all': ['Ti', 'O']
        }
    },
    'pagination': {
        'page_size': 10
    },
    'required': {
        'include': ['results.material.chemical_formula_hill']
    }
}

formulas = set()

while len(formulas) < 100:
    response = requests.post(f'{base_url}/entries/query', json=json_body)
    response_json = response.json()

    for data in response_json['data']:
        formulas.add(data['results']['material']['chemical_formula_hill'])

    next_value = response_json['pagination'].get('next_page_after_value')
    if not next_value:
        break
    json_body['pagination']['page_after_value'] = next_value

print(formulas)
