# type: ignore

import requests
import json

response = requests.post(
    'http://nomad-lab.eu/prod/rae/api/v1/entries/query', json={
        'query': {
            'and': [
                {
                    'dft.code_name': 'VASP',
                },
                {
                    'not': {
                        'atoms': {
                            'any': ["H", "C", "Li", "Na", "K", "Rb", "Cs"]
                        }
                    }
                }
            ]
        },
        'required': {
            'include': [
                'formula',
                'encyclopedia.material.formula',
                'encyclopedia.material.formula_reduced'
            ]
        },
        'pagination': {
            'page_size': 10,
            'page_after_value': '----9KNOtIZc9bDFEWxgjeSRsJrC'
        }
    })

print(json.dumps(response.json(), indent=2))
