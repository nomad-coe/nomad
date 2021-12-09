# type: ignore

import requests
import json

response = requests.post(
    'http://nomad-lab.eu/prod/rae/api/v1/entries/archive/query', json={
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
        'pagination': {
            'page_size': 10,
            'page_after_value': '----9KNOtIZc9bDFEWxgjeSRsJrC'
        },
        'required': {
            'section_run': {
                'section_single_configuration_calculation[-1]': {
                    'energy_total': '*'
                },
                'section_system[-1]': {
                    'chemical_composition_bulk_reduced': '*'
                }
            }
        }
    })

print(json.dumps(response.json(), indent=2))
