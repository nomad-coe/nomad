"""
This is a brief example on how use requests with authentication to talks to the NOMAD API.
"""

import requests

from nomad.config import config
from nomad.client import Auth

nomad_url = config.client.url
user = 'yourusername'
password = 'yourpassword'

# create an auth object
auth = Auth(user=user, password=password)

# simple search request to print number of user entries
response = requests.get(f'{nomad_url}/v1/entries', params=dict(owner='user'), auth=auth)
print(response.json()['data'])
