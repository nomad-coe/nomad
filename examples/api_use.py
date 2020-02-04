"""
This is a brief example on how to use the public nomad@FAIRDI API.
"""

from bravado.client import SwaggerClient

nomad_url = 'http://repository.nomad-coe.eu/app/api'

# create the bravado client
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url)

# simple search request to print number of public entries
print(client.repo.search().response().result.pagination.total)
