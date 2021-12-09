"""
This is a brief example on how to use the public nomad@FAIRDI API.
"""

from bravado.client import SwaggerClient

nomad_url = 'http://nomad-lab.eu/prod/rae/api'

# create the bravado client
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url)

# perform the search request to print number of public entries
data = client.repo.search(atoms=['Si', 'O']).response().result
# print the total ammount of search results
print(data.pagination.total)
# print the data of the first result
print(data.results[0])
