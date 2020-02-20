# Archive API tutorial

This contains the tutorials to use the new archive query functionality.
It uses the new metainfo definition for the archive data. In addition, the archive data
can now be filtered through the new api. The archive are now also stored using a new binary
format msgpack which in principle makes querying faster.

## Archive API
First, we look at how to use the new archive query api. Here we use the python requests
library.
```python
import requests

data = {
    'atoms': 'Fe', 'scroll': True, 'per_page': 10,
    'results': [{"section_run": {"section_single_configuration_calculation[-1]": {"energy_total": None}}}]}
response = requests.post('https://repository.nomad-coe.eu/app/api', data=data)

data = response.json

results = data.get('results')
```
To query the archive, we use the post method where we provide the usual query parameters
in a dictionary. In addition, we provide a schema for the archive data ala graphQL, i.e.
a heirarchical dictionary with null values for each of the property we would like to query.
In the example, we would like to return only the total energy for the last image. It is
important to point out that this schema uses the key 'results' and is a list since
this will be filled with a list of archive data with this schema.

## Archive and the new metainfo
A wrapper for the archive query api is implemented in ArchiveQuery.
```python
from nomad.archive_library.filedb import ArchiveQuery

q = ArchiveQuery(
    atoms=Fe, scroll=True, per_page=10, archive_data={
        "section_run": {"section_single_configuration_calculation[-1]": {"energy_total": None}}})
metainfo = q.query()

for calc in metainfo:
    print(calc.section_run.section_single_configuration_calculation[-1].energy_total)
```
Similarly, we provide query parameters and also the schema which in this case is 'archive_data'.
When we invoke query, a recursive api request is made until all the data matching our
parameters are downloaded. The results are then expressed in the new metainfo scheme
which offers auto-completion feature, among others.

## Msgpack container
The archive data are now stored in a binary format called msgpack. To create a msgpack database
from the archive data and query it, one uses ArchiveFileDB.
```python
from nomad.archive_library.filedb import ArchiveFileDB

db = ArchiveFileDB('archive.msg', mode='w', entry_toc_depth=2)
db.add_data({'calc1':{'secA': {'subsecA': {'propA': 1.0}}, 'secB': {'propB': 'X'}}})
db.add_data({'calc2':{'secA': {'subsecA': {'propA': 2.0}}, 'secB': {'propB': 'Y'}}})
db.close()

db = ArchiveFileDB('archive.msg')
db.query({'calc1':{'secA': None}})
```
In the example, we first create a database in 'archive.msg', and data which are added
will be fragmented down to subsections. We reload it for reading and query all entries
under 'secA' of 'calc1'.
