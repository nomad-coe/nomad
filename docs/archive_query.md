# Query Processed Data

The `ArchiveQuery` allows you to search for entries and access their parsed and processed *archive* data
at the same time. Furthermore, all data is accessible through a convenient Python interface
based on the [NOMAD metainfo](schema/python.md#wrap-data-with-python-schema-classes) rather than plain JSON.

## Basic Usage

To define a query, one can, for example, write

```python
from nomad.client.archive import ArchiveQuery

query = ArchiveQuery(query={}, required={}, page_size=10, results_max=10000)
```

Although the above query object has an empty query.

The query object is constructed only. To access the desired data, users need to perform two operations manually.

### Synchronous Interface

#### Fetch

The fetch process is carried out **synchronously**. Users can call

```python
# number_of_entries = query.fetch(1000) # fetch 1000 entries
number_of_entries = query.fetch()  # fetch at most results_max entries
```

to perform the fetch process to fetch up to `results_max` entries. An indicative number `n` can be provided `fetch(n)`
. Given that each upload may contain various numbers of entries, the fetch process guarantees at least `n` entries
will be fetched. The exact number is determined by `page_size`, indicating how many uploads in each page. However, this
would be limited to the `results_max`. The exact qualified number of entries will be returned. Meanwhile, the qualified
upload list would be populated with their IDs. To check all qualified upload IDs, one can call `upload_list()` method to
return the full list.

```python
print(query.upload_list())
```

If applicable, it is possible to fetch a large number of entries first and then perform a second fetch by using some
upload ID in the first fetch result as the `after` argument so that some middle segment can be downloaded.

#### Download

After fetching the qualified uploads, the desired data can be downloaded **asynchronously**. One can call

```python
# results = query.download(1000) # download 1000 entries
results = query.download()  # download all fetched entries if fetched otherwise fetch and download up to `results_max` entries
```

to download up to `results_max` entries. The downloaded results are returned as a list. Alternatively, it is possible to
just download a portion of previously fetched entries at a single time. For example,

```python
# previously fetched for example 1000 entries
# but only download the first 100 (approx.) entries
results = query.download(100)
```

The same `download(n)` method can be called repeatedly. If there are no sufficient entries, new entries will be
automatically fetched. If there are no more entries, the returned result list is empty. For example,

```python
total_results = []
while True:
    result = query.download(100)
    if len(result) == 0:
        break
    total_results.extend(result)
```

There is no retry mechanism in the download process. If any uploads fail to be downloaded due to server error, it is
kept in the list otherwise removed.

### Asynchronous Interface

Some applications, such as Jupyter Notebook, may run a global/top level event loop. To query data in those environments,
one can use the asynchronous interface.

```python
number_of_entries = await query.async_fetch()  # indicative number n applies: async_fetch(n)
results = await query.async_download()  # indicative number n applies: async_download(n)
```

Alternatively, if one wants to use the asynchronous interface, it is necessary to patch the global event loop to allow
nested loops.

To do so, one can add the following at the beginning of the notebook.

```python
import nest_asyncio

nest_asyncio.apply()
```

## A Complete Rundown

Here we show a valid query and acquire data from server.

We first define the desired query and construct the object. We limit the maximum number of entries to be 10000 and 10
uploads per page.

```python
from nomad.client.archive import ArchiveQuery

required = {
    'workflow': {
        'calculation_result_ref': {
            'energy': '*',
            'system_ref': {
                'chemical_composition_reduced': '*'
            }
        }
    }
}

query = {
    'results.method.simulation.program_name': 'VASP',
    'results.material.elements': ['Ti']
}

query = ArchiveQuery(query=query, required=required, page_size=10, results_max=10000)
```

Let's fetch some entries.

```python
query.fetch(1000)
print(query.upload_list())
```

If we print the upload list, it would be

```text
[('-19NlAwxTCCXb6YT9Plifw', 526), ('-2ewONNGTZ68zuTQ6zrRZw', 4), ('-3LrFBvFQtCtmEp3Hy15EA', 12), ('-3ofEqLvSZiqo59vtf-TAQ', 4), ('-4W-jogwReafpdva4ELdrw', 32), ('-BLVfvlJRWawtHyuUWvP_g', 68), ('-Dm30DqRQX6pZUbJYwHUmw', 320), ('-Jfjp-lZSjqyaph2chqZfw', 6), ('-K2QS7s4QiqRg6nMPqzaTw', 82), ('-Li36ZXhQPucJvkd8yzYoA', 10)]
```

So upload `-19NlAwxTCCXb6YT9Plifw` has 526 qualified entries, upload `-2ewONNGTZ68zuTQ6zrRZw` has 4 qualified entries,
and so on. The summation of the above entries gives 1064 entries in total as shown in terminal message.

Now data can be downloaded.

```python
result = query.download(100)
print(f'Downloaded {len(result)} entries.')  # Downloaded 526 entries.
```

Since the first upload has 526 entries, they will be downloaded in this call. The list would have the first upload
removed as it has been downloaded.

```text
[('-2ewONNGTZ68zuTQ6zrRZw', 4), ('-3LrFBvFQtCtmEp3Hy15EA', 12), ('-3ofEqLvSZiqo59vtf-TAQ', 4), ('-4W-jogwReafpdva4ELdrw', 32), ('-BLVfvlJRWawtHyuUWvP_g', 68), ('-Dm30DqRQX6pZUbJYwHUmw', 320), ('-Jfjp-lZSjqyaph2chqZfw', 6), ('-K2QS7s4QiqRg6nMPqzaTw', 82), ('-Li36ZXhQPucJvkd8yzYoA', 10)]
```

It is possible to download more data.

```python
result = query.download(300)
print(f'Downloaded {len(result)} entries.')  # Downloaded 440 entries.
```

The first six uploads will be downloaded to meet the required (at least) 300 entries, which total to 440 entries. What's
left in the list would be

```text
[('-Jfjp-lZSjqyaph2chqZfw', 6), ('-K2QS7s4QiqRg6nMPqzaTw', 82), ('-Li36ZXhQPucJvkd8yzYoA', 10)]
```

We perform one more download call to illustrate that fetch process will be automatically performed.

```python
result = query.download(100)
print(f'Downloaded {len(result)} entries.')  # Downloaded 102 entries.
```

In the above, we request additional 100 entries, however, the list contains only `6+82+10=98` entries, fetch process
will be called to fetch new entries from server. You will see the following message in terminal.

```text
Fetching remote uploads...
787 entries are qualified and added to the download list.
Downloading required data...
Downloaded 102 entries.
[('-NiRWNGjS--JtFoEnYrCfg', 8), ('-OcPUKZtS6u3lXlkWBM4qg', 129), ('-PA35e2ZRsq4AdDfBU4M_g', 14), ('-TG77dGiSTyrDAFNqTKa6Q', 366), ('-VzlPYtnS4q1tSl3NOmlCw', 178), ('-XeqzVqwSMCJFwhvDqWs8A', 14), ('-Y7gwnleQI6Q61jp024fXQ', 16), ('-Zf1RO1MQXegYVTbFybtQQ', 8), ('-Zm4S9VGRdOX-kbF1J5lOA', 50)]
```

## Argument List

The following arguments are acceptable for `ArchiveQuery`.

- `owner` : `str` The scope of data to access. Default: `'visible'`
- `query` : `dict` The API query. There are no validations of any means carried out by the class, users shall make sure
  the provided query is valid. Otherwise, server would return error message.
- `required` : `dict` The required quantities.
- `url` : `str` The database url. It can be the one of your local database. The official NOMAD database is used be
  default if no valid one defined. Default: `http://nomad-lab.eu/prod/v1/api`
- `after` : `str` It can be understood that the data is stored in a sequential list. Each upload has a unique ID,
  if `after` is not provided, the query always starts from the first upload. One can choose to query the uploads in the
  middle of storage by assigning a proper value of `after`.
- `results_max` : `int` Determine how many entries to download. Note each upload may have multiple entries.
- `page_size` : `int` Page size.
- `username` : `str` Username for authentication.
- `password` : `str` Password for authentication.
- `retry` : `int` In the case of server errors, the fetch process is automatically retried every `sleep_time` seconds.
  This argument limits the maximum times of retry.
- `sleep_time` : `float` The interval of fetch retry.