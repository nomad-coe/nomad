# How to use the API

This guide is about using NOMAD's REST APIs directly, e.g. via Python's *request*.

To access the processed data with our client library `nomad-lab` follow
[How to access processed data](archive_query.md). You can also watch our
[video tutorial on the API](../../tutorial/access_api.md#access-data-via-api).

## Different options to use the API

NOMAD offers all its functionality through application
programming interfaces (APIs). More specifically [RESTful HTTP APIs](https://en.wikipedia.org/wiki/Representational_state_transfer){:target="_blank"} that allows you
to use NOMAD as a set of resources (think data) that can be uploaded, accessed, downloaded,
searched for, etc. via [HTTP requests](https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol){:target="_blank"}.

You can get an overview on all NOMAD APIs on the [API page]({{ nomad_url() }}../../gui/analyze/apis).
We will focus here on NOMAD's main API (v1). In fact, this API is also used by
the web interface and should provide everything you need.

There are different tools and libraries to use the NOMAD API that come with different
trade-offs between expressiveness, learning curve, and convenience.

<h4>You can use your browser</h4>

For example to see the metadata for all entries with elements *Ti* and *O* go here:
[{{ nomad_url() }}/v1/entries?elements=Ti&elements=O]({{ nomad_url() }}/v1/entries?elements=Ti&elements=O)

<h4>Use *curl* or *wget*</h4>

REST API's use resources located via URLs. You access URLs with *curl* or *wget*. Same
*Ti*, *O* example as before:
```sh
curl "{{ nomad_url() }}/v1/entries?results.material.elements=Ti&results.material.elements=O" | python -m json.tool
```

<h4> Use Python and requests</h4>

Requests is a popular Python library to use the internets HTTP protocol that is used to
communicate with REST APIs. Install with `pip install requests`.
See [the initial example](#using-request).

<h4>Use our dashboard</h4>

The NOMAD API has an [OpenAPI dashboard]({{ nomad_url() }}/v1). This is an interactive documentation of all
API functions that allows you to try these functions in the browser.

<h4>Use NOMAD's Python package</h4>

Install the [NOMAD Python client library](./pythonlib.md) and use it's `ArchiveQuery`
functionality for a more convenient query based access of archive data following the
[How-to access the processed data](archive_query.md) guide.

## Using request

If you are comfortable with REST APIs and using Pythons `requests` library, this example
demonstrates the basic concepts of NOMAD's main API. You can get more documentation and
details on all functions from the [API dashboard]({{ nomad_url() }}/v1/).

The following issues a search query for all entries that have both *Ti* and *O*
among the elements of their respective materials. It restricts the results to
one entry and only returns the `entry_id`.

```py
import requests
import json

base_url = 'http://nomad-lab.eu/prod/v1/api/v1'

response = requests.post(
    f'{base_url}/entries/query',
    json={
        'query': {
            'results.material.elements': {
                'all': ['Ti', 'O']
            }
        },
        'pagination': {
            'page_size': 1
        },
        'required': {
            'include': ['entry_id']
        }
    })
response_json = response.json()
print(json.dumps(response.json(), indent=2))
```

This will give you something like this:
```json
{
  "owner": "public",
  "query": {
    "name": "results.material.elements",
    "value": {
      "all": [
        "Ti",
        "O"
      ]
    }
  },
  "pagination": {
    "page_size": 1,
    "order_by": "entry_id",
    "order": "asc",
    "total": 17957,
    "next_page_after_value": "--SZVYOxA2jTu_L-mSxefSQFmeyF"
  },
  "required": {
    "include": [
      "entry_id"
    ]
  },
  "data": [
    {
      "entry_id": "--SZVYOxA2jTu_L-mSxefSQFmeyF"
    }
  ]
}
```

The `entry_id` is a unique identifier for, well, entries. You can use it to access
other entry data. For example, you want to access the entry's archive. More
precisely, you want to gather the formula and energies from the main workflow result.
The following requests the archive based on the `entry_id` and only requires some archive sections.

```py
first_entry_id = response_json['data'][0]['entry_id']
response = requests.post(
    f'{base_url}/entries/{first_entry_id}/archive/query',
    json={
        'required': {
            'workflow': {
                'calculation_result_ref': {
                    'energy': '*',
                    'system_ref': {
                        'chemical_composition': '*'
                    }
                }
            }
        }
    })
response_json = response.json()
print(json.dumps(response_json, indent=2))
```

The result will look like this:
```json
{
  "required": {
    "workflow": {
      "calculation_result_ref": {
        "energy": "*",
        "system_ref": {
          "chemical_composition": "*"
        }
      }
    }
  },
  "entry_id": "--SZVYOxA2jTu_L-mSxefSQFmeyF",
  "data": {
    "entry_id": "--SZVYOxA2jTu_L-mSxefSQFmeyF",
    "upload_id": "YXUIZpw5RJyV3LAsFI2MmQ",
    "parser_name": "parsers/fhi-aims",
    "archive": {
      "run": [
        {
          "system": [
            {
              "chemical_composition": "OOSrTiOOOSrTiOOOSrTiOFF"
            }
          ],
          "calculation": [
            {
              "energy": {
                "fermi": -1.1363378335891879e-18,
                "total": {
                  "value": -5.697771591896252e-14
                },
                "correlation": {
                  "value": -5.070133798617076e-17
                },
                "exchange": {
                  "value": -2.3099755059272454e-15
                },
                "xc": {
                  "value": -2.360676843913416e-15
                },
                "xc_potential": {
                  "value": 3.063766944960246e-15
                },
                "free": {
                  "value": -5.697771595558439e-14
                },
                "sum_eigenvalues": {
                  "value": -3.3841806795825544e-14
                },
                "total_t0": {
                  "value": -5.697771593727346e-14
                },
                "correction_entropy": {
                  "value": -1.8310927833270112e-23
                },
                "correction_hartree": {
                  "value": -4.363790430157292e-17
                },
                "correction_xc": {
                  "value": -2.3606768439090564e-15
                }
              },
              "system_ref": "/run/0/system/0"
            }
          ]
        }
      ],
      "workflow": [
        {
          "calculation_result_ref": "/run/0/calculation/0"
        }
      ]
    }
  }
}
```

You can work with the results in the given JSON (or respective Python dict/list) data already.
If you have [NOMAD's Python library](./pythonlib.md) installed ,
you can take the archive data and use the Python interface.
The [Python interface](../plugins/schema_packages.md#wrap-data-with-python-schema-classes) will help with code-completion (e.g. in notebook environments),
resolve archive references (e.g. from workflow to calculation to system), and allow unit conversion:
```py
from nomad.datamodel import EntryArchive
from nomad.metainfo import units

archive = EntryArchive.m_from_dict(response_json['data']['archive'])
result = archive.workflow[0].calculation_result_ref
print(result.system_ref.chemical_composition)
print(result.energy.total.value.to(units('eV')))
```

This will give you an output like this:
```
OOSrTiOOOSrTiOOOSrTiOFF
-355626.93095025205 electron_volt
```

## Different kinds of data

We distinguish between different kinds of NOMAD data and there are different functions in
the API:

- Entry metadata, a summary of extracted data for an entry.
- Raw files, the files as they were uploaded to NOMAD.
- Archive data, all of the extracted data for an entry.

There are also different entities (see also [Datamodel](../../explanation/basics.md)) with different functions in the API:

- Entries
- Uploads
- Datasets
- Users

The API URLs typically start with the entity, followed by the kind of data. Examples
are:

- `entries/query` - Query entries for metadata
- `entries/archive/query` - Query entries for archive data
- `entries/{entry-id}/raw` - Download raw data for a specific entry
- `uploads/{upload-id}/raw/path/to/file` - Download a specific file of an upload

## Common concepts

The [initial example](#getting-started) above, showed how to execute a basic search.
This includes some fundamental concepts that can be applied to many parts of the API.
Let's discuss some of the common concepts.

### Response layout

Functions that have a JSON response, will have a common layout. First, the response will contain all keys and values of the request. The request is not repeated verbatim, but
in a normalized form. Abbreviations in search queries might be expanded, default values for optional parameters are added, or additional response specific information
is included. Second, the response will contain the results under the key `data`.

### Owner

All functions that allow a query will also allow to specify the `owner`. Depending on
the API function, its default value will be mostly `visible`. Some values are only
available if you are [logged in](#authentication).

{{ doc_snippet('owner')}}
### Queries

{{ doc_snippet('query') }}

### Pagination

When you issue a query, usually not all results can be returned. Instead, an API returns
only one *page*. This behavior is controlled through pagination parameters,
like `page_site`, `page`, `page_offset`, or `page_after_value`.

Let's consider a search for entries as an example.
```py
response = requests.post(
    f'{base_url}/entries/query',
    json={
        'query': {
            'results.material.elements': {
                'all': ['Ti', 'O']
            }
        },
        'pagination': {
            'page_size': 10
        }
    }
)
```

This will only result in a response with a maximum of 10 entries. The response will contain a
`pagination` object like this:
```json
{
    "page_size": 10,
    "order_by": "entry_id",
    "order": "asc",
    "total": 17957,
    "next_page_after_value": "--SZVYOxA2jTu_L-mSxefSQFmeyF"
}
```

In this case, the pagination is based on *after values*. This means that the search can
be continued with a follow up request at a certain point characterized by the
`next_page_after_value`. If you follow up with:

```py
response = requests.post(
    f'{base_url}/entries/query',
    json={
        'query': {
            'results.material.elements': {
                'all': ['Ti', 'O']
            }
        },
        'pagination': {
            'page_size': 10,
            'page_after_value': '--SZVYOxA2jTu_L-mSxefSQFmeyF'
        }
    }
)
```
You will get the next 10 results.

Here is a full example that collects the first 100 formulas from entries that match
a certain query by paginating.

```python
--8<-- "examples/docs/api/pagination.py"
```

### Authentication

Most of the API operations do not require any authorization and can be freely used
without a user or credentials. However, to upload, edit, or view your own and potentially unpublished data, the API needs to authenticate you.

The NOMAD API uses OAuth and tokens to authenticate users. We provide simple operations
that allow you to acquire an *access token* via username and password:

```py
import requests

response = requests.get(
    '{{ nomad_url() }}/v1/auth/token', params=dict(username='myname', password='mypassword'))
token = response.json()['access_token']

response = requests.get(
    '{{ nomad_url() }}/v1/uploads',
    headers={'Authorization': f'Bearer {token}'})
uploads = response.json()['data']
```

If you have the [NOMAD Python package](./pythonlib.md) installed. You can use its `Auth`
implementation:

```py
import requests
from nomad.client import Auth

response = requests.get(
    '{{ nomad_url() }}/v1/uploads',
    auth=Auth(user='myname or email', password='mypassword'))
uploads = response.json()['data']
```

To use authentication in the dashboard, simply use the Authorize button. The
dashboard GUI will manage the access token and use it while you try out the various
operations.

#### App token

If the short-term expiration of the default *access token* does not suit your needs,
you can request an *app token* with a user-defined expiration. For example, you can
send the GET request `/auth/app_token?expires_in=86400` together with some way of
authentication, e.g. header `Authorization: Bearer <access token>`. The API will return
an app token, which is valid for 24 hours in subsequent request headers with the format
`Authorization: Bearer <app token>`. The request will be declined if the expiration is
larger than the maximum expiration defined by the API config.

!!! warning
    Despite the name, the app token is used to impersonate the user who requested it.
    It does not discern between different uses and will only become invalid once it
    expires (or when the API's secret is changed).

## Search for entries

See [getting started](#getting-started) for a typical search example. Combine the [different
concepts](#common-concepts) above to create the queries that you need.

Searching for entries is typically just an initial step. Once you know what entries exist
you'll probably want to do one of the following things.

## Download raw files
You can use [queries](#queries) to download raw files, but typically you don't want to
download file-by-file or entry-by-entry. Therefore, we allow to download a large set of
files in one big zip-file. Here, you might want to use a program like *curl* to download
directly from the shell:

```
curl "{{ nomad_url() }}/v1/entries/raw?results.material.elements=Ti&results.material.elements=O" -o download.zip
```

## Access processed data (archives)
Above under [getting started](#getting started), you've already learned how to access
archive data. A special feature of the archive API functions is that you can define what is `required`
from the archives.

```py
response = requests.post(
    f'{base_url}/entries/archive/query',
    json={
        'query': ...,
        'pagination': ...,
        'required': {
            'workflow': {
                'calculation_result_ref': {
                    'energy': '*',
                    'system_ref': {
                        'chemical_composition': '*'
                    }
                }
            }
        }
    })
```

{{ doc_snippet('archive-required') }}

{{ metainfo_data() }}

## Limits

The API allows you to ask many requests in parallel and to put a lot of load
on NOMAD servers. Since this can accidentally or deliberately reduce the service
quality for other, we have to enforce a few limits.

- *rate limit*: you can only run a certain amount of requests at the same time
- *rate limit*: you can only run a certain amount of requests per second
- *api limit*: many API endpoints will enforce a maximum page size

If you get responses with an HTTP code **503 Service Unavailable**, you are hitting
a rate limit and you cannot use the service until you fall back into our limits. Consider,
to ask fewer requests in a larger time frame.

Rate limits are enforced based on your IP address. Please note that when you or your
colleagues are sharing a single external IPs from within a local network, e.g.
via [NAT](https://en.wikipedia.org/wiki/Network_address_translation),
you are also sharing the rate limits.
Depending on the NOMAD installation, these limits can be as low as 30 requests per second
or 10 concurrent requests.

Consider to use endpoints that allow you to retrieve full
pages of resources, instead of endpoints that force you to access resources one at a time.
See also the sections on [types of data](#different-kinds-of-data) and [pagination](#pagination).

However, pagination also has its limits and you might ask for pages that are too large.
If you get responses in the 400 range, e.g. **422 Unprocessable Content** or **400 Bad request**,
you might hit an api limit. Those responses are typically accompanied by an error message
in the response body that will inform you about the limit, e.g. the maximum allowed
page size.