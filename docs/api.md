# How to use the APIs

The NOMAD Repository and Archive offers all its functionality through an application
programming interface (API). More specifically a [RESTful HTTP API](https://en.wikipedia.org/wiki/Representational_state_transfer) that allows you
to use NOMAD as a set of resources (think data) that can be uploaded, accessed, downloaded,
searched for, etc. via [HTTP requests](https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol).

There are different tools and libraries to use the NOMAD API that come with different
trade-offs between expressiveness, learning curve, and convinience:

- use an HTTP program like *curl* or *wget* to directly use NOMAD from within a shell
- use a generic Python HTTP library like [requests](https://requests.readthedocs.io/en/master/)
- directly in the browser via our generated [OpenAPI dashboard](../api/v1)
- use the NOMAD Python client library, which offers custom and more powerful
  implementations for certain tasks (currently only for accessing the NOMAD Archive)

This set of tutorials provides a few examples for common NOMAD tasks using the various
options.

## Using *curl* (or *wget*)

Terminal programs like *curl* act as an HTTP client and allow you to send requests and
display or store the respective responses. HTTP basically allows you to GET, POST, PUT,
and DELETE "resources" on a remote server. These resources are identified via URLs (=uniform
resource locator). URLs usually consists of a protocol (e.g. HTTP), a domain (our servers),
a path (a place on our servers), and query parameters (additional options).

NOMAD provides three main set of resources: **repo** (i.e. the NOMAD Repository), **raw**
(raw uploaded files), **archive** (i.e. the NOMAD Archive). Within all these resource sets
you have endpoints that either allow you to directly locate a NOMAD entry (i.e. an
uploaded code run) or to ask a query to locate many NOMAD entries at the same time. Here,
the **repo** will return the repository metadata for said entries, **archive** the archive
data, ...

Let's say you want to see the repository metadata (i.e. the information that you see in
our gui) for entries that fit search criteria, like compounds having atoms *Si* and *O* in
it:

```sh
curl -X GET "http://nomad-lab.eu/prod/rae/api/repo/?atoms=Si&atoms=O"
```

Here we used curl to send an HTTP GET request to return the resource located by the given URL.
In practice you can omit the `-X GET` (which is the default) and you might want to format
the output:

```sh
curl "http://nomad-lab.eu/prod/rae/api/repo/?atoms=Si&atoms=O" | python -m json.tool
```

You'll see the the metadata of the first 10 entries that match your criteria. There
are various other query parameters. You find a full list in the generated [swagger dashboard
of our API](https://nomad-lab.eu/prod/rae/api/).

Besides search criteria you can determine how many results (`per_page`) and what page of
results should be returned (`page`). If you want to go beyond the first 10.000 results
you can use our *scroll* API (`scroll=true`, `scroll_after`). You can limit what properties
should be returned (`include`, `exclude`). See the the generated [swagger dashboard
of our API](https://nomad-lab.eu/prod/rae/api/) for more parameters.

If you use the [NOMAD Repository and Archive search interface](https://nomad-lab.eu/prod/rae/gui/search)
and create a query, you can click th a **<>**-button (right and on top of the result list).
This will give you some code examples with URLs for your search query.

Similar functionality is offered to download archive or raw data. Let's say you have
identified an entry (given via a `upload_id`/`calc_id`, see the query output), and
you want to download it:

```sh
curl "http://nomad-lab.eu/prod/rae/api/raw/calc/JvdvikbhQp673R4ucwQgiA/k-ckeQ73sflE6GDA80L132VCWp1z/*" -o download.zip
```

With `*` you basically requests all the files under an entry or path..
If you need a specific file (that you already know) of that calculation:

```sh
curl "http://nomad-lab.eu/prod/rae/api/raw/calc/JvdvikbhQp673R4ucwQgiA/k-ckeQ73sflE6GDA80L132VCWp1z/INFO.OUT"
```

You can also download a specific file from the upload (given a `upload_id`), if you know
the path of that file:

```sh
curl "http://nomad-lab.eu/prod/rae/api/raw/JvdvikbhQp673R4ucwQgiA/exciting_basis_set_error_study/monomers_expanded_k8_rgkmax_080_PBE/72_Hf/INFO.OUT"
```

If you have a query
that is more selective, you can also download all results. Here all compounds that only
consist of Si, O, bulk material simulations of cubic systems (currently ~100 entries):

```sh
curl "http://nomad-lab.eu/prod/rae/api/raw/query?only_atoms=Si&only_atoms=O&system=bulk&crystal_system=cubic" -o download.zip
```

Here are a few more examples for downloading the raw data of based on DOI or dataset.
You will have to encode non URL safe characters in potential dataset names (e.g. with a service like [www.urlencoder.org](https://www.urlencoder.org/)):

```sh
curl "http://nomad-lab.eu/prod/rae/api/raw/query?datasets.doi=10.17172/NOMAD/2020.03.18-1" -o download.zip
curl "http://nomad-lab.eu/prod/rae/api/raw/query?dataset=Full%20ahnarmonic%20stAViC%20approach%3A%20Silicon%20and%20SrTiO3" -o download.zip
```

In a similar way you can see the archive of an entry:

```sh
curl "http://nomad-lab.eu/prod/rae/api/archive/f0KQE2aiSz2KRE47QtoZtw/6xe9fZ9xoxBYZOq5lTt8JMgPa3gX" | python -m json.tool
```

Or query and display the first page of 10 archives:

```sh
curl "http://nomad-lab.eu/prod/rae/api/archive/query?only_atoms=Si&only_atoms=O" | python -m json.tool
```

## Using Python's *request* library

Similar to *curl* in the shell, you can use *requests* in Python. Its a generic HTTP
client library that allows you to send requests:

```python
import requests
import json

response = requests.get("http://nomad-lab.eu/prod/rae/api/archive/query?only_atoms=Si&only_atoms=O")
data = response.json()
print(json.dumps(data), indent=2)
```

## NOMAD's Python client library

This library is part devevloped by NOMAD. It is supposed to provide more powerful
access to common yet complex tasks. It currently only support access to the NOMAD
Archive. It has its separate documentation [here](pythonlib).