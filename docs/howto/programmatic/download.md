A common use-case for the NOMAD API is to download large amounts of NOMAD data.
In this how-to guide, we use curl and API endpoints
that stream .zip files to download many resources with a single request directly from
the command line.

## Prerequisites

Here is some background information to understand the examples better.

### curl

To download resources from a REST API using curl, you can utilize the powerful command-line tool
to send HTTP requests and retrieve the desired data. Curl provides a simple and efficient
way to interact with RESTful APIs, allowing you to specify the necessary headers, parameters,
and authentication details. Whether you need to download files, retrieve JSON data, or access
other resources, curl offers a flexible and widely supported solution for programmatically
fetching data from REST APIs.

### Raw files vs processed data

We are covering two types of resources: *raw files* and *processed data*.
The former is organized into uploads and sub directory. The organization depends
on how the author was providing the files.
The later is organized by entries. Each NOMAD entry has corresponding structured data.

Endpoints that target raw files typically contain `raw`, e.g. `uploads/<id>/raw`
or `entries/raw/query`. Endpoints that target processed data contain `archive`
(because we call the entirety of all processed data the NOMAD Archive), e.g.
`entries/<id>/archive` or `entries/archive/query`.

### Entry vs upload

API endpoints for data download either target *entries* or *uploads*. For both types
of entities, endpoints for raw files and processed data (as well as searchable metadata)
exist. API endpoint paths start with the entity, e.g. `uploads/<id>/raw` or `entries/<id>/raw`.

## Download a whole upload

Let's assume you want to download an entire upload. In this example the upload id is
`wW45wJKiREOYTY0ARuknkA`.

```sh
curl -X GET "{{ nomad_url() }}/v1/uploads/wW45wJKiREOYTY0ARuknkA/raw" -o download.zip
```

This will create a `download.zip` file in the current folder. The zip file will contain
the raw file directory of the upload.

The used `uploads/<id>/raw` endpoint is only available for published uploads. For those,
all raw files have already been
packed into a zip file and this endpoint simply lets you download it. This is the simplest
and most reliable download implementation.

Alternatively, you can download specific files or sub-directories. This method is available
for all uploads. Including un-published uploads.

```sh
curl -X GET "{{ nomad_url() }}/v1/uploads/wW45wJKiREOYTY0ARuknkA/raw/?compress=true" -o download.zip
```

This endpoint looks very similar, but is implemented very differently. Note that we
put an empty path `/` to the end of the URL, plus a query parameter `compress=true`.
The path can be replaced with any directory or file path in the upload; `/` would denote the
whole upload. The query parameter says that we want to download the whole directory
as a zip file, instead of an individual file. This traverses through all files and
creates a zip file on the fly.

## Download a whole dataset

Now let's assume that you want to download all raw files that are associated with
all the entries of an entire dataset. In this example the dataset DOI is
`10.17172/NOMAD/2023.11.17-2`.

```sh
curl -X POST "{{ nomad_url() }}/v1/entries/raw/query" \
-H 'Content-Type: application/json' \
-d '{
    "query": {
        "datasets.doi": "10.17172/NOMAD/2023.11.17-2"
    }
}' \
-o download.zip
```

This time, we use the `entries/raw/query` endpoint that is based on entries and not on uploads. Here, we
select entries with a query. In the example, we query for the dataset DOI, but you
can replace this with any NOMAD search query (look out for the `<>` symbol on the
[search interface]({{ nomad_url() }}/../gui/search/entries)). The zip file will contain all raw files from all the
directories that have the mainfile of one of the entries that match the queries.

This might not necessarily download all uploaded files. Alternatively, you can use a query to
get all upload ids and then use the method from the previous section:

```sh
curl -X POST "{{ nomad_url() }}/v1/entries/query" \
-H 'Content-Type: application/json' \
-d '{
    "query": {
        "datasets.doi": "10.17172/NOMAD/2023.11.17-2"
    },
    "pagination": {
        "page_size": 0
    },
    "aggregations": {
        "upload_ids": {
            "terms": {
                "quantity": "upload_id"
            }
        }
    }
}'
```

The last command will print JSON data that contains all the upload ids. It uses
the `entries/query` endpoint that allows you to query NOMAD's search.
It does not return any results (`page_size: 0`),
but performs an aggregation over all search results and collects the upload ids
from all entries.

## Download some processed data for a whole dataset

Similar to raw files, you can also download processed data. This is also an
entry based operation based on a query. This time we also specify a `required`
to explain which parts of the processed data, we are interested in:

```sh
curl -X POST "{{ nomad_url() }}/v1/entries/archive/download/query" \
-H 'Content-Type: application/json' \
-d '{
    "query": {
        "datasets.doi": "10.17172/NOMAD/2023.11.17-2"
    },
    "required": {
        "metadata": {
            "entry_id": "*",
            "mainfile": "*",
            "upload_id": "*"
        },
        "results": {
            "material": "*"
        },
        "run": {
            "system[-1]": {
                "atoms": "*"
            }
        }
    }
}' \
-o download.zip
```

Here we use the `entries/archive/download/query` endpoint. The result is a zip file
with one json file per entry. There are no directories and the files are named
`<entry-id>.json`. To associate the json files with entries, you should require
information that tells you more about the entries, e.g. `required.metadata.mainfile`.

See also the [How to access processed data](./archive_query.md) how-to guide.

