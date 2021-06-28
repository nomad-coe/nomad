# How to use the APIs

The NOMAD Repository and Archive offers all its functionality through an application
programming interface (API). More specifically a [RESTful HTTP API](https://en.wikipedia.org/wiki/Representational_state_transfer) that allows you
to use NOMAD as a set of resources (think data) that can be uploaded, accessed, downloaded,
searched for, etc. via [HTTP requests](https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol).

There are different tools and libraries to use the NOMAD API that come with different
trade-offs between expressiveness, learning curve, and convinience:

- use an HTTP program like *curl* or *wget* to directly use NOMAD from within a shell
- use a generic Python HTTP library like [requests](https://requests.readthedocs.io/en/master/)
- use more specific Python libraries like [bravado](https://github.com/Yelp/bravado) that turn HTTP requests into NOMAD
  specific function calls based on an [OpenAPI spec](https://swagger.io/specification/) that NOMAD offers and that describes our API
- directly in the browser via our generated [swagger dashboard](../api/)
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

## Using bravado and our OpenAPI spec

The Python library *bravado* is also an HTTP client, but instead of generic *GET URL*
style functions, it takes a formal specification of the NOMAD API and provides NOMAD
specific functions for you.

```python
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
```

Read on and learn how to install bravado and perform various tasks, like:

- upload some data
- publish the data
- find it
- download it again

### Python packages
We do not assume many specific python packages. Only the *bravado* package (available
via pipy) is required. It allows us to use the nomad ReST API in a more friendly and
pythonic way. You can simply install it the usual way

Optionally, if you need to access your private data, the package *python-keycloak* is
required to conveniently acquire the necessary tokens to authenticate your self towards
NOMAD.

```sh
pip install bravado
pip install python-keycloak
```

For the following code snippets, we need the following imports:

```python
from bravado.requests_client import RequestsClient
from bravado.client import SwaggerClient
from bravado.exception import HTTPNotFound
from urllib.parse import urlparse
import time
import os.path
import sys
```

And optionally:
```python
from bravado.requests_client import RequestsClient, Authenticator
from keycloak import KeycloakOpenID
```

### An example file
Lets assume you have an example upload file ready. Its a `.zip` (`.tgz` would also work)
with some *VASP* data from a single run at `/example/AcAg/vasprun.xml`, `/example/AcAg/OUTCAR`, ...
Lets keep the filename in a variable:

```python
upload_file = 'example.zip'
```

### Nomad
We need to know the nomad installation to use and its respective API URL. To upload
data you also need an account (email, password). The toy account used here, should be
available on most nomad installations:

```python
nomad_url = 'https://nomad-lab.eu/prod/rae/api'
user = 'leonard.hofstadter@nomad-fairdi.tests.de'
password = 'password'
```

### Using bravado
Bravado reads a ReST API's definition from a `swagger.json` as it is provided by
many APIs, including nomad's of course.

```python
host = urlparse(nomad_url).netloc.split(':')[0]
http_client = RequestsClient()
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)
```

Bravado also allows to use authentication, if required. The following would be a typical setup:

```python
class KeycloakAuthenticator(Authenticator):
    """ A bravado authenticator for NOMAD's keycloak-based user management. """
    def __init__(self, user, password):
        super().__init__(host=urlparse(nomad_url).netloc.split(':')[0])
        self.user = user
        self.password = password
        self.token = None
        self.__oidc = KeycloakOpenID(
            server_url='https://nomad-lab.eu/fairdi/keycloak/auth/',
            realm_name='fairdi_nomad_prod',
            client_id='nomad_public')

    def apply(self, request):
        if self.token is None:
            self.token = self.__oidc.token(username=self.user, password=self.password)
            self.token['time'] = time()
        elif self.token['expires_in'] < int(time()) - self.token['time'] + 10:
            try:
                self.token = self.__oidc.refresh_token(self.token['refresh_token'])
                self.token['time'] = time()
            except Exception:
                self.token = self.__oidc.token(username=self.user, password=self.password)
                self.token['time'] = time()

        request.headers.setdefault('Authorization', 'Bearer %s' % self.token['access_token'])

        return request

http_client = RequestsClient()
http_client.authenticator = KeycloakAuthenticator(user=user, password=password)
client = SwaggerClient.from_url('%s/swagger.json' % nomad_url, http_client=http_client)
```

### Uploading data
Now, we can look at actually using the nomad API. The API is divided into several
modules: *uploads*, *repo*, *archive*, *raw*, etc. Each provided functionality for
a certain aspect of nomad.

The *uploads* endpoints can be used to, you guessed it, upload your data. But they
also allow to get process on the upload processing; inspect, delete, and publish uploads;
and get details about the uploaded data, which code input/output files where found, etc.

#### Uploading a file

Its simple, since bravado supports uploading files:

```python
with open(upload_file, 'rb') as f:
    upload = client.uploads.upload(file=f).response().result
```

If you already have you file on the nomad servers, e.g. under `/nomad/my_files/example.zip`,
you can skip the actual upload and say:
```python
upload = client.uploads.upload(local_path='/nomad/my_files/example.zip').response().result
```

#### Supervising the processing

When files are added to an upload, nomad will initiate a *process* to extract/update the
files, identify code data, parse and normalize the data.

You can continuously pull the API, to get an update on the processing and check if the
processing has been completed.

```python
while upload.process_running:
    upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
    time.sleep(5)
    print('processed: %d, failures: %d' % (upload.processed_calcs, upload.failed_calcs))
```

Once the process completes, you can check if your upload was a success. If it
was not successful, you can also delete the upload again:
```python
if upload.process_status != 'SUCCESS':
    print('something went wrong')
    print('errors: %s' % str(upload.errors))

    # delete the unsuccessful upload
    client.uploads.delete_upload(upload_id=upload.upload_id).response().result
    sys.exit(1)
```

Of course, you can also visit the nomad GUI
([https://nomad-lab.eu/prod/rae/gui/uploads](https://nomad-lab.eu/prod/rae/gui/uploads))
to inspect your uploads. (You might click reload, if you had the page already open.)


#### Publishing your upload
The uploaded data is only visible to you. We call this *staging*. After the processing
was successful and you are satisfied with our processing, you have to publish the upload.
This also allows you to add additional meta-data to your upload (e.g. comments, references, coauthors, etc.).
Here you also determine, if you want an *embargo* on your data.

Once the data was published, you cannot delete it anymore. You can skip this step, but
the reset of the tutorial, will only work for you, because the data is only visible to you.

To initiate the publish and provide further data:
```python
client.uploads.exec_upload_operation(upload_id=upload.upload_id, payload={
    'operation': 'publish',
    'metadata': {
        'comment': 'Data from a cool external project',
        'references': ['http://external.project.eu'],
        # 'coauthors': ['sheldon.cooper@ucla.edu'],  this does not yet work with emails
        # 'external_id': 'external_id'  this does also not work, but we could implement something like this
    }
})
```

Publishing, also might take a while. You can inspect this analog to the upload processing:

```python
while upload.process_running:
    try:
        upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
        time.sleep(1)
    except HTTPNotFound:
        # upload gets deleted from the upload staging area once published
        break
```

This time we needed some exception handling, since the upload will be removed from the
staging area, and you will get a 404 on the `uploads` endpoint.


### Searching for data
The *repo* part of the API contains a *search* endpoint that support many different
quantities to search for. These include `formula` (e.g. *AcAg*), `system` (e.g. *bulk/2D/atom*), `spacegroup`, `authors`, `code` (e.g. *VASP*), etc.
In the following example, we search for the specific path segment `AcAg`.

```python
result = client.repo.search(paths='AcAg').response().result
if result.pagination.total == 0:
    print('not found')
elif result.pagination.total > 1:
    print('my ids are not specific enough, bummer ... or did I uploaded stuff multiple times?')
calc = result.results[0]
print(calc)
```

The result of a search always contains the key `pagination` with pagination data (`total`, `page`, `per_page`) and `results` with an array of the search result. The search results depend on
the type of search and their is no formal swagger model for it, therefore you get plain
dictionaries.


### Downloading data
The *raw* api allows to download data. You can do that either via bravado:
```python
client.raw.get(upload_id=calc['upload_id'], path=calc['mainfile']).response()
```

In case of published data, you can also create plain URLs and use a tool like *curl*:
```python
print('%s/raw/%s/%s' % (nomad_url, calc['upload_id'], calc['mainfile']))
print('%s/raw/%s/%s/*' % (nomad_url, calc['upload_id'], os.path.dirname(calc['mainfile'])))
```

There are different options to download individual files, or zips with multiple files.

### Using *curl* to access the API

The shell tool *curl* can be used to call most API endpoints. Most endpoints for searching
or downloading data are only **GET** operations controlled by URL parameters. For example:

Downloading data:
```sh
curl http://nomad-lab.eu/prod/rae/api/raw/query?upload_id=<your_upload_id> -o download.zip
```

It is a litle bit trickier, if you need to authenticate yourself, e.g. to download
not yet published or embargoed data. All endpoints support and most require the use of
an access token. To acquire an access token from our usermanagement system with curl:
```sh
curl --data 'grant_type=password&client_id=nomad_public&username=<your_username>&password=<your password>' \
    https://nomad-lab.eu/fairdi/keycloak/auth/realms/fairdi_nomad_prod/protocol/openid-connect/token
```

You can use the access-token with:
```sh
curl -H 'Authorization: Bearer <you_access_token>' \
    http://nomad-lab.eu/prod/rae/api/raw/query?upload_id=<your_upload_id> -o download.zip
```

### Conclusions
This was just a small glimpse into the nomad API. You should checkout our
[swagger-ui](nomad-lab.eu/prod/rae/api/)
for more details on all the API endpoints and their parameters. You can explore the
API via the swagger-ui and even try it in your browser.


## NOMAD's Python client library

This library is part devevloped by NOMAD. It is supposed to provide more powerful
access to common yet complex tasks. It currently only support access to the NOMAD
Archive. It has its separate documentation [here](archive.html).