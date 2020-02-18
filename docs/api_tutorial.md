# API Tutorial

This tutorial assumes that you want to

- upload some data
- publish the data
- find it
- download it again

## Prequisites

### Python
The tutorial was tested with Python 3, but it might as well work with Python 2.

### Python packages
We do not assume many specific python packages. Only the *bravado* package (available
via pipy) is required. It allows us to use the nomad ReST API in a more friendly and
pythonic way. You can simply install it the usual way

Optionally, if you need to access your private data, the package *python-keycloak* is
required to conveniently acquire the necessary tokens to authenticate your self towards
NOMAD.

```
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
nomad_url = 'https://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/latest/api'
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
            server_url='https://repository.nomad-coe.eu/fairdi/keycloak/auth/',
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

## Uploading data
Now, we can look at actually using the nomad API. The API is divided into several
modules: *uploads*, *repo*, *archive*, *raw*, etc. Each provided functionality for
a certain aspect of nomad.

The *uploads* endpoints can be used to, you guessed it, upload your data. But they
also allow to get process on the upload processing; inspect, delete, and publish uploads;
and get details about the uploaded data, which code input/output files where found, etc.

### Uploading a file

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

### Supervising the processing

Once uploaded, nomad will extract the file, identify code data, parse and normalize the
data. We call this *processing* and *processing* consists of *tasks* (uploading, extracting, parsing).
You can consistently pull the API, to get an update on the processing and check if all
tasks have completed.

```python
while upload.tasks_running:
    upload = client.uploads.get_upload(upload_id=upload.upload_id).response().result
    time.sleep(5)
    print('processed: %d, failures: %d' % (upload.processed_calcs, upload.failed_calcs))
```

Once there are no more tasks running, you can check if your upload was a success. If it
was not successful, you can also delete the upload again:
```python
if upload.tasks_status != 'SUCCESS':
    print('something went wrong')
    print('errors: %s' % str(upload.errors))

    # delete the unsuccessful upload
    client.uploads.delete_upload(upload_id=upload.upload_id).response().result
    sys.exit(1)
```

Of course, you can also visit the nomad GUI
([https://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/latest/gui/uploads](https://labdev-nomad.esc.rzg.mpg.de/fairdi/nomad/latest/gui/uploads))
to inspect your uploads. (You might click reload, if you had the page already open.)


### Publishing your upload
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


## Searching for data
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


## Downloading data
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

## Using *curl* to access the API

The shell tool *curl* can be used to call most API endpoints. Most endpoints for searching
or downloading data are only **GET** operations controlled by URL parameters. For example:

Downloading data:
```
curl http://repository.nomad-coe.eu/app/api/raw/query?upload_id=<your_upload_id> -o download.zip
```

It is a litle bit trickier, if you need to authenticate yourself, e.g. to download
not yet published or embargoed data. All endpoints support and most require the use of
an access token. To acquire an access token from our usermanagement system with curl:
```
curl --data 'grant_type=password&client_id=nomad_public&username=<your_username>&password=<your password>' \
    https://repository.nomad-coe.eu/fairdi/keycloak/auth/realms/fairdi_nomad_prod/protocol/openid-connect/token
```

You can use the access-token with:
```
curl -H 'Authorization: Bearer <you_access_token>' \
    http://repository.nomad-coe.eu/app/api/raw/query?upload_id=<your_upload_id> -o download.zip
```

## Conclusions
This was just a small glimpse into the nomad API. You should checkout our [swagger-ui](https://repository.nomad-coe.eu/app/api/) for more details on all the API endpoints and their parameters. You can explore the
API via the swagger-ui and even try it in your browser.
