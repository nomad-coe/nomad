# How to publish data using python

## Uploading, changing metadata, and publishing via python API

The [NOMAD API](https://nomad-lab.eu/prod/rae/docs/api.html){:target="_blank"} allows uploading, publishing, etc. using a local python environment, as an alternative to the NOMAD GUI. An overview of all API functionalities is provided in [How to use the API](api.md)

We have prepare some simple python functions to facilitate use of this API. For use as demonstrated below, copy the following code into a file called NOMAD_API.py:

```python
import requests

def get_authentication_token(nomad_url, username, password):
    '''Get the token for accessing your NOMAD unpublished uploads remotely'''
    try:
        response = requests.get(
            nomad_url + 'auth/token', params=dict(username=username, password=password), timeout=10)
        token = response.json().get('access_token')
        if token:
            return token

        print('response is missing token: ')
        print(response.json())
        return
    except Exception:
        print('something went wrong trying to get authentication token')
        return


def create_dataset(nomad_url, token, dataset_name):
    '''Create a dataset to group a series of NOMAD entries'''
    try:
        response = requests.post(
            nomad_url + 'datasets/',
            headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
            json={"dataset_name": dataset_name},
            timeout=10
            )
        dataset_id = response.json().get('dataset_id')
        if dataset_id:
            return dataset_id

        print('response is missing dataset_id: ')
        print(response.json())
        return
    except Exception:
        print('something went wrong trying to create a dataset')
        return

def upload_to_NOMAD(nomad_url, token, upload_file):
    '''Upload a single file for NOMAD upload, e.g., zip format'''
    with open(upload_file, 'rb') as f:
        try:
            response = requests.post(
                nomad_url + 'uploads',
                headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
                data=f, timeout=30)
            upload_id = response.json().get('upload_id')
            if upload_id:
                return upload_id

            print('response is missing upload_id: ')
            print(response.json())
            return
        except Exception:
            print('something went wrong uploading to NOMAD')
            return

def check_upload_status(nomad_url, token, upload_id):
    '''
    # upload success => returns 'Process publish_upload completed successfully'
    # publish success => 'Process publish_upload completed successfully'
    '''
    try:
        response = requests.get(
            nomad_url + 'uploads/' + upload_id,
            headers={'Authorization': f'Bearer {token}'}, timeout=30)
        status_message = response.json().get('data').get('last_status_message')
        if status_message:
            return status_message

        print('response is missing status_message: ')
        print(response.json())
        return
    except Exception:
        print('something went wrong trying to check the status of upload' + upload_id)
        # upload gets deleted from the upload staging area once published...or in this case something went wrong
        return

def edit_upload_metadata(nomad_url, token, upload_id, metadata):
    '''
    Example of new metadata:
    upload_name = 'Test_Upload_Name'
    metadata = {
        "metadata": {
        "upload_name": upload_name,
        "references": ["https://doi.org/xx.xxxx/xxxxxx"],
        "datasets": dataset_id,
        "embargo_length": 0,
        "coauthors": ["coauthor@affiliation.de"],
        "comment": 'This is a test upload...'
        },
    }
    '''

    try:
        response = requests.post(
            nomad_url+'uploads/' + upload_id + '/edit',
            headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
            json=metadata, timeout=30)
        return response
    except Exception:
        print('something went wrong trying to add metadata to upload' + upload_id)
        return

def publish_upload(nomad_url, token, upload_id):
    '''Publish an upload'''
    try:
        response = requests.post(
            nomad_url+'uploads/' + upload_id + '/action/publish',
            headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
            timeout=30)
        return response
    except Exception:
        print('something went wrong trying to publish upload: ' + upload_id)
        return
```

Now, we will demonstrate how to use these functions. Within a notebook or python script, import the above functions:

```python
from Nomad_API import *`
```

Define the following user information:
```python
username = 'nomad_email@affiliation.edu'
password = 'password'
```

Define the NOMAD API endpoint:
```python
# nomad_url = 'https://nomad-lab.eu/prod/v1/api/v1/'  # production nomad
nomad_url = 'https://nomad-lab.eu/prod/v1/test/api/v1/'  # test nomad (deleted occassionally)
```

Get a token for accessing your unpublished uploads:

```python
token = get_authentication_token(nomad_url, username, password)
```

Create a dataset for grouping uploads that belong to, e.g., a publication:

```python
dataset_id = create_dataset(nomad_url, token, 'Test_Dataset')
```

Upload some test data to NOMAD:

```python
upload_id = upload_to_NOMAD(nomad_url, token, 'test_data.zip')
```

Check the status to make sure the upload was processed correctly:

```python
last_status_message = check_upload_status(nomad_url, token, upload_id)
print(last_status_message)
```

The immediate result may be:

    'Waiting for results (level 0)'

After some time you will get:

    'Process process_upload completed successfully'

??? tip

    Some data, e.g., large systems or molecular dynamics trajectories, take some time to process. In this case, you can call the above function intermittantly, e.g., in a while loop with a sleep call in between, waiting for `last_status_message` to be "Process process_upload completed successfully"


Now that the upload processing is complete, we can add coauthors, references, and other comments, as well as link to a dataset and provide a proper name for the upload:

```python
metadata = {
    "metadata": {
    "upload_name": 'Test_Upload',
    "references": ["https://doi.org/xx.xxxx/x.xxxx"],
    "datasets": dataset_id,
    "embargo_length": 0,
    "coauthors": ["coauthor@affiliation.de"],
    "comment": 'This is a test upload...',
},
}
response = edit_upload_metadata(nomad_url, token, upload_id, metadata)
```

Check the upload again to make sure that the metadata was changed:

```python
last_status_message = check_upload_status(nomad_url, token, upload_id)
print(last_status_message)
```

    'Process edit_upload_metadata completed successfully'


Now, we are ready to publish:

```python
response = publish_upload(nomad_url, token, upload_id)
```

Once again check the status:

```python
last_status_message = check_upload_status(nomad_url, token, upload_id)
print(last_status_message)
```

    'Process publish_upload completed successfully'

