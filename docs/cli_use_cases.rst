.. _cli_use_cases:

Use cases
*********

Mirroring data between production environments
""""""""""""""""""""""""""""""""""""""""""""""
Sometimes you would wish to transfer data between separate deployments of the
NOMAD infrastructure. This use case covers the situation when the deployments
are up and running and both have access to the underlying file storage, part of
which is mounted inside each container under :code:`.volumes/fs`.

With both the source and target deployment running, you can use the
:code::ref:`cli_ref:mirror` command to transfer the data from source to target. The
mirror will copy everything: i.e. the raw data, archive data and associated
metadata in the database.

The data to be mirrored is specified by using a query API path. For example to
mirror the upload from source deployment to target deployment, you would use
the following CLI command inside the target deployment:

.. code-block:: sh

    nomad client -n <api_url> -u <username> -w <password> mirror <query_json> --source-mapping <target_docker_path>:<shared_path>

Here is a breakdown of the different arguments:

  * :code:`-n <url>`: Url to the API endpoint in the source deployment. This API will
    be queried to fetch the data to be mirrored. E.g.
    http://repository.nomad-coe.eu/api
  * :code:`-u <username>`: Your username that is used for authentication in the API call.
  * :code:`-w <password>`: Your password that is used for authentication in the API call.
  * :code:`mirror <query>`: Your query as a JSON dictionary. See the documentation for
    available keywords. E.g. "{"upload_id: "<upload_id>"}"
  * :code:`--source-mapping <mapping>`: The deployments use a separate folder to store
    the archive and raw data. To correctly find the data that should be
    mirrored, the absolute path on the filesystem that is shared between the
    deployments needs to be provided. E.g. *.volumes/fs:/nomad/fairdi/prod/fs*.
    The first part of this mapping indicates a docker volume path
    (*.volumes/fs* in this example) that should be mapped to the second
    filepath on the shared filesystem (*/nomad/fairdi/prod/fs* in this example).

Updating the AFLOW prototype information
""""""""""""""""""""""""""""""""""""""""
NOMAD uses the `AFLOW prototype library
<http://www.aflowlib.org/CrystalDatabase/>`_ to link bulk crystal entries with
prototypical structures based on their symmetry. The
:ref:`cli_ref:prototypes-update` subcommand can be used to update this
database from the online information provided by AFLOW. The command produces a
prototype dataset as a python module.

The dataset should be recreated if the AFLOW dataset has been updated or if the
symmetry matching routine used within NOMAD is updated (e.g. the symmetry
tolerance is modified). To produce a new dataset run the following command:



.. code-block:: sh

    nomad admin ops prototypes-update <module_path>

Here is a breakdown of the different arguments:

  * :code:`<module_name>`: Name of the python module in which the data should
    be stored. If the file does not exist it will be created. The prototype
    data used by NOMAD is under the path:
    *nomad/normalizing/data/aflow_prototypes.py*

The command also provides a :code:`--matches-only` flag for only updating the
dataset entry that is used for matching the prototypes. This means that the
online information from AFLOW is not queried. This makes the process faster
e.g. in the case when you want only to update the matches after modifying the
symmetry routines.
