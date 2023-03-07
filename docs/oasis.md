# Operating an OASIS

Originally, NOMAD Central Repository is a service run at Max-Planck's computing facility in Garching, Germany.
However, the NOMAD software is Open-Source, and everybody can run it. Any service that
uses NOMAD software independently is called a *NOMAD OASIS*. A *NOMAD OASIS* does not
need to be fully isolated. For example, you can publish uploads from your OASIS to the
central NOMAD installation.

## Quick-start

- Find a linux computer.
- Make sure you have [docker](https://docs.docker.com/engine/install/) installed.
Docker nowadays comes with `docker compose` build in. Prior, you needed to
install the stand alone [docker-compose](https://docs.docker.com/compose/install/).
- Download our basic configuration files [nomad-oasis.zip](assets/nomad-oasis.zip)
- Run the following commands (skip `chown` on MacOS and Windows computers)

```sh
unzip nomad-oasis.zip
cd nomad-oasis
sudo chown -R 1000 .volumes
docker compose pull
docker compose up -d
curl localhost/nomad-oasis/alive
```

- Open [http://localhost/nomad-oasis](http://localhost/nomad-oasis) in your browser.

To run NORTH (the NOMAD Remote Tools Hub), the `hub` container needs to run docker and
the container has to be run under the docker group. You need to replace the default group
id `991` in the `docker-compose.yaml`'s `hub` section with your systems docker group id.
Run `id` if you are a docker user, or `getent group | grep docker` to find our your
systems docker gid. The user id 1000 is used as the nomad user inside all containers.

This is good as a quick test. We strongly recommend to read the following instructions
carefully and adapt the configuration files accordingly. The following might also
include meaningful help, if you run into problems.

## Before you start

### Hardware considerations

Of course this depends on how much data you need to manage and process. Data storage is
the obvious aspect here. NOMAD keeps all files that it manages as they are. The files
that NOMAD processes in addition (e.g. through parsing) are typically smaller than
the original raw files. Therefore, you can base your storage requirements based on the
size of the data files that you expect to manage. The additional mongo database and
elasticsearch index is comparatively small.

Storage speed is another consideration. You can work with NAS systems. All that NOMAD
needs is a "regular" POSIX filesystem as an interface. So everything you can (e.g. docker host)
mount should be fine. For processing data obviously relies on read/write speed, but
this is just a matter of convenience. The processing is designed to run as managed asynchronous
tasks. Local storage might be favorable for mongodb and elasticsearch operation, but it
is not a must.

The amount of compute resource (e.g. processor cores) is also a matter of convenience
(and amount of expected users). Four cpu-cores are typically enough to support a
research group and run application, processing, and databases in parallel. Smaller systems
still work, e.g. for testing.

There should be enough RAM to run databases, application, and processing at the same
time. The minimum requirements here can be quite low, but for processing the metadata
for individual files is kept in memory. For large DFT geometry-optimizations this can
add up quickly, especially if many CPU cores are available for processing entries in
parallel. We recommend at least 2GB per core and a minimum of 8GB. You also need to consider
RAM and CPU for running tools like jupyter, if you opt to use NOMAD NORTH.

### Using the central user management

Our recommendation is to use the central user management provided by nomad-lab.eu. We
simplified its use and you can use it out-of-the-box. You can even run your system
from `localhost` (e.g. for initial testing). The central user management system is not
communicating with your OASIS directly. Therefore, you can run your OASIS without
exposing it to the public internet.

There are two requirements. First, your users must be able to reach the OASIS. If a user is
logging in, she/he is redirected to the central user management server and after login,
she/he is redirected back to the OASIS. These redirects are executed by your user's browser
and do not require direct communication.

Second, your OASIS must be able to request (via HTTP) the central user management and central NOMAD
installation. This is necessary for non JWT-based authentication methods and to
retrieve existing users for data-sharing features.

The central user management will make future synchronizing data between NOMAD installations easier
and generally recommend to use the central system.
But in principle, you can also run your own user management. See the section on
[your own user management](#provide-and-connect-your-own-user-management).

## Docker and docker compose

We recommend the installation via docker and docker-compose. It is the most documented, simplest,
easiest to update, and generally the most frequently chosen option.

### Pre-requisites

NOMAD software is distributed as a set of docker containers and there are also other services required that can be run with docker.
Further, we use docker-compose to setup all necessary containers in the simplest way possible.

You will need a single computer, with **docker** and **docker-compose** installed. Refer
to the official [docker](https://docs.docker.com/engine/install/) (and [docker-compose](https://docs.docker.com/compose/install/))
documentation for installation instructions. Newer version of docker have a re-implementation
of docker-compose integrated as the `docker compose` sub-command. This should be fully
compatible and you might chose to can replace `docker compose` with `docker-compose` in this tutorial.

The following will run all necessary services with docker. These comprise: a **mongo**
database, an **elasticsearch**, a **rabbitmq** distributed task queue, the NOMAD **app**,
NOMAD **worker**, and NOMAD **gui**. In this [introduction](index.md#architecture),
you will learn what each service does and why it is necessary.

### Configuration

All docker containers are configured via docker-compose and the respective `docker-compose.yaml` file.
Further, we will need to mount some configuration files to configure the NOMAD services within their respective containers.

There are three files to configure:

- `docker-compose.yaml`
- `configs/nomad.yaml`
- `configs/nginx.conf`

In this example, we have all files in the same directory (the directory we are also working in).
You can download minimal example files [here](assets/nomad-oasis.zip).

#### docker-compose.yaml

The most basic `docker-compose.yaml` to run an OASIS looks like this:

```yaml
--8<-- "ops/docker-compose/nomad-oasis/docker-compose.yaml"
```

Changes necessary:

- The group in the value of the hub's user parameter needs to match the docker group
on the host. This should ensure that the user which runs the hub, has the rights to access the host's docker.

- On Windows or MacOS computers you have to run the `app` and `worker` container without `user: '1000:1000'` and the `north` container with `user: root`.

A few things to notice:

- The app, worker, and north service use the NOMAD docker image. Here we use the `latest` tag, which
gives you the latest *beta* version of NOMAD. You might want to change this to `stable`,
a version tag (format is `vX.X.X`, you find all releases [here](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/-/tags)), or a specific [branch tag](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/-/branches).
- All services use docker volumes for storage. This could be changed to host mounts.
- It mounts two configuration files that need to be provided (see below): `nomad.yaml`, `nginx.conf`.
- The only exposed port is `80` (proxy service). This could be changed to a desired port if necessary.
- The NOMAD images are pulled from our gitlab at MPCDF, the other services use images from a public registry (*dockerhub*).
- All containers will be named `nomad_oasis_*`. These names can be used later to reference the container with the `docker` cmd.
- The services are setup to restart `always`, you might want to change this to `no` while debugging errors to prevent indefinite restarts.
- Make sure that the `PWD` environment variable is set. NORTH needs to create bind mounts that require absolute paths and we need to pass the current working directory to the configuration from the PWD variable (see hub service in the `docker-compose.yaml`).
- The `hub` service needs to run docker containers. We have to use the systems docker group as a group. You might need to replace `991` with your
systems docker group id.


#### nomad.yaml

NOMAD app and worker read a `nomad.yaml` for configuration.

```yaml
--8<-- "ops/docker-compose/nomad-oasis/configs/nomad.yaml"
```

You should change the following:

- Replace `localhost` with the hostname of your server. I user-management will redirect your
users back to this host. Make sure this is the hostname, your users can use.
- Replace `deployment`, `deployment_url`, and `maintainer_email` with representative values.
The `deployment_url` should be the url to the deployment's api (should end with `/api`).
- You can change `api_base_path` to run NOMAD under a different path prefix.
- You should generate your own `north.jupyterhub_crypt_key`. You can generate one
with `openssl rand -hex 32`.
- On Windows or MacOS, you have to add `hub_connect_ip: 'host.docker.internal'` to the `north` section.

A few things to notice:

- Under `mongo` and `elastic` you can configure database and index names. This might
be useful, if you need to run multiple NOMADs with the same databases.
- All managed files are stored under `.volumes` of the current directory.

#### nginx.conf

The GUI container serves as a proxy that forwards requests to the app container. The
proxy is an nginx server and needs a configuration similar to this:

```none
--8<-- "ops/docker-compose/nomad-oasis/configs/nginx.conf"
```

A few things to notice:

- It configures the base path (`nomad-oasis`). It needs to be changed, if you use a different base path.
- You can use the server for additional content if you like.
- `client_max_body_size` sets a limit to the possible upload size.

You can add an additional reverse proxy in front or modify the nginx in the docker-compose.yaml
to [support https](http://nginx.org/en/docs/http/configuring_https_servers.html).
If you operate the GUI container behind another proxy, keep in mind that your proxy should
not buffer requests/responses to allow streaming of large requests/responses for `api/v1/uploads` and `api/v1/.*/download`.
An nginx reverse proxy location on an additional reverse proxy, could have these directives
to ensure the correct http headers and allows the download and upload of large files:
```nginx
client_max_body_size 35g;
proxy_set_header Host $host;
proxy_pass_request_headers on;
proxy_buffering off;
proxy_request_buffering off;
proxy_set_header X-Real-IP $remote_addr;
proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
proxy_pass http://<your-oasis-host>/nomad-oasis;
```

### Running NOMAD

If you prepared the above files, simply use the usual `docker compose` commands to start everything.

To make sure you have the latest docker images for everything, run this first:
```sh
docker compose pull
```

In the beginning and to simplify debugging, it is recommended to start the services separately:
```sh
docker compose up -d mongo elastic rabbitmq
docker compose up app worker gui
```

The `-d` option runs container in the background as *daemons*. Later you can run all at once:
```sh
docker compose up -d
```

Running all services also contains NORTH. When you use a tool in NORTH for the first time,
your docker needs to pull the image that contains this tool. Be aware that this might take longer
than timeouts allow and starting a tool for the very first time might fail.

You can also use docker to stop and remove faulty containers that run as *daemons*:
```sh
docker stop nomad_oasis_app
docker rm nomad_oasis_app
```

You can wait for the start-up with curl using the apps `alive` "endpoint":
```sh
curl http://<your host>/nomad-oasis/alive
```

If everything works, the gui should be available under:
```none
http://<your host>/nomad-oasis/gui/
```

If you run into problems, use the dev-tools of your browser to check the javascript logs
or monitor the network traffic for HTTP 500/400/404/401 responses.

To see if at least the api works, check
```none
http://<your host>/nomad-oasis/alive
http://<your host>/nomad-oasis/api/info
```

To see logs or 'go into' a running container, you can access the individual containers
with their names and the usual docker commands:

```sh
docker logs nomad_oasis_app
```

```sh
docker exec -ti nomad_oasis_app /bin/bash
```

If you want to report problems with your OASIS. Please provide the logs for

- nomad_oasis_app
- nomad_oasis_worker
- nomad_oasis_gui

### Provide and connect your own user management

NOMAD uses [keycloak](https://www.keycloak.org/) for its user management. NOMAD uses
keycloak in two ways. First, the user authentication uses the OpenID Connect/OAuth interfaces provided by keycloak.
Second, NOMD uses the keycloak realm-management API to get a list of existing users.
Keycloak is highly customizable and numerous options to connect keycloak to existing
identity providers exist.

This tutorial assumes that you have some understanding of what keycloak is and
how it works.

Start with the regular docker-compose installation above. Now you need to modify the
`docker-compose.yaml` to add a keycloak service. You need to modify the `nginx.conf` to add
another location for keycloak. You need to modify the `nomad.yaml` to tell nomad to
use your and not the official NOMAD keycloak.

```yaml
--8<-- "ops/docker-compose/nomad-oasis-with-keycloak/docker-compose.yaml"
```

A few notes:

- You have to change the `KEYCLOAK_FRONTEND_URL` variable to match your host and set a path prefix.
- The environment variables on the keycloak service allow to use keycloak behind the nginx proxy with a path prefix, e.g. `keycloak`.
- By default, keycloak will use a simple H2 file database stored in the given volume. Keycloak offers many other options to connect SQL databases.
- We will use keycloak with our nginx proxy here, but you can also host-bind the port `8080` to access keycloak directly.

Second, we add a keycloak location to the nginx config:
```nginx
--8<-- "ops/docker-compose/nomad-oasis-with-keycloak/configs/nginx.conf"
```

A few notes:

- Again, we are using `keycloak` as a path prefix. We configure the headers to allow
keycloak to pick up the rewritten url.

Third, we modify the keycloak configuration in the `nomad.yaml`:
```yaml
--8<-- "ops/docker-compose/nomad-oasis-with-keycloak/configs/nomad.yaml"
```

You should change the following:

- There are two urls to configure for keycloak. The `server_url` is used by the nomad
services to directly communicate with keycloak within the docker network. The `public_server_url`
is used by the UI to perform the authentication flow. You need to replace `localhost`
in `public_server_url` with `<yourhost>`.

A few notes:

- The particular `admin_user_id` is the Oasis admin user in the provided example realm
configuration. See below.

If you open `http://<yourhost>/keycloak/auth` in a browser, you can access the admin
console. The default user and password are `admin` and `password`.

Keycloak uses `realms` to manage users and clients. A default NOMAD compatible realm
is imported by default. The realm comes with a test user and password `test` and `password`.

A few notes on the realm configuration:

- Realm and client settings are almost all default keycloak settings.
- You should change the password of the admin user in the nomad realm.
- The admin user in the nomad realm has the additional `view-users` client role for `realm-management`
assigned. This is important, because NOMAD will use this user to retrieve the list of possible
users for managing co-authors and reviewers on NOMAD uploads.
- The realm has one client `nomad_public`. This has a basic configuration. You might
want to adapt this to your own policies. In particular you can alter the valid redirect URIs to
your own host.
- We disabled the https requirement on the default realm for simplicity. You should change
this for a production system.

### Backups

To backup your Oasis at least the file data and mongodb data needs to be saved. You determined the path to your file data (your uploads) during the installation. By
default all data is stored in a directory called `.volumes` that is created in the
current working directory of your installation/docker-compose. This directory can be backed up like any other file backup (e.g. rsync).

To backup the mongodb, please refer to the official [mongodb documentation](https://docs.mongodb.com/manual/core/backups/). We suggest a simple mongodump export that is backed up alongside your files. The default configuration mounts `.volumes/mongo` into the mongodb container (as `/backup`) for this purpose. You can use this to export the NOMAD mongo database. Combine this with rsync on the `.volumes` directory and everything should be set. To create a new mongodump run:

```sh
docker exec nomad_oasis_mongo mongodump -d nomad_oasis_v1 -o /backup
```

The elasticsearch contents can be reproduced with the information in the files and the mongodb.

To create a new oasis with the backup data, create the `.volumes` directory from your backup. Start the new oasis. Use mongorestore:

```sh
docker exec nomad_oasis_mongo mongorestore /backup
```

Now you still have to recreate the elasticsearch index:
```sh
docker exec nomad_oasis_app python -m nomad.cli admin uploads index
```

### A few useful NOMAD commands

The NOMAD command line interface (CLI) provides a few useful administrative functions. To use the NOMAD CLI, open a shell into the app container and run it from there:

```sh
docker exec -ti nomad_oasis_app bash
```

For example you can ls or remove uploads:
```sh
nomad admin uploads ls
nomad admin uploads rm -- <upload_id> <upload_id>
```

You can also reset the processing (of "stuck") uploads and reprocess:
```sh
nomad admin uploads reset -- <upload_id>
nomad admin uploads process -- <upload_id>
```

You can also use the CLI to wipe the whole installation:
```sh
nomad admin reset --i-am-really-sure
```

## Base Linux (without docker)

### Pre-requisites

We will run NOMAD from the *nomad-lab* Python package. This package contains all the necessary
code to run the processing, api, and gui. But, it is missing the necessary databases. You might be
able to run NOMAD in user space.

You will need:

- preferably a linux computer, which Python 3.9, preferable a Python virtual environment
- elasticsearch 7.x, running without users and authentication, preferable on the default settings
- mongodb 5.x, running without users and authentication, preferable on the default settings
- rabbitmq 3.x, running without users and authentication, preferable on the default settings
- nginx
- an empty directory to work in

### Install the NOMAD Python package

You should install everything in a virtual environment. For example like this:
```sh
virtualenv -p `which python3` nomadpyenv
source nomadpyenv/bin/activate
```

You can simply install the Python package from pypi:
```sh
pip install nomad-lab[all]
```

If you need the latest version, you can also download the latest package from our
"beta" installation.
```sh
curl "https://nomad-lab.eu/prod/v1/staging/dist/nomad-lab.tar.gz" -o nomad-lab.tar.gz
pip install nomad-lab.tar.gz[all]
```

### nomad.yaml

The `nomad.yaml` is our central config file. You should write a `nomad.yaml` like this:

```yaml
client:
  url: 'http://<your-host>/nomad-oasis/api'

services:
  api_base_path: '/nomad-oasis'
  admin_user_id: '<your admin user id>'

keycloak:
  realm_name: fairdi_nomad_prod
  username: '<your admin username>'
  password: '<your admin user password>'
  oasis: true

mongo:
    db_name: nomad_v0_8

elastic:
    index_name: nomad_v0_8
```

You need to change the following:

- Replace `your-host` and admin credentials respectively.
- `api_base_path` defines the path under which the app is run. It needs to be changed, if you use a different base path.

A few things to notice:

- Be secretive about your admin credentials; make sure this file is not publicly readable.

### nginx

You can generate a suitable `nginx.conf` with the `nomad` command line command that
comes with the *nomad-lab* Python package. If you server other content but NOMAD with
your nginx, you need to incorporate the config accordingly.

If you have a standard installation of nginx, this might work. Adapt to your set-up:
```sh
nomad admin ops nginx-conf > /etc/nginx/conf.d/default.conf
nginx -t
nginx -s reload
```

If you want to run nginx in docker, this might work. Adapt to your set-up:
```sh
nomad admin ops nginx-conf --host host.docker.internal > nginx.conf
docker run --rm -v `pwd`/nginx.conf:/etc/nginx/conf.d/default.conf -p 80:80 nginx:stable nginx -g 'daemon off;' &
```

### Running NOMAD

Before you start, we need to transfer your `nomad.yaml` config values to the GUI's
javascript. You need to repeat this, if you change your `nomad.yaml`. You can do this by running:
```
nomad admin ops gui-config
```

To run NOMAD, you must run two services. One is the NOMAD app, it serves the API and GUI:
```sh
--8<-- "run.sh"
```

the second is the NOMAD worker, that runs the NOMAD processing.
```
celery -A nomad.processing worker -l info -Q celery
```

This should give you a working OASIS at `http://<your-host>/<your-path-prefix>`.

## Kubernetes

*This is not yet documented.*

## Migrating from an older version (0.8.x to 1.x)

Between versions 0.10.x and 1.x we needed to change how archive and metadata data is stored
internally in files and databases. This means you cannot simply start a new version of
NOMAD on top of the old data. But there is a strategy to adapt the data. This should
work for data based on NOMAD >0.8.0 and <= 0.10.x.

The overall strategy is to create a new mongo database, copy all information, and then
reprocess all data for the new version.

First, shutdown the OASIS and remove all old containers.
```sh
docker compose stop
docker compose rm -f
```

Update your config files (`docker-compose.yaml`, `nomad.yaml`, `nginx.conf`) according
to the latest documentation (see above). Make sure to use index and database names that
are different. The default values contain a version number in those names, if you don't
overwrite those defaults, you should be safe.

Make sure you get the latest images and start the OASIS with the new version of NOMAD:
```sh
docker compose pull
docker compose up -d
```

If you go to the GUI of your OASIS, it should now show the new version and appear empty,
because we are using a different database and search index now.

To migrate the data, we created a command that you can run within your OASIS' NOMAD
application container. This command takes the old database name as an argument, it will copy
all data from the old mongodb to the current one. The default v8.x database name
was `nomad_fairdi`, but you might have changed this to `nomad_v0_8` as recommended by
our old Oasis documentation.

```sh
docker exec -ti nomad_oasis_app bash -c 'nomad admin upgrade migrate-mongo --src-db-name nomad_v0_8'
docker exec -ti nomad_oasis_app bash -c 'nomad admin uploads reprocess'
```

Now all your data should appear in your OASIS again. If you like, you can remove the
old index and database:

```sh
docker exec nomad_oasis_elastic bash -c 'curl -X DELETE http://elastic:9200/nomad_fairdi'
docker exec nomad_oasis_mongo bash -c 'mongo nomad_fairdi --eval "printjson(db.dropDatabase())"'
```

## Restricting access to your Oasis

An Oasis works exactly the same way the official NOMAD works. It is open and everybody
can access published data. Everybody with an account can upload data. This might not be
what you want.

Currently there are two ways to restrict access to your Oasis. First, you do not
expose the Oasis to the public internet, e.g. you only make it available on an intra-net or
through a VPN.

Second, we offer a simple white-list mechanism. As the Oasis administrator you provide a
list of accounts as part of your Oasis configuration. To use the Oasis, all users have to
be logged in and be on your white list of allowed users. To enable white-listing, you
can provide a list of NOMAD account email addresses in your `nomad.yaml` like this:

```
oasis:
    allowed_users:
        - user1@gmail.com
        - user2@gmail.com
```

## Performance considerations

If you run the OASIS on a single computer, like described here (either with docker or bare
linux), you might run into problems with processing large uploads. If the NOMAD worker
and app are run on the same computer, the app might become unresponsive, when the worker
consumes all system resources.

By default, the worker container might have as many worker processes as the system as CPU cores.
In addition, each worker process might spawn additional threads and consume
more than one CPU core.

There are multiple ways to restrict the worker's resource consumption:

- limit the number of worker processes and thereby lower the number of used cores
- disable or restrict multithreading
- limit available CPU utilization of the worker's docker container with docker

#### Limit the number of worker processes

The worker uses the Python package celery. Celery can be configured to use less than the
default number of worker processes (which equals the number of available cores). To use only
a single core only, you can alter the worker service command in the `docker-compose.yml` and
add a `--concurrency` argument:

```
command: python -m celery -A nomad.processing worker -l info --concurrency=1 -Q celery
```

See also the [celery documentation](https://docs.celeryproject.org/en/stable/userguide/workers.html#id1).

#### Limiting the use of threads

You can also reduce the usable threads that Python packages based on OpenMP could use to
reduce the threads that might be spawned by a single worker process. Simply set the `OMP_NUM_THREADS`
environment variable in the worker container in your `docker-compose.yml`:

```
services:
    worker:
        ...
        environment:
            ...
            OMP_NUM_THREADS: 1
```

#### Limit CPU with docker

You can add a `deploy.resources.limits` section to the worker service in the `docker-compose.yml`:

```
services:
    worker:
        ...
        deploy:
            resources:
                limits:
                    cpus: '0.50'
```

The number refers to the percentage use of a single CPU core.
See also the [docker-compose documentation](https://docs.docker.com/compose/compose-file/compose-file-v3/#resources).


## NOMAD Oasis FAQ

#### Why use an Oasis?
There are three reasons:

- You want to manage data in private,
- you want local availability of public NOMAD data without network restrictions,
- or you want to manage large amounts of low quality and preliminary data that is not intended for publication.

#### How to organize data in NOMAD Oasis?

##### How can I categorize or label my data?
Currently, NOMAD supports the following mechanism to organize data:

- Data always belongs to one upload, which has a main author, and is assigned various metadata, like the time the upload was created, a custom human friendly upload name, etc.
- Data can be assigned to multiple independent datasets.
- Data can hold a proprietary id called “external_id”.
- Data can be assigned to multiple coauthors in addition to the main author.
- Data can be filtered and searched in various ways.

#####  Is there some rights-based visibility?
An upload can be viewed at any time by the main author, the upload coauthors, and the
upload reviewers. Other users will not be able to see the upload until it is published and
the embargo (if requested) has been lifted. The main author decides which users to register as coauthors and reviewers, when to publish and if it should be published with an embargo.
The main author and upload coauthors are the only users who can edit the upload while it is in staging.

#### How to share data with the central NOMAD?

Keep in mind, it is not entirely clear, how we will do this.

##### How to designate Oasis data for publishing to NOMAD?
Currently, you should use one of the organizational mechanism to designate data for being published to NOMAD. You can use a dedicated dataset for publishable data.

##### How to upload?

We will probably provide functionality in the API of the central NOMAD to upload data from an Oasis to the central NOMAD and provide the necessary scripts and detailed instructions. Most likely the data that is uploaded to the central NOMAD can be selected via a search query. Therefore, using a dedicated dataset, would be an easy to select criteria.

#### How to maintain an Oasis installation?

##### How to install a NOMAD Oasis?
Follow our guide: https://nomad-lab.eu/prod/v1/docs/ops.html#operating-a-nomad-oasis

##### How do version numbers work?
There are still a lot of things in NOMAD that are subject to change. Currently, changes in the minor version number (0.x.0) designate major changes that require data migration. Changes in the patch version number (0.7.x) just contain minor changes and fixes and do not require data migration. Once we reach 1.0.0, NOMAD will use the regular semantic versioning conventions.

##### How to upgrade a NOMAD Oasis?
When we release a new version of the NOMAD software, it will be available as a new Docker image with an increased version number. You simply change the version number in your docker-compose.yaml and restart.

##### What about major releases?
Going from NOMAD 0.7.x to 0.8.x will require data migration. This means the layout of the data has changed and the new version cannot be used on top of the old data. This requires a separate installation of the new version and mirroring the data from the old version via NOMAD’s API. Detailed instructions will be made available with the new version.

##### How to move data between installations?
With the release of 0.8.x, we will clarify how to move data between installations. (See last question)