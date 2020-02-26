# Operating a NOMAD OASIS

The following describes the simplest way to run your own NOMAD.

## What is an OASIS

Originally NOMAD is a service run at Max-Planck's compute facility in Garching, Germany.
However, the NOMAD software is Open-Source and everybody can run it. We call any service that
uses NOMAD software independently a *NOMAD OASIS*.

While there are several use cases that require different setups, this documentations
describes the simples NOMAD OASIS setup possible. It will allow you to use NOMAD to
manage research data locally, while using NOMAD's central user-management and its users.

## Pre-requisites

NOMAD software is distributed as a set of docker containers. Further, other services
that can be run with docker are required. Further, we use docker-compose to setup
all necessary container in the simples possible manner.

You will need a single computer, with **docker** and **docker-compose** installed.

The following will run all necessary services with docker. These comprise: a **mongodb**
database, an **elasticsearch**, a **rabbitmq** distributed task queue, the NOMAD **app**,
NOMAD **worker**, and NOMAD **gui**. Refer to this [introduction](/app/docs/introduction.html#architecture)
to learn what each service does and why it is necessary.

There is also some information you need to configure your NOMAD OASIS:
- The hostname for the machine you run NOMAD on. This is important for redirects between
your OASIS and the central NOMAD user-management and to allow your users to upload files (via GUI or API).
Your machine needs to be accessible under this hostname from the public internet. The host
name needs to be registered with the central NOMAD in order to configure the central user-
management correctly.
- A NOMAD account that acts as an admin account for your OASIS. This account must be declared
to the central NOMAD as an OASIS admin in order to give it the necessary rights in the central user-
management.

## Configuration

All docker container are configured via docker-compose an the respective `docker-compose.yaml` file.
Further, we will need to mount some configuration files to configure the NOMAD services within
their respective containers.

Please [write us](mailto:webmaster@nomad-coe.eu) to register your NOMAD account as an OASIS
admin and to register your hostname. Please replace the indicated configuration items with
the right information.

There are three files to configure:
- `docker-compose.yaml`
- `nomad.yaml`
- `env.js`
- `nginx.conf`

In this example, we have all files in the same directory (the directory we also work from).
You can download examples files from
[here](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/tree/master/ops/docker-compose/nomad-oasis/).

### Docker compose

The most basic `docker-compose.yaml` to run an OASIS looks like this:

```yaml
version: '3.4'

x-common-variables: &nomad_backend_env
    NOMAD_RABBITMQ_HOST: rabbitmq
    NOMAD_ELASTIC_HOST: elastic
    NOMAD_MONGO_HOST: mongo

services:
    # broker for celery
    rabbitmq:
        restart: always
        image: rabbitmq:3.7.17
        container_name: nomad_oasis_rabbitmq
        environment:
            - RABBITMQ_ERLANG_COOKIE=SWQOKODSQALRPCLNMEQG
            - RABBITMQ_DEFAULT_USER=rabbitmq
            - RABBITMQ_DEFAULT_PASS=rabbitmq
            - RABBITMQ_DEFAULT_VHOST=/
        volumes:
            - nomad_oasis_rabbitmq:/var/lib/rabbitmq

    # the search engine
    elastic:
        restart: always
        image: docker.elastic.co/elasticsearch/elasticsearch:6.3.2
        container_name: nomad_oasis_elastic
        volumes:
            - nomad_oasis_elastic:/usr/share/elasticsearch/data

    # the user data db
    mongo:
        restart: always
        image: mongo:4
        container_name: nomad_oasis_mongo
        environment:
            - MONGO_DATA_DIR=/data/db
            - MONGO_LOG_DIR=/dev/null
        volumes:
            - nomad_oasis_mongo:/data/db
        command: mongod --logpath=/dev/null # --quiet

    # nomad worker (processing)
    worker:
        restart: always
        image: gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair:stable
        container_name: nomad_oasis_worker
        environment:
            <<: *nomad_backend_env
            NOMAD_SERVICE: nomad_oasis_worker
        links:
            - rabbitmq
            - elastic
            - mongo
        volumes:
            - nomad_oasis_files:/app/.volumes/fs
            - ./nomad.yaml:/app/nomad.yaml
        command: python -m celery worker -l info -A nomad.processing -Q celery,calcs,uploads

    # nomad app (api)
    app:
        restart: always
        image: gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair:stable
        container_name: nomad_oasis_app
        environment:
            <<: *nomad_backend_env
            NOMAD_SERVICE: nomad_oasis_app
        links:
            - rabbitmq
            - elastic
            - mongo
        volumes:
            - nomad_oasis_files:/app/.volumes/fs
            - ./nomad.yaml:/app/nomad.yaml
        command: python -m gunicorn.app.wsgiapp -w 4 --log-config ops/gunicorn.log.conf -b 0.0.0.0:8000 --timeout 300 nomad.app:app

    # nomad gui (serving the js)
    gui:
        restart: always
        image: gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair/frontend:stable
        container_name: nomad_oasis_gui
        command: nginx -g 'daemon off;'
        volumes:
            - ./nginx.conf:/etc/nginx/conf.d/default.conf
            - ./env.js:/app/nomad/env.js
        links:
            - app
        ports:
            - 80:80
        command: ["./run.sh", "/example-nomad"]

volumes:
    nomad_oasis_mongo:
    nomad_oasis_elastic:
    nomad_oasis_rabbitmq:
    nomad_oasis_files:
```

There are no mandatory changes necessary.

A few things to notice:
- All services use docker volumes for storage. This could be changed to host mounts.
- It mounts three configuration files that need to be provided (see below): `nomad.yaml`, `nginx.conf`, `env.js`.
- The only exposed port is `80`. This could be changed to a desired port if necessary.
- The NOMAD images are pulled from our gitlab in Garching, the other services use images from a public registry (*dockerhub*).
- All container will be named `nomad_oasis_*`. These names can be used to later reference the container with the `docker` cmd.
- The NOMAD images we use are tagged `stable`. This could be replaced with concrete version tags.
- The services are setup to restart `always`, you might want to change this to `no` while debugging errors to prevent
indefinite restarts.

### nomad.yaml

NOMAD app and worker read a `nomad.yaml` for configuration.

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
```

You need to change:
- Replace `your-host` and admin credentials respectively.
- `api_base_path` defines the path under with the app is run. It needs to be changed, if you use a different base path.

A few things to notice:
- Be secretive about your admin credentials; make sure this file is not publicly readable.

### env.js

The GUI also has a config file, called `env.js` with a similar function than `nomad.yaml`.

```js
window.nomadEnv = {
  'appBase': '/nomad-oasis/',
  'keycloakBase': 'https://repository.nomad-coe.eu/fairdi/keycloak/auth/',
  'keycloakRealm': 'fairdi_nomad_prod',
  'keycloakClientId': 'nomad_public',
  'debug': false,
};
```

You need to change:
- `appBase` defines the base path again. It needs to be changed, if you use a different base path.

### nginx.conf

The GUI container also serves as a proxy that forwards request to the app container. The
proxy is an nginx server and needs a configuration similar to this:

```
server {
    listen        80;
    server_name   <your-host>;

    location /nomad-oasis {
        proxy_set_header Host $host;
        proxy_pass_request_headers on;
        proxy_pass http://app:8000;
    }

    location /nomad-oasis/gui {
        root      /app/;
        rewrite ^/nomad-oasis/gui/(.*)$ /nomad/$1 break;
        try_files $uri /nomad-oasis/gui/index.html;
    }

    location /nomad-oasis/gui/service-worker.js {
        add_header Last-Modified $date_gmt;
        add_header Cache-Control 'no-store, no-cache, must-revalidate, proxy-revalidate, max-age=0';
        if_modified_since off;
        expires off;
        etag off;
        root      /app/;
        rewrite ^/nomad-oasis/gui/service-worker.js /nomad/service-worker.js break;
    }

    location /nomad-oasis/api/uploads {
        client_max_body_size 35g;
        proxy_request_buffering off;
        proxy_set_header Host $host;
        proxy_pass_request_headers on;
        proxy_pass http://app:8000;
    }

    location /nomad-oasis/api/raw {
        proxy_buffering off;
        proxy_set_header Host $host;
        proxy_pass_request_headers on;
        proxy_pass http://app:8000;
    }
}
```

You need to change:
- Replace `<your-host>`

A few things to notice:
- It configures the base path (`nomad-oasis`) at multiple places. It needs to be changed, if you use a different base path.
- You can use the server to server additional content if you like.
- `client_max_body_size` sets a limit to the possible upload size.
- If you operate the GUI container behind another proxy, keep in mind that your proxy should not buffer requests/responses to allow streaming of large requests/responses for `../api/uploads` and `../api/raw`.

## Starting and stopping

If you prepared the above files, simply use the usual `docker-compose` commands to start everything.
In the beginning and for debugging problems, it is recommended to start services separately:
```
docker-compose up -d mongodb elastic rabbitmq
docker-compose up app worker gui
```

The `-d` option runs container in the background as *daemons*. Later you can run all at once:
```
docker-compose up -d
```

You can also use docker to stop and remove faulty containers that run as *daemons*:
```
docker stop nomad_oasis_app
docker rm nomad_oasis_app
```

If everything works, the gui should be available under:
```
http://<your host>/nomad-oasis/gui/
```

If you run into troubles, use the dev-tools of you browser to check the javascript logs
or monitor the network traffic for HTTP 500/400/404/401 responses.

To see if at least the api works, check
```
http://<your host>/nomad-oasis/alive
http://<your host>/nomad-oasis/api/info
```

To see logs or 'go into' a running container, you can access the individual containers
with their names and the usual docker commands:

```
docker logs nomad_oasis_app
```

```
docker exec -ti nomad_oasis_app /bin/bash
```

If you want to report problems with your OASIS. Please provide the logs for
- nomad_oasis_app
- nomad_oasis_worker
- nomad_oasis_gui