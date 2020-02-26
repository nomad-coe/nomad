## Single Node NOMAD Deployment, Using Docker Compose

You can run NOMAD with [docker-compose](https://docs.docker.com/compose/) on a single node
that supports docker and docker-compse. The docker-compose files are part of the
[NOMAD source code](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR)
and can be found under `ops/docker-compose/nomad`.

### How we use docker-compose

You can use docker-compose to run all necessary databases, the API, the worker, and
gui with one single docker-compose configuration. The `docker-compose.yml` defines
all the different container.

We use docker-compose overrides to extend a base configuration for different scenarios.
Example docker-compose usage:

```
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d app
```

The different overrides are:
- *.prod.yml, production (to run the necessary databases for kubenetes deployments)
- *.override.yml, development (development configuration, will be automatically used by docker-compose)
- *.develk.yml, like development but also runs ELK

Within the overrides you can configure the NOMAD containers (app, worker, gui).

### Configuring the API and worker

The API and worker can be configured through a nomad.yaml file or with environment
variables. A nomad.yaml files needs to be mounted to `/app/nomad.yaml`. There are
[several options](https://docs.docker.com/compose/environment-variables/) to set
environment variables in docker-compose.

### Configuring the GUI

This encompasses to parts a `env.js` file for the client side GUI code. And an
`nginx.conf` for the web server running the GUI (and reverse proxying the API).

Both need to be mounted to the respective containers. The `env.js` under `/app/nomad/env.js`
and the `nginx.conf` under `/etc/nginx/conf.d/default.conf`. Usually the only configuration
item for the `nginx.conf` is the desired URL path prefix. You can generate an `nginx.conf`
with:

```
nomad admin nginx-conf --prefix /example-nomad
```
