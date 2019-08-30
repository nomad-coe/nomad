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
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d api
```

The different overrides are:
- *.prod.yml, production (to run the necessary databases for kubenetes deployments)
- *.override.yml, development (development configuration, will be automatically used by docker-compose)
- *.develk.yml, like development but also runs ELK
- *.example.yml, an example production configuration

To run your own NOMAD mirror or oasis, you should override the `example.yml` to fit your needs.
Within the overrides you can configure the NOMAD containers (api, worker, gui).

### Configuring the API and worker

The API and worker can be configured through a nomad.yaml file or with environment
variables. A nomad.yaml files needs to be mounted to `/app/nomad.yaml`. There are
[several options](https://docs.docker.com/compose/environment-variables/) to set
environment variables in docker-compose.

An example file can be found under `ops/docker-compose/nomad/example/nomad.yaml`.

### Configuring the GUI

This encompasses to parts a `env.js` file for the client side GUI code. And an
`nginx.conf` for the web server running the GUI (and reverse proxying the API).

Both need to be mounted to the respective containers. The `env.js` under `/app/nomad/env.js`
and the `nginx.conf` under `/etc/nginx/conf.d/default.conf`.

Example files can be found here `ops/docker-compose/nomad/example/`.


### Example *.prod.yml override:

Finally the full override for the example deployment based on the example config files
can be found here: `ops/docker-compose/nomad/docker-compose.example.yml`.

To simply run a fully fresh nomad:

```
git clone
cd nomad/ops/docker-compose/nomad/
docker login
docker-compose -f docker-compose.yml -f docker-compose.example.yml up
```

If everything goes well, NOMAD should be available at `http://your-host/my_example_nomad/gui/`.

We recommend to either change the nginx.conf to use SSL or put this behind a reverse-proxy
that supports SSL. If you use a reverse-proxy, you should disable any buffering to support
large downloads/uploads. We recommend to change the volume `nomad_data` to a bind mount.