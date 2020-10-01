## Run (dev) infrastructure components with docker compose

You can run all necessary databases and other infrastructure with [docker-compose](https://docs.docker.com/compose/)
on a single node/computer that supports docker and docker-compse.

### How we use docker-compose

You can use docker-compose to run all necessary databases with one single docker-compose configuration.
 The `docker-compose.yml` defines all the different container.

To run the infrastructure for a typical development environment (where you need mongodb,
elastic, and rabbitmq), simply run:

```
docker-compose up -d mongo elastic rabbitmq
```

We use docker-compose overrides to extend a base configuration for different scenarios.
Example docker-compose usage for starting the production infrastructure:

```
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d mongo rabbitmq elastic keycloak elk
```

The different overrides are:
- .prod.yml, production (to run the necessary databases for kubenetes deployments)
- .override.yml, development (development configuration, will be automatically used by docker-compose)
- .develk.yml, like development but also runs ELK

To run nomad on top use the nomad command:

```
nomad admin run appworker
```
