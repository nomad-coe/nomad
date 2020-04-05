## Run (dev) infrastructure with docker compose

You can all necessary databases and other infrastructure with [docker-compose](https://docs.docker.com/compose/)
on a single node/computer that supports docker and docker-compse.

### How we use docker-compose

You can use docker-compose to run all necessary databases with one single docker-compose configuration.
 The `docker-compose.yml` defines all the different container.

We use docker-compose overrides to extend a base configuration for different scenarios.
Example docker-compose usage:

```
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d mongo rabbitmq elastic
```

The different overrides are:
- *.prod.yml, production (to run the necessary databases for kubenetes deployments)
- *.override.yml, development (development configuration, will be automatically used by docker-compose)
- *.develk.yml, like development but also runs ELK

To run nomad on top use the nomad command:
```
nomad admin run appworker
```