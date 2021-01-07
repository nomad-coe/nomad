## Run (dev) infrastructure components with docker compose

You can run all necessary databases and other infrastructure with [docker-compose](https://docs.docker.com/compose/)
on a single node/computer that supports docker and docker-compse.

### How we use docker-compose

You can use docker-compose to run all necessary databases with one single docker-compose configuration.
 The `docker-compose.yml` defines all the different container.

To run the infrastructure for a typical development environment (where you need mongodb,
elastic, and rabbitmq), simply run:

```sh
docker-compose up -d mongo elastic rabbitmq
```

To run nomad on top use the nomad command:

```sh
nomad admin run appworker
```
