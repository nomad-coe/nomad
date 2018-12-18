Contains all files to run nomad with docker-compose on a single docker, both
for development and production.

We use docker-compose overrides to modify config for development and production. Example:
```
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d api
```

The different overrides are:
- *.prod.yml, production (currently on enc-preprocessing-nomad.esc)
- *.override.yml, development (will be automatically used by docker-compose)
- *.develk.yml, like development but also runs ELK

The .env file contains some additional config and secrets. The development secrets do
not matter and are in the git (.env_development) and are replaced by real secret on
the production machine.