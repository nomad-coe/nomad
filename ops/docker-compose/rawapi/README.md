Contains all files to run the rawapi with docker-compose on a single docker, both
for development and production.

We use docker-compose overrides to modify config for development and production. Example:
```
docker-compose -f docker-compose.yml -f docker-compose.prod.yml up -d
```

The different overrides are:
- *.prod.yml, production (currently on enc-preprocessing-nomad.esc)

You have to create a .env file with the variable `RAW_FILE_DIR` to determine the
directory with the rawfiles: `<RAW_FILE_DIR>/data/R*.zip` .