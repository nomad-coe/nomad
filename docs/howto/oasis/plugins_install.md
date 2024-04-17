# How to mount a plugin into a NOMAD Oasis
[Plugins](../customization/plugins.md#register-your-plugin) allow the customization of a
NOMAD deployment in terms of which parsers, schemas, and normalizers are included or excluded.
In the following we will show to how to mount specific plugins in a NOMAD Oasis.

The NOMAD docker image adds the folder `/app/plugins` to the `PYTHONPATH`. You simply have
to add the *plugin metadata* to your Oasis' `nomad.yaml` and mount your code into the `/app/plugins`
directory via the volumes section of the `app` and `worker` services in your `docker-compose.yaml`.

For example, you can do this by adding an extension to the `docker-compose.yaml`, e.g. a file called
`docker-compose.plugins.yaml`. Assuming you cloned the example plugins above into the Oasis folder as
`./nomad-schema-plugin-example`, `./nomad-parser-plugin-example` and `./nomad-normalizer-plugin-example`,
your `docker-compose.plugins.yaml` should look like this:

```yaml
services:
  worker:
    volumes:
      - ./nomad-schema-plugin-example/nomadschemaexample:/app/plugins/nomadschemaexample
      - ./nomad-parser-plugin-example/nomadparserexample:/app/plugins/nomadparserexample
      - ./nomad-normalizer-plugin-example/nomadparserexample:/app/plugins/nomadparserexample
  app:
    volumes:
      - ./nomad-schema-plugin-example/nomadschemaexample:/app/plugins/nomadschemaexample
      - ./nomad-parser-plugin-example/nomadparserexample:/app/plugins/nomadparserexample
      - ./nomad-normalizer-plugin-example/nomadparserexample:/app/plugins/nomadparserexample
```

You have to tell docker that there are now two compose files. This can be done via the
`COMPOSE_FILE` environment variable. This is how you can start the Oasis with the plugins:

```sh
export COMPOSE_FILE=docker-compose.yaml:docker-compose.plugins.yaml
docker compose up -d
```

Here is a complete Oasis setup [nomad-oasis-with-plugins.zip](../../assets/nomad-oasis-with-plugins.zip).
Simply download, extract, and start like any other Oasis:

```sh
unzip nomad-oasis-with-plugins.zip
cd nomad-oasis-with-plugins
sudo chown -R 1000 .volumes
sudo chown -R 1000 nomad-schema-plugin-example
sudo chown -R 1000 nomad-parser-plugin-example
sudo chown -R 1000 nomad-normalizer-plugin-example
export COMPOSE_FILE=docker-compose.yaml:docker-compose.plugins.yaml
docker compose pull
docker compose up -d
curl localhost/nomad-oasis/alive
```

!!! warning "Attention"
    It is important to set up the correct user rights for your volumes and
    plugins. Our default `docker-compose` setup uses the user `1000` in group
    `1000` to run the services, this is the reason for the `chown` commands
    above that ensure that the processes have access to the data stored in
    volumes and in the plugins. If you use another user/group to run the docker
    services, update the commands accordingly.

Read the [Oasis install guide](install.md) for more details.

### Install PyPI/pip package

If the plugin is published on PyPI, you can simply install it with pip. If the
plugin was published to our MPCDF GitLab registry, you have to use the `--index-url`
parameter:

```
pip install nomad-example-schema-plugin --index-url https://gitlab.mpcdf.mpg.de/api/v4/projects/2187/packages/pypi/simple
```

Installing via pip works for NOMAD developers, but how to pip install into an Oasis?
The package could either be installed when NOMAD is started or via
a customized docker image.

!!! warning "Attention"
    We still need to implement that configured plugins, if not already installed,
    get automatically installed during NOMAD start.

You can build a custom NOMAD docker image that has your packages already installed.
Here is an example `Dockerfile`:

```Dockerfile
--8<-- "examples/plugins/schema/Dockerfile"
```

The image can be build like this:

```
docker build -t nomad-with-plugins .
```
