# How to install plugins into a NOMAD Oasis

[Plugins](../plugins/plugins.md) allow the customization of a NOMAD deployment in terms of which apps, normalizers, parsers and schema packages are available. In the following we will show to how to install plugins into a NOMAD Oasis.

## Option 1: Mount the plugin code

The NOMAD docker image adds the folder `/app/plugins` to the `PYTHONPATH`. This means that you can mount your code into the `/app/plugins` directory via the volumes section of the `app` and `worker` services in your `docker-compose.yaml`.

For example, you can do this by adding an extension to the `docker-compose.yaml`, e.g. a file called `docker-compose.plugins.yaml`. Assuming you have cloned three plugins into the Oasis folder as `./nomad-schema-plugin-example`, `./nomad-parser-plugin-example` and `./nomad-normalizer-plugin-example`,
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

You have to tell docker that there are now two compose files. This can be done via the `COMPOSE_FILE` environment variable. This is how you can start the Oasis with the plugins:

```sh
export COMPOSE_FILE=docker-compose.yaml:docker-compose.plugins.yaml
docker compose up -d
```

Here is a complete Oasis setup [nomad-oasis-with-plugins.zip](../../assets/nomad-oasis-with-plugins.zip). Simply download, extract, and start like any other Oasis:

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

## Option 2: Create a derived Docker image with plugin installed via `pip`

Instead of mounting the code into an existing image, you can also create a new, derived image which has your plugin installed as a `pip` package. For this you will to create a new `Dockerfile`, which runs the installation step. The basic idea is that your Dockerfile looks something like this:

```Dockerfile
--8<-- "examples/plugins/schema/Dockerfile"
```

The image can then be build like this:

```
docker build -t nomad-with-plugins .
```

Depending on how your plugin code is distributed, you have several options for the actual install steps:

1. Plugin published in PyPI:

    ```sh
    RUN pip install <package_name>`
    ```

2. Plugin code available in GitHub:

    ```sh
    RUN pip install git+https://<repository_url>
    ```

3. Plugin published in MPCDF GitLab registry:

    ```sh
    RUN pip install nomad-example-schema-plugin --index-url https://gitlab.mpcdf.mpg.de/api/v4/projects/2187/packages/pypi/simple
    ```

4. Copy plugin code from host machine:

    ```sh
    RUN pip install build

    COPY \
        nomadschemaexample \
        tests \
        README.md \
        LICENSE \
        pyproject.toml \
        .

    RUN python -m build --sdist

    RUN pip install dist/nomad-schema-plugin-example-*.tar.gz
    ```
