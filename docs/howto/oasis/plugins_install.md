# How to install plugins into a NOMAD Oasis

[Plugins](../plugins/plugins.md) allow the customization of a NOMAD deployment in terms of which apps, normalizers, parsers and schema packages are available. In order for these customization to be activated, they have to be installed into an Oasis.

Oasis is controlled and run through a `docker-compose.yaml` file, which specifies the different software services and how they interact. Some of these services are using a Docker image that contains the actual NOMAD software. It is in this image where we will need to install any additional plugins with `pip install`.

The following sections contain some alternatives for achieving this, with the first option being the preferred one.

## Option 1: Create a new customized NOMAD Oasis distribution with your plugins

When initially starting to create a customized NOMAD Oasis distribution, it is strongly advised that you create a GitHub repository to persist your work, collaborate with coworkers and also to automate the building and distribution of your custom image. To streamline this process, we have created a [GitHub template repository](https://github.com/FAIRmat-NFDI/nomad-distribution-template) that helps with all of this. It can do the following for you:

- Plugins are controlled with a simple `plugins.txt` file where it is easy to install plugins from PyPI, Git repositories, local files, etc.
- The automatic pipeline will create a new Docker image for your Oasis. This image will also be stored on GitHub servers.
- Initial modifications to the `docker-compose.yaml` are done automatically so you can boot up the software directly.

To learn more, head over to the [template repository](https://github.com/FAIRmat-NFDI/nomad-distribution-template) and follow the instructions there.

## Option 2: Only create a customized Docker image

If you already have an existing NOMAD Oasis setup, or do not wish to use the template, you can also just create a new Docker image which has your plugin installed as a `pip` package. For this approach, you need to create a new `Dockerfile`, which runs the installation step on top of our default image. The basic idea is that your Dockerfile looks something like this:

```Dockerfile
FROM gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair:latest

# Install your plugin here, e.g.:
RUN pip install git+https://<repository_url>
```

Depending on how your plugin code is distributed, you have several options for the actual install steps:

1. Plugin published in PyPI:

    ```sh
    RUN pip install <package_name>
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

The customized image can then be built like this:

```
docker build -t nomad-with-plugins . --build-arg SETUPTOOLS_SCM_PRETEND_VERSION=<insert-nomad-version>
```

This will create a new image with the tag `nomad-with-plugins`, which you can use in your `docker-compose.yaml` file:

```yaml
#image: gitlab-registry.mpcdf.mpg.de/nomad-lab/nomad-fair:latest
image: nomad-with-plugins
```

## Option 3 (deprecated): Mount the plugin code directly into the container

!!! warning "Attention"
    This option only works with the old plugin mechanism that is based on `nomad_plugin.yaml` files instead of [Python entry points](https://setuptools.pypa.io/en/latest/userguide/entry_point.html).

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
