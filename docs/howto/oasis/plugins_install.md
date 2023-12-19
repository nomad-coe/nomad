# How to install plugins

Plugins allow you to add Python-based functionality to NOMAD without a custom build
NOMAD image or release. Plugins can be installed at NOMAD start-up time. Therefore, you can
configure each NOMAD (Oasis) with a different custom set of plugins or disable unnecessary
plugins.

We support different kinds of plugins:

- Python **schema**
- **parser**
- **normalizer**
- additional custom **APIs** (coming soon...)

## Develop a plugin

We provide template projects on GitHub. You can fork these projects and follow the
instructions in their `README.md`. These instructions will give you everything you
need to run and test your plugin as a plugin developer.

The following sections explain how to add plugins to a NOMAD installation.<br />
Dedicated Explanation sections provide more background information on [what is a schema](../../explanation/data.md#schema) and [what is a parser](../../explanation/processing.md#schemas-parsers-plugins)

- [schema plugin](https://github.com/nomad-coe/nomad-schema-plugin-example){:target="_blank"}
- [parser plugin](https://github.com/nomad-coe/nomad-parser-plugin-example){:target="_blank"}
- [normalizer plugin](https://github.com/nomad-coe/nomad-normalizer-plugin-example.git){:target="_blank"}


## Plugin anatomy

A plugin usually consist of the *plugin code* (a Python package) and
*plugin metadata*. The installation independent *plugin metadata* (e.g. name, description, python package, etc.)
can be defined in a `nomad_plugin.yaml` that is part of the *plugin code*.
The installation dependent *plugin metadata* (e.g. plugin key, order and priority, parser matching rules, etc.)
is added to the [`nomad.yaml` file](../develop/setup.md#nomadyaml) of the NOMAD installation.

Here is the project layout of the schema example:

```
my-nomad-schema
├── nomadschemaexample
│   ├── __init__.py
│   ├── nomad_plugin.yaml
│   └── schema.py
├── tests
│   ├── data
│   │   └── test.archive.yaml
│   └── test_schema.py
├── LICENSE
├── README.md
├── nomad.yaml
└── requirements.txt
```

## Plugin code

The directory `nomadschemaexample` is our Python package *plugin code*. In this case,
it contains a simple `schema.py`. Read the [Schema plugin documentation](../customization/plugins_dev.md#develop-a-schema-plugin)
for more details:

```python
{{ file_contents('examples/plugins/schema/nomadschemaexample/schema.py') }}
```

## Plugin metadata

The file `nomad_plugin.yaml` contains the installation independent *plugin metadata*:

```yaml
{{ file_contents('examples/plugins/schema/nomadschemaexample/nomad_plugin.yaml') }}
```

The metadata contains the `plugin_type` (e.g. `schema` or `parser`). The rest of the
yaml will depend on the type and the underlying metadata model. For schemas there are only
descriptive metadata like `name` or `description` as schemas do not contain any technical
metadata that is necessary to use them. See below for a reference of the *plugin metadata*
models.

The file `nomad.yaml` shows how to add the plugin to a nomad installation. As a plugin
developer you have [installed our Python package](../programmatic/pythonlib.md) and can run the `nomad parse`
command as your "installation" to try your schema:

```yaml
{{ file_contents('examples/plugins/schema/nomad.yaml') }}
```

Plugins are defined under the `plugins` key. This consists of `include` (or `exclude`) to
select a subset of all plugins defined under `options`. The `options` given in the
`nomad.yaml` will be merged with all the default "plugins" that come with NOMAD (mostly the NOMAD parsers at the moment).

In this example, we disable all default plugins by just including our `schemas/example`.
The `options` field can be used to add define the *code* and installation independent *metadata*
via the `python_package` key. But you can overwrite or add more *metadata* keys here as well:
each `options` entry with `python_package` will be merged with the data in the package's
`nomad_plugin.yaml`.

Please note that `python_package` is the name of a Python package and not a path to the
code. This also means that the package has to be in your `PYTHONPATH` (see below).


Now follow the instructions for one of our examples and try for yourself:

- [schema plugin](https://github.com/nomad-coe/nomad-schema-plugin-example){:target="_blank"}
- [parser plugin](https://github.com/nomad-coe/nomad-parser-plugin-example){:target="_blank"}
- [normalizer plugin](https://github.com/nomad-coe/nomad-normalizer-plugin-example){:target="_blank"}


# Publish a plugin

!!! warning "Attention"
    The processes around publishing plugins and using plugins of others are still
    worked on. The "best" practices mentioned here are preliminary.

## Create a (GitHub) project

If you forked from our examples, you already have a GitHub project. Otherwise, you
should create one. This allows others to get your plugin sources or initiate communication
via issues or pull requests.

These are good names for plugin projects, depending on if you maintain one or more
plugins in a project (a project can contain multiple modules with multiple
`nomad-plugin.yaml` files and therefore multiple plugins):

- nomad-<yourname\>-plugin
- nomad-<yourname\>-plugins

!!! note
    If you develop a plugin in the context of **FAIRmat** or the **NOMAD CoE**, put your
    plugin projects in the respective GitHub organization for [FAIRmat](https://github.com/fairmat-nfdi){:target="_blank"}
    and the [NOMAD CoE](https://github.com/nomad-coe){:target="_blank"}. Here, the naming convention above is binding.

Your plugin projects should follow the layout of our example projects.

## Different forms of plugin distribution

- **source code**: Mounting plugin code into a NOMAD (Oasis) installation. This is described above and only
the plugin source code is needed.
- **built-in**: Plugins that are directly maintained by NOMAD as distributed as part of
the NOMAD docker images. The Python code for those plugins is already installed, you only need
to configure NOMAD to use the plugins (or not).
- **PyPI/pip package**: Plugin projects can be published as PyPI/pip packages. Those
packages can then be installed either during NOMAD start-up (not implemented yet) or
when building a customized docker images (see [below](#pypipip-package)).

Independent of the form of distribution, you'll still need to add the plugin to
your configuration as explained above.

## PyPI/pip package

Learn from the PyPI documentation how to [create a package for PyPI](https://packaging.python.org/en/latest/tutorials/packaging-projects/){:target="_blank"}.
We recommend to use the `pyproject.toml`-based approach. Here is an example `pyproject.toml` file:

```toml
--8<-- "examples/plugins/schema/pyproject.toml"
```

The package can be build like this:
```
pip install build
python -m build --sdist
```

Learn from the PyPI documentation how to [publish a package to PyPI](https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives){:target="_blank"}.
If you have access to the MPCDF GitLab and NOMAD's presence there, you can also
use the `nomad-FAIR` registry:

```
pip install twine
twine upload \
    -u <username> -p <password> \
    --repository-url https://gitlab.mpcdf.mpg.de/api/v4/projects/2187/packages/pypi \
    dist/nomad-example-schema-plugin-*.tar.gz
```

## Register your plugin

!!! warning "Attention"
    This is work in progress. We plan to provide a plugin registry that allows you to
    publish your plugin's *metadata*. This will then be used to simplify plugin management
    within a NOMAD installation.

    The built-in plugins can already be found in the [documentation reference](../../reference/plugins.md).

# Add a plugin to your NOMAD

Adding a plugin, depends on the form of plugin distribution and how you run NOMAD.
Eventually, you need to add the *plugin metadata* to `nomad.yaml` (see above) and you need
to add the *plugin code* to the `PYTHONPATH`. The `nomad.yaml` needs to be
edited manually in the usual ways. There are several ways to add *plugin code*.

## Built-in plugins

Those are already part of the NOMAD sources or NOMAD docker images. You only need
to configure them in your `nomad.yaml`.

## Add to Python path

When you run NOMAD as a developer, simply add the plugin directory to the `PYTHONPATH` environment variable.
When you start the application (e.g. `nomad admin run appworker`), Python will find your code when NOMAD
imports the `python_package` given in the `plugins.options` of your `nomad.yaml`.

## Mount into a NOMAD Oasis

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

## Install PyPI/pip package

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
