Plugins allow you to add Python-based functionality to NOMAD without a custom build
NOMAD image or release. Plugins can be installed at NOMAD start-up time. Therefore, you can
configure each NOMAD (Oasis) with a different custom set of plugins or disable unnecessary
plugins.

We support different kinds of plugins:

- Python **schema**, read also [Python schema documentation](schema/python.md).
- **parser**, read also [parser development documentation](develop/parser.md).
- **normalizer** (coming soon...)
- additional custom **APIs** (coming soon...)

## Develop a plugin

We provide template projects on GitHub. You can fork these projects and follow the
instructions in their `README.md`. These instructions will give you everything you
need to run and test your plugin as a plugin developer.
The following sections here contain more background information and explain how to
add plugins to a NOMAD installation.

- [schema plugin](https://github.com/nomad-coe/nomad-schema-plugin-example)
- [parser plugin](https://github.com/nomad-coe/nomad-parser-plugin-example)

### Plugin anatomy

A plugin usually consist of the *plugin code* (a Python package) and
*plugin metadata*. The installation independent *plugin metadata* (e.g. name, description, python package, etc.)
can be defined in a `nomad_plugin.yaml` that is part of the *plugin code*.
The installation dependent *plugin metadata* (e.g. plugin key, order and priority, parser matching rules, etc.)
is added to the `nomad.yaml` of the NOMAD installation.

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

### Plugin code

The directory `nomadschemaexample` is our Python package *plugin code*. In this case,
it contains a simple `schema.py`. Read the [Python schema documentation](schema/python.md)
for more details:

```python
{{ file_contents('examples/plugins/schema/nomadschemaexample/schema.py') }}
```

### Plugin metadata

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
developer you have [installed our Python package](./pythonlib.md) and can run the `nomad parse`
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


{{pydantic_model('nomad.config.plugins.Schema', heading='### Parser plugin metadata')}}
{{pydantic_model('nomad.config.plugins.Parser', heading='### Schema plugin metadata')}}

Now follow the instructions for one of our examples and try for yourself:

- [schema plugin](https://github.com/nomad-coe/nomad-schema-plugin-example)
- [parser plugin](https://github.com/nomad-coe/nomad-parser-plugin-example)

## Add a plugin to your NOMAD

To add a plugin, you need to add the *plugin metadata* to `nomad.yaml` (see above) and you need
to add the *plugin code* to the `PYTHONPATH` of your NOMAD. The `nomad.yaml` needs to be
edited manually in the usual way. There are several ways to
add *plugin code* to a NOMAD installation.

### Development setup of NOMAD

Simply add the plugin directory to the `PYTHONPATH` environment variable. When you start
the application (e.g. `nomad admin run appworker`), Python will find your code when NOMAD
imports the `python_package` given in the `plugins.options` of your `nomad.yaml`.

### NOMAD Oasis

The NOMAD docker image adds the folder `/app/plugins` to the `PYTHONPATH`. You simply have
to add the *plugin metadata* to your Oasis' `nomad.yaml` and mount your code into the `/app/plugins`
directory via the volumes section of the `app` and `worker` services in your `docker-compose.yaml`.

For example, you can do this by adding an extension to the `docker-compose.yaml`, e.g. a file called
`docker-compose.plugins.yaml`. Assuming you cloned the example plugins above into the Oasis folder as
`./nomad-schema-plugin-example` and `./nomad-parser-plugin-example`,
your `docker-compose.plugins.yaml` should look like this:

```yaml
services:
  worker:
    volumes:
      - ./nomad-schema-plugin-example/nomadschemaexample:/app/plugins/nomadschemaexample
      - ./nomad-parser-plugin-example/nomadparserexample:/app/plugins/nomadparserexample
  app:
    volumes:
      - ./nomad-schema-plugin-example/nomadschemaexample:/app/plugins/nomadschemaexample
      - ./nomad-parser-plugin-example/nomadparserexample:/app/plugins/nomadparserexample
```

You have to tell docker that there are now two compose files. This can be done via the
`COMPOSE_FILE` environment variable. This is how you can start the Oasis with the plugins:

```sh
export COMPOSE_FILE=docker-compose.yaml:docker-compose.plugins.yaml
docker compose up -d
```

Here is a complete Oasis setup [nomad-oasis-with-plugins.zip](assets/nomad-oasis-with-plugins.zip).
Simply download, extract, and start like any other Oasis:

```sh
unzip nomad-oasis-with-plugins.zip
cd nomad-oasis-with-plugins
sudo chown -R 1000 .volumes
export COMPOSE_FILE=docker-compose.yaml:docker-compose.plugins.yaml
docker compose pull
docker compose up -d
curl localhost/nomad-oasis/alive
```

Read the [Oasis documentation](oasis.md) for more details.

### Other means

- via python packages (coming soon...)
- via github projects (coming soon...)

## Publish a plugin

coming soon...

We plan to provide a plugin registry that allows you to publish your plugin's *metadata*.
This can then be used to simplify plugin management within a NOMAD installation.