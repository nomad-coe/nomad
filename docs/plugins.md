Plugins allow you to add Python-based functionality to NOMAD without a custom build
NOMAD image or package. We support different kinds of plugins:

- Python **schema**, read also [Python schema documentation](schema/python.md).
- **parser**, read also [parser development documentation](develop/parser.md).
- **normalizer** (coming soon...)
- additional custom **APIs** (coming soon...)

## Develop a plugin

We provide template projects on GitHub. You can fork these projects and follow the
instructions in their `README.md`. These instructions will give you everything you
need to run your plugin as a plugin developer:

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

The file `nomad_plugin.yaml` contains the installation independent plugin metadata:

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

Now follow the instructions for one of our examples and try for yourself:

- [schema plugin](https://github.com/nomad-coe/nomad-schema-plugin-example)
- [parser plugin](https://github.com/nomad-coe/nomad-parser-plugin-example)

## Add a plugin to your NOMAD (Oasis)

To add a plugin, you need to add the *plugin metadata* to `nomad.yaml` and you need
to add the *plugin code* to the `PYTHONPATH` of your NOMAD. The `nomad.yaml` needs to be
edited manually in the usual way. There are several ways to
add *plugin code* to a NOMAD installation:

- adding it to the `PYTHONPATH`, if you run nomad directly, e.g. as a developer
- mounting plugin sources into NOMAD (Oasis) containers, e.g. as an Oasis admin (coming soon...)
- via python packages (coming soon...)
- via github projects (coming soon...)

## Publish a plugin

coming soon...

We plan to provide a plugin registry that allows you to publish your plugin's *metadata*.
This can then be used to simplify plugin management within a NOMAD installation.