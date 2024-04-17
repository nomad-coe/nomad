# How to write a plugin

The following sections explain how to write a plugin and add it to a NOMAD installation.
Dedicated explanation sections provide more background information on the types of plugins:
[schema](../../explanation/data.md#schema), [parser](../../explanation/processing.md#schemas-parsers-plugins)
and [normalizer](../../explanation/processing.md#normalizing).

## Plugin anatomy

A plugin usually consists of the *plugin code* (a Python package) and
*plugin metadata*. The installation **independent** *plugin metadata* (e.g. name, description, python package, etc.)
can be defined in a `nomad_plugin.yaml` that is part of the *plugin code*.
The installation **dependent** *plugin metadata* (e.g. plugin key, order and priority, parser matching rules, etc.)
are added to the [`nomad.yaml` file](../develop/setup.md#nomadyaml) of the NOMAD installation.

Here is the project layout of the schema example:
```
nomad-schema-plugin-example
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
<!-- TODO pyproject.toml, MANIFEST.in, setup.py are missing. requirements.txt is no longer there Additionally, we could adopt following a src structure, src/nomadschemaexample. -->

## Plugin code

The directory `nomadschemaexample` is our Python package *plugin code*, and contains `schema.py`:

```python
{{ file_contents('examples/plugins/schema/nomadschemaexample/schema.py') }}
```

Read the [schema plugin documentation](schemas.md#develop-a-schema-plugin)
for more details.

### Code Quality and Linting

While developing NOMAD plugins, we highly recommend using a Python linter, such as [Ruff](https://docs.astral.sh/ruff), to analyze and enforce coding standards in your plugin projects. This also ensures smoother integration and collaboration. Ruff is also included in the templates provided on Github.

## Plugin metadata

The file `nomad_plugin.yaml` contains the installation **independent** *plugin metadata*. The following is for the schema plugin example:

```yaml
{{ file_contents('examples/plugins/schema/nomadschemaexample/nomad_plugin.yaml') }}
```

The metadata contains the `plugin_type` (e.g. `schema`, `parser` or `normalizer`). The rest
will depend on the type and the underlying metadata model. For schemas there are only
descriptive metadata like `name` or `description` as schemas do not contain any technical
metadata that are necessary to use them. Refer to the *plugin metadata* models for
[schema](schemas.md#schema-plugin-metadata), [parser](parsers.md#parser-plugin-metadata)
and [normalizer](normalizers.md#normalizer-plugin-metadata).

One can specify which plugin to enable in a nomad installation in the `nomad.yaml` file:

```yaml
{{ file_contents('examples/plugins/schema/nomad.yaml') }}
```

Plugins are defined under the `plugins` key. This consists of `include` (or `exclude`) to
select (or ignore) a subset of all plugins. In this example, we disable all [built-in plugins](#different-forms-of-plugin-distribution) by only including `schemas/example` under `plugins`.
The `options` field can be used to define the plugin metadata.
This allows one to overwrite the metadata in the package's `nomad_plugin.yaml`.

Please note that `python_package` is the name of a Python package and not a path to the
code. This also means that the package has to be in your `PYTHONPATH` (see [Add a plugin to your NOMAD](#add-a-plugin-to-your-nomad)).


As a plugin developer you have [installed the NOMAD Python package](../programmatic/pythonlib.md)
and can run the `nomad parse <mainfile>` command to make sure installation is successful.
Now follow the instructions for [one of our examples](#develop-a-plugin) and try for yourself!

# Publish a plugin

!!! warning "Attention"
    The standard processes for publishing plugins and using plugins from other developers are still being worked out. The "best" practices mentioned in the following are preliminary.

## Create a (GitHub) project

If you forked from our examples, you already have a GitHub project. Otherwise, you
should create one. This allows others to get your plugin sources or initiate communication
via issues or pull requests.

!!! tip "Important"
    If you create a project from scratch, you should still follow the layout of our example projects.

We suggest the following naming convention for plugin projects:

- `nomad-<projectname>-plugin` (a single plugin)
- `nomad-<projectname>-plugins` (multiple plugins)

A project can contain multiple plugins if it has multiple modules with corresponding
`nomad-plugin.yaml` files.

!!! note
    If you develop a plugin in the context of [FAIRmat](https://github.com/fairmat-nfdi) or
    the [NOMAD CoE](https://github.com/nomad-coe), put your plugin projects in the
    corresponding GitHub organization. In these cases, the naming convention above is required.

## Different forms of plugin distribution

- **source code**: [Mounting plugin code into a NOMAD (Oasis) installation](../oasis/plugins_install.md#mount-plugin-into-a-nomad-oasis) -- only the plugin source code is needed.
- **built-in**: Plugins that are directly maintained by NOMAD as distributed as part of
the NOMAD docker images. The Python code for those plugins is already installed, you only need
to configure NOMAD to use the plugins (or not).
- **PyPI/pip package**: Plugin projects can be published as PyPI/pip packages. Those
packages can then be installed either during NOMAD start-up (not implemented yet) or
when building a customized docker images (see [PyPI/pip package](#pypipip-package) below).

Independent of the form of distribution, you will still need to add the plugin to
your configuration as explained in [previously](#plugin-metadata).

## PyPI/pip package

Learn from the PyPI documentation how to [create a package for PyPI](https://packaging.python.org/en/latest/tutorials/packaging-projects/){:target="_blank"}.
We recommend to use the `pyproject.toml`-based approach. Here is an example `pyproject.toml` file:

```toml
--8<-- "examples/plugins/schema/pyproject.toml"
```

The package can be built like this:
```
pip install build
python -m build --sdist
```

The PyPI documentation provides further information about how to [publish a package to PyPI](https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives){:target="_blank"}.
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

Adding a plugin depends on the type of plugin distribution and how you run NOMAD.
However, in all cases you will need to add the *plugin metadata* to `nomad.yaml` (see [above](#plugin-metadata)) and include the *plugin code* to the `PYTHONPATH`. There are several ways to add *plugin code*.

## Add to Python path

When you run NOMAD as a developer, simply add the plugin directory to the `PYTHONPATH` environment variable.
When you [run NOMAD](../develop/setup.md#run-nomad) (e.g. `nomad admin run appworker`), Python will find your code when NOMAD imports the `python_package` given in the `plugins.options` of your `nomad.yaml`.

