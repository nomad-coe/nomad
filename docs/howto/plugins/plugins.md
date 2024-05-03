# Get started with plugins

The main way to customize a NOMAD installation is through the use of **plugins**. A NOMAD Plugin is a Git repository that contains a Python package that an administrator can install into a NOMAD deployment to add custom features. This page contains the basics of how to create, develop and publish a NOMAD Plugin.

## Plugin anatomy

!!! tip
    We provide a [template repository](https://github.com/FAIRmat-NFDI/nomad-plugin-template) which you can use to create the initial plugin layout for you.

Plugin Git repositories should roughly follow this layout:

```txt
├── nomad-example
│   ├── src
|   │   ├── nomad_example
|   |   │   ├── apps
|   |   │   │   ├── __init__.py
|   |   │   ├── normalizers
|   |   │   │   ├── mynormalizer.py
|   |   │   │   ├── __init__.py
|   |   │   ├── schema_packages
|   |   │   │   ├── mypackage.py
|   |   │   │   ├── __init__.py
|   |   │   ├── parsers
|   |   │   │   ├── myparser.py
|   |   │   │   ├── __init__.py
│   ├── docs
│   ├── tests
│   ├── pyproject.toml
│   ├── LICENSE.txt
│   ├── README.md
```

We suggest using the following convention for naming the repository name and the plugin package:

 - repository name: `nomad-<plugin name>`
 - package name: `nomad_<plugin name>`

In the folder structure you can see that a single plugin can contain multiple types of customizations: apps, parsers, schema packages and normalizers. These are called a **plugin entry points** and you will learn more about them next.

## Plugin entry points

Plugin entry points represent different types of customizations that can be added to a NOMAD installation. The following plugin entry point types are currently supported:

 - [Apps](./apps.md)
 - [Normalizers](./parsers.md)
 - [Parsers](./parsers.md)
 - [Schema packages](./schema_packages.md)

Entry points contain **configuration**, but also a separate **resource**, which should live in a separate Python module. This split enables lazy-loading: the configuration can be loaded immediately, while the resource is loaded later when/if it is required. This can significantly improve startup times, as long as all time-consuming initializations are performed only when loading the resource. This split also helps to avoid cyclical imports between the plugin code and the `nomad-lab` package.

For example the entry point instance for a parser is contained in `.../parsers/__init__.py` and it contains e.g. the name, version and any additional entry point-specific parameters that control its behaviour. The entry point has a `load` method than can be called lazily to return the resource, which is a `Parser` instance defined in `.../parsers/myparser.py`.

In `pyproject.toml` you can expose plugin entry points for automatic discovery. E.g. to expose an app and a package, you would add the following to `pyproject.toml`:

```toml
[project.entry-points.'nomad.plugin']
myapp = "nomad_example.parsers:myapp"
mypackage = "nomad_example.schema_packages:mypackage"
```

Here it is important to use the `nomad.plugin` group name in the `project.entry-points` header. The plugin name used on the left side (`mypackage`) can be arbitrary, what matters is that the key (`"nomad_example.schema_packages:mypackage"`) is a path pointing to a plugin entry point instance inside the python code. This unique key will be used to identify the plugin entry point when e.g. accessing it to read some of it's configuration values.

You can read more about how to write different types of entry points in their dedicated documentation pages or learn more about the [Python entry point mechanism](https://setuptools.pypa.io/en/latest/userguide/entry_point.html).

### Controlling loading of plugin entry points

By default, plugin entry points are automatically loaded, and as an administrator you only need to install the Python package. You can, however, control which entry points to load by explicitly including/excluding them in your `nomad.yaml`. For example, if a plugin has the following `pyproject.toml`:

```toml
[project.entry-points.'nomad.plugin']
myparser = "nomad_example.parsers:myparser"
```

You could disable the parser entry point in your `nomad.yaml` with:

```yaml
plugins:
  entry_points:
    exclude: ["nomad_plugin.parsers:myparser"]
```

### Extending and using the entry point

The plugin entry point is an instance of a [`pydantic`](https://docs.pydantic.dev/1.10/) model. This base model may already contain entry point-specific fields (such as the file extensions that a parser plugin will match) but it is also possible to extend this model to define additional fields that control your plugin behaviour.

To specify new configuration options, you can add new `pydantic` fields to the subclass. For example, if we wanted to add a new configuration option for a parser, we could do the following:

```python
from pydantic import Field
from nomad.config.models.plugins import ParserEntryPoint


class MyParserEntryPoint(ParserEntryPoint):
    parameter: int = Field(0, description='Config parameter for this parser.')
```

where we have defined a new subclass of `ParserEntryPoint` and added a new configuration field `parameter`. The plugin users can then control these settings in their `nomad.yaml` using `plugins.entry_points.options`:

```yaml
plugins:
  entry_points:
    options:
      "nomad_example.parsers:myparser":
        parameter: 47
```

Note that the model will also validate the values coming from `nomad.yaml`, and you should utilize the validation mechanisms of `pydantic` to provide users with helpful messages about invalid configuration.

In your code, you can then access the whole entry point by loading it with `config.get_plugin_entry_point`:

```python
from nomad.config import config

configuration = config.get_plugin_entry_point('nomad_example.parsers:myparser')
print(f'The parser parameter is: {configuration.parameter}')
```

## Plugin development guidelines

### Linting and formatting

While developing NOMAD plugins, we highly recommend using a Python linter, such as [Ruff](https://docs.astral.sh/ruff), to analyze and enforce coding standards in your plugin projects. This also ensures smoother integration and collaboration. If you have used our [template repository](https://github.com/FAIRmat-NFDI/nomad-plugin-template), you will automatically have `ruff` defined as a development dependency with suitable defaults set in `pyproject.toml` together with a GitHub actions that runs the linting and formatting checks on each push to the Git repository.

### Testing

For testing, you should use [pytest](https://docs.pytest.org/), and a folder structure that mimics the package layout with test modules named after the tested module. For example, if you are developing a parser in `myparser.py`, the test folder structure should look like this:

```txt
├── nomad-example-plugin
│   ├── src
|   │   ├── nomad_example
|   |   │   ├── parsers
|   |   │   │   ├── myparser.py
|   |   │   │   ├── __init__.py
│   ├── tests
|   │   ├── parsers
|   |   │   ├── test_myparser.py
|   |   │   ├── conftest.py
|   │   ├── conftest.py
```

Any shared test utilities (such as `pytest` fixtures) should live in `conftest.py` modules placed at the appropriate level in the folder hierarchy, i.e. utilities dealing with parsers would live in `tests/parsers/conftest.py`, while root level utilities would live in `tests/conftest.py`. If you have used our [template repository](https://github.com/FAIRmat-NFDI/nomad-plugin-template), you will automatically have an initial test folder structure, `pytest` defined as a development dependency in `pyproject.toml` and a GitHub action that runs the test suite on each push to the Git repository.

In the `pytest` framework, test cases are created by defining functions with the `test_` prefix, which perform assertions. A typical test case could look like this:

```python
def test_parse_file():
    parser = MyParser()
    archive = EntryArchive()
    parser.parse('tests/data/example.out', archive, logging)

    sim = archive.data
    assert len(sim.model) == 2
    assert len(sim.output) == 2
    assert archive.workflow2.x_example_magic_value == 42
```

You can run all the tests in the `tests/` directory with:

```shell
python -m pytest -svx tests
```

### Documentation

As your plugin matures, you should also think about documenting its usage. We recommend using [`mkdocs`](https://www.mkdocs.org/) to create your documentation as a set of markdown files. If you have used our [template repository](https://github.com/FAIRmat-NFDI/nomad-plugin-template), you will automatically have an initial documentation folder structure, `mkdocs` defined as a development dependency in `pyproject.toml` and a GitHub action that builds the docs to a separate `gh-pages` branch each push to the Git repository. Note that if you wish to host the documentation using [GitHub pages](https://pages.github.com/), you need to [enable](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site#publishing-from-a-branch) this in the repository settings.

## Publishing a plugin

!!! warning "Attention"
    The standard processes for publishing plugins and using plugins from other developers are still being worked out. The "best" practices mentioned in the following are preliminary. We aim to set up a dedicated plugin registry that allows you to publish your plugin and find plugins from others.

### GitHub repository

The simplest way to publish a plugin is to have it live in a publicly shared Git
repository. The package can then be installed with:

```sh
pip install git+https://<repository_url>
```

!!! note
    If you develop a plugin in the context of [FAIRmat](https://github.com/fairmat-nfdi) or the [NOMAD CoE](https://github.com/nomad-coe), put your plugin repositories in the corresponding GitHub organization.

### PyPI/pip package

You may additionally publish the plugin package in PyPI. Learn from the PyPI documentation how to [create a package for PyPI](https://packaging.python.org/en/latest/tutorials/packaging-projects/){:target="_blank"}. We recommend to use the `pyproject.toml`-based approach.

The PyPI documentation provides further information about how to [publish a package to PyPI](https://packaging.python.org/en/latest/tutorials/packaging-projects/#uploading-the-distribution-archives){:target="_blank"}. If you have access to the MPCDF GitLab and NOMAD's presence there, you can also
use the `nomad-FAIR` package registry:

```
pip install twine
twine upload \
    -u <username> -p <password> \
    --repository-url https://gitlab.mpcdf.mpg.de/api/v4/projects/2187/packages/pypi \
    dist/nomad-example-plugin-*.tar.gz
```
