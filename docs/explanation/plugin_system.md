# NOMAD plugin system

## Introduction

NOMAD is used by many research communities with their specific data, workflows, and analysis tools. NOMAD plugins are key
to adopt NOMAD to these heterogeneous environments.
You can think of plugins as “add-ons” that provide additional capabilities.
Each plugin is a small independant software project that integrates with the core NOMAD and provides features without modifictions to the core NOMAD itself.
Some key advantages of using plugins:

- **Modularity**: You can pick and choose which features or functions to add, rather than having everything baked into the core NOMAD.

- **Customizability**: Users can add their own plugins to address specific use cases, without changing the official NOMAD software.

- **Easy updates**: If a feature needs to be updated or improved, it can be done at the plugin level, without having to release a new NOMAD version.

- **Collaboration**: Since plugins are independent, multiple developers can work on different features in parallel and with different release cycles without interfering with each other.

## Architecture

There are three core components to the plugin system:

- **Distributions** define lists of plugins and their version. A distribution is a small
  Git and Python project that maintains a list of plugin dependencies in its `pyproject.toml`.
- **Plugins** are Git and Python projects that contain one or many _entry points_.
  We provide a [template repository](https://github.com/FAIRmat-NFDI/nomad-plugin-template)
  for a quick start into plugin development.
- **Entry points** are individual contributions (e.g. parsers, schemas, or apps)
  which are defined using a feature of Python called [_entry points_](https://setuptools.pypa.io/en/latest/userguide/entry_point.html).

<figure markdown style="width: 100%">
  ``` mermaid
  %%{init:{'flowchart':{'nodeSpacing': 25, 'subGraphTitleMargin': {'top': 5, 'bottom': 10}, 'padding': 10}}}%%
  graph LR
    subgraph NOMAD Distribution
      subgraph NOMAD Plugin C
        ro7(Entry point: Schema 1)
        ro8(Entry point: Schema 2)
      end
      subgraph NOMAD Plugin B
        ro4(Entry point: Schema)
        ro5(Entry point: App 1)
        ro6(Entry point: App 2)
      end
      subgraph NOMAD Plugin A
        ro1(Entry point: Schema)
        ro2(Entry point: Parser 1)
        ro3(Entry point: Parser 2)
      end
    end
  ```
  <figcaption>Relation between NOMAD distributions, plugins and entry points.</figcaption>
</figure>

This architecture allows plugin developers to freely choose a suitable granularity for their use case: they may create a single plugin package that contains everything that e.g. a certain lab needs: schemas, parsers and apps. Alternatively they may also develop multiple plugins, each containing a single entry point. Reduction in package scope can help in developing different parts indepedently and also allows plugin users to choose only the parts that they need.

## Entry point discovery

Entry points are like pre-defined connectors or hooks that allow the main system to recognize and load the code from plugins without needing to hard-code them directly into the platform. This mechanism enables the automatic discovery of plugin code. Plugin entry points represent different types of customizations that can be added to a NOMAD installation. The following plugin entry point types are currently supported:

- [Apps](../howto/plugins/apps.md)
- [Example uploads](../howto/plugins/example_uploads.md)
- [Normalizers](../howto/plugins/parsers.md)
- [Parsers](../howto/plugins/parsers.md)
- [Schema packages](../howto/plugins/schema_packages.md)

## Loading plugins

Entry points contain **configuration** (the entry point), but also a separate **resource** (the implementation). This split enables lazy-loading: the configuration can be loaded immediately, while the resource is loaded later when/if it is required. This can significantly improve startup times, as long as all time-consuming initializations are performed only when loading the resource. This split also helps to avoid cyclical imports between the plugin code and the `nomad-lab` package. The following diagram illustrates how NOMAD interacts with the entry points in a plugin:

<figure markdown style="width: 100%">
  ``` mermaid
  %%{init:{'sequence':{'mirrorActors': false}}}%%
  sequenceDiagram
    autonumber
    rect
      NOMAD->>Plugin: Request entry points matching the nomad.plugin group
      Plugin-->>NOMAD: Return all entry point configurations
    end
      Note over NOMAD: Other tasks
    rect
      NOMAD->>Plugin: Request a specific entry point resource
      opt
        Plugin->>NOMAD: Request configuration overrides from nomad.yaml
        NOMAD-->>Plugin: Return final configuration for entry point
      end
      Plugin-->>NOMAD: Return fully initialized entry point resource
    end
  ```
  <figcaption>NOMAD interaction with a plugin.</figcaption>
</figure>

1.  When NOMAD starts, it scans for plugin entry points defined under the `nomad.plugin` group in all of the Python packages that have been installed.
2.  The plugin returns all entry points that it has registered in `pyproject.toml` under the `nomad.plugin` group. This only loads the configuration, but does not yet load the resource, i.e. main Python implementation.
3.  When NOMAD needs to load the actual resource for an entry point (e.g. a parser), loads it by using the configuration instance.
4.  When the resource is being loaded, the entry point may ask for any configuration overrides that may have been set in `nomad.yaml`.
5.  NOMAD will return the final validated configuration that contains the default values and possible overrides.
6.  The plugin loads and returns the resource using the final configuration. This typically involves creating an instance of a specific class, e.g. `Parser` in the case of parser entry points.

## Learn how to write plugins

You can learn more about plugin development in the [Get started with plugins](../howto/plugins/plugins.md) -page.
