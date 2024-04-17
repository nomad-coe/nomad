# Plugins

Plugins allow one to add Python-based functionality to NOMAD without a custom build
NOMAD image or release. Plugins can be installed at NOMAD start-up time. Therefore, a NOMAD
installation or [Oasis](../howto/oasis/install.md) can be configured with a different
custom set of plugins or disable unnecessary plugins.

NOMAD support different kinds of plugins:

- Python **schema**
- **parser**
- **normalizer**
- additional custom **APIs** (coming soon...)

!!! note

    You might also want to read [the plugin how-tos](../howto/customization/plugins.md)

## Types of plugins

We provide template projects on GitHub for each kind plugin in NOMAD.

- [schema plugin](https://github.com/nomad-coe/nomad-schema-plugin-example){:target="_blank"}
- [parser plugin](https://github.com/nomad-coe/nomad-parser-plugin-example){:target="_blank"}
- [normalizer plugin](https://github.com/nomad-coe/nomad-normalizer-plugin-example.git){:target="_blank"}

You can fork these projects and follow the instructions in their `README.md`. These
instructions will give you everything you need to run and test your plugin.

## Built-in plugins
This is a list of all built-in plugins:

{{ plugin_list() }}
