# Plugins

Plugins allow one to add Python-based functionality to NOMAD without a custom NOMAD image or release. Plugins can be installed at NOMAD start-up time. Therefore, a NOMAD installation or [Oasis](../howto/oasis/install.md) can be configured with a different custom set of plugins or disable unnecessary plugins.

!!! note

    You might also want to read [the how-to guide on plugins](../howto/plugins/plugins.md)

## Plugin entry point reference

This is a list of the available plugin entry point configuration models.

{{ pydantic_model('nomad.config.models.plugins.AppEntryPoint') }}
{{ pydantic_model('nomad.config.models.plugins.NormalizerEntryPoint') }}
{{ pydantic_model('nomad.config.models.plugins.ParserEntryPoint') }}
{{ pydantic_model('nomad.config.models.plugins.SchemaPackageEntryPoint') }}

## Default plugin entry points

This is a list of the plugin entry points that are activated by default:

{{ plugin_entry_point_list() }}
