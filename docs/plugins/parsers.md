NOMAD uses parsers to convert raw code input and output files into NOMAD's common Archive format. This is the documentation on how to develop such a parser.

## Getting started

Fork and clone the [parser example project](https://github.com/nomad-coe/nomad-parser-plugin-example) as described in [before](plugins.md). Follow the original [How-to write a parser](../develop/parsers.md) documentation.

{{pydantic_model('nomad.config.plugins.Parser', heading='### Parser plugin metadata', hide=['code_name','code_category','code_homepage','metadata'])}}