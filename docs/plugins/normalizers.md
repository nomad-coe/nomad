A normalizer can be any Python algorithm that takes the archive of an entry as input
and manipulates (usually expands) the given archive. This way, a normalizer can add
additional sections and quantities based on the information already available in the
archive.

All normalizers are executed after parsing. Normalizers are run for each entry (i.e. each
set of files that represent a code run). Normalizers are run in a particular order, and
you can make assumptions about the availability of data created by other normalizers.
A normalizer is run in any case, but it might choose not to do anything. A normalizer
can perform any operation on the archive, but in general it should only add more
information, not alter existing information.

## Getting started

Fork and clone the [normalizer example project](https://github.com/nomad-coe/nomad-normalizer-plugin-example) as described in [before](plugins.md). Follow the original [how-to on writing a parser](../develop/normalizers.md).

{{pydantic_model('nomad.config.plugins.Normalizer', heading='### Normalizer plugin metadata')}}