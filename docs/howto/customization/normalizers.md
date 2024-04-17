# How to write a normalizer

## Introduction

A normalizer can be any Python algorithm that takes the archive of an entry as input
and manipulates (usually expands) the given archive. This way, a normalizer can add
additional sections and quantities based on the information already available in the
archive.

All normalizers are executed after parsing. Normalizers are run for each entry (i.e. each
of the parsed mainfiles). Normalizers are run in a particular order, and
you can make assumptions about the availability of data created by other normalizers.
A normalizer is run in any case, but it might choose not to do anything. A normalizer
can perform any operation on the archive, but in general it should only add more
information, not alter existing information.

## Getting started
A normalizer can be developed either as part of the NOMAD package or as a plugin. In the
case of the former, the code can be added to the `nomad.normalizing` module. Developing
it as a plugin is discussed [How to develop a normalizer plugin](#developing-a-normalizer-plugin). In the
following, we illustrate the structure of a normalizer and describe how to interface it
to NOMAD in each of the two cases.

## Starting example

This is an example for a very simple normalizer that computes the unit cell volume from
a given lattice and adds it to each of the member of the system in the `run` section of the archive.

```python
from nomad.normalizing import Normalizer
from nomad.atomutils import get_volume

class UnitCellVolumeNormalizer(Normalizer):

    normalizer_level = 1

    def normalize(self):
        for system in self.entry_archive.run[-1].system:
            system.unit_cell_volume = get_volume(system.lattice_vectors.magnitude)

            self.logger.debug('computed unit cell volume', system_index=system.m_parent_index)
```

A normalizer simply inherits `Normalizer` and implements the `normalize` method. The
`entry_archive` is available as a field. There is also a logger on the object that can be used.
Be aware that the processing will already report the run of the normalizer, log its
execution time and any exceptions that might been thrown.

## Implementing a system normalizer

There is a special base class for normalizing systems that allows to run the normalization
on all (or only the resulting) `representative` systems:

```python
from nomad.normalizing import SystemBasedNormalizer
from nomad.atomutils import get_volume

class UnitCellVolumeNormalizer(SystemBasedNormalizer):
    def _normalize_system(self, system, is_representative):
        system.unit_cell_volume = get_volume(system.lattice_vectors.magnitude)
```

For `SystemBasedNormalizer`, we implement the `_normalize_system` method.
The parameter `is_representative` will be true for the `representative` systems. The
representative system refers to the system that corresponds to the calculation result.
It is determined by scanning the archive sections starting with `workflow2` until
the system fitting the criterion is found. For example, it refers to the final step in a
geometry optimization or other workflow.

Of course, if you add new information to the archive, this also needs to be defined in the
schema (see [How-to extend the schema](schemas.md#extending-existing-sections)). For example you could extend the section system with a special system definition
that extends the existing section system definition:

```python
import numpy as np
from nomad.datamodel.metainfo import runschema
from nomad.metainfo import Section, Quantity

class UnitCellVolumeSystem(runschema.system.System):
    m_def = Section(extends_base_section=True)
    unit_cell_volume = Quantity(np.dtype(np.float64), unit='m^3')
```

Here, we used the schema definition for the `run` section defined in this [plugin](schemas.md#pre-defined-schemas-in-nomad).

## Adding a normalizer to the processing

For any new normalizer class to be recognized by the processing, the normalizer class
needs to be added to the list of normalizers in the config file `nomad/config/__init__.py`.
For [normalizer plugins](#developing-a-normalizer-plugin) one needs to include it in
the `nomad.yaml` file (see [Adding a plugin to NOMAD](plugins.md#add-a-plugin-to-your-nomad)).
By default, the execution order of the normalizers during processing is determined by the
order of the normalizers in the list. One can specify the order of a normalizer relative
to the other normalizers by specifying the `normalizer_level` field. The following lists
the normalizers used in NOMAD in order of execution:

| Normalizer class             | Path/Project url                                                            |
| ---------------------------- | --------------------------------------------------------------------------- |
| SimulationWorkflowNormalizer | <https://github.com/nomad-coe/nomad-schema-plugin-simulation-workflow.git>  |
| SystemNormalizer             | <https://github.com/nomad-coe/nomad-normalizer-plugin-system.git>           |
| SoapNormalizer               | <https://github.com/nomad-coe/nomad-normalizer-plugin-soap.git>             |
| SpectraNormalizer            | <https://github.com/nomad-coe/nomad-normalizer-plugin-spectra.git>          |
| OptimadeNormaizer            | nomad.normalizing.optimade.OptimadeNormalizer                               |
| MetainfoNormaizer            | nomad.normalizing.optimade.MetainfoNormalizer                               |
| DosNormalizer                | <https://github.com/nomad-coe/nomad-normalizer-plugin-dos.git>              |
| BandStructureNormalizer      | <https://github.com/nomad-coe/nomad-normalizer-plugin-bandstructure.git>    |
| ResultsNormalizer            | nomad.normalizing.results.ResultsNormalizer                                 |

In the future, all normalizers will be migrated to plugins therefore new normalizers should be developed as plugins.

## Testing a normalizer

To simply try out a normalizer, you could use the CLI and run the parse command:

```shell
nomad --debug parse --show-archive <path-to-example-file>
```

But eventually you need to add a more formal test. Place your `pytest` tests in
`tests/normalizing/test_unitcellvolume.py` similar to the existing tests or in the plugin
tests for the case of plugins. Necessary test data can be added to `tests/data/normalizers`.

## Developing a Normalizer plugin
Fork and clone the [normalizer example project](https://github.com/nomad-coe/nomad-normalizer-plugin-example)
as described in [How to develop and publish plugins](plugins.md). The normalizer class is
defined in `nomadnormalizerexample/normalizer.py`

```python
from nomad.normalizing import Normalizer
from nomad.datamodel.metainfo.workflow import Workflow


class ExampleNormalizer(Normalizer):

    domain = None

    def normalize(self, logger):
        super().normalize(logger)
        logger.info('ExampleNormalizer called')

        self.entry_archive.workflow2 = Workflow(name='Example workflow')
```

In this simple example, we create a the workflow2 section in the archive.

## Normalizer plugin metadata

{{pydantic_model('nomad.config.models.plugins.Normalizer')}}