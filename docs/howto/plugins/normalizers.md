# How to write a normalizer

A normalizer takes the archive of an entry as input and manipulates (usually expands) the given archive. This way, a normalizer can add additional sections and quantities based on the information already available in the archive. All normalizers are executed in the order [determined by their `level`](#control-normalizer-execution-order) after parsing, but the normalizer may decide to not do anything based on the entry contents.

This documentation shows you how to write a plugin entry point for a normaliser. You should read the [documentation on getting started with plugins](./plugins.md) to have a basic understanding of how plugins and plugin entry points work in the NOMAD ecosystem.

## Getting started

You can use our [template repository](https://github.com/FAIRmat-NFDI/nomad-plugin-template) to create an initial structure for a plugin containing a normalizer. The relevant part of the repository layout will look something like this:

```txt
nomad-example
   ├── src
   │   ├── nomad_example
   │   │   ├── normalizers
   │   │   │   ├── __init__.py
   │   │   │   ├── mynormalizer.py
   ├── LICENSE.txt
   ├── README.md
   └── pyproject.toml
```

See the documentation on [plugin development guidelines](./plugins.md#plugin-development-guidelines) for more details on the best development practices for plugins, including linting, testing and documenting.

## Normalizer entry point

The entry point defines basic information about your normalizer and is used to automatically load the normalizer code into a NOMAD distribution. It is an instance of a `NormalizerEntryPoint` or its subclass and it contains a `load` method which returns a `nomad.normalizing.Normalizer` instance that will perform the actual normalization. You will learn more about the `Normalizer` class in the next sections. The entry point should be defined in `*/normalizers/__init__.py` like this:

```python
from pydantic import Field
from nomad.config.models.plugins import NormalizerEntryPoint


class MyNormalizerEntryPoint(NormalizerEntryPoint):

    def load(self):
        from nomad_example.normalizers.mynormalizer import MyNormalizer

        return MyNormalizer(**self.dict())


mynormalizer = MyNormalizerEntryPoint(
    name = 'MyNormalizer',
    description = 'My custom normalizer.',
)
```

Here you can see that a new subclass of `NormalizerEntryPoint` was defined. In this new class you can override the `load` method to determine how the `Normalizer` class is instantiated, but you can also extend the `NormalizerEntryPoint` model to add new configurable parameters for this normalizer as explained [here](./plugins.md#extending-and-using-the-entry-point).

We also instantiate an object `mynormalizer` from the new subclass. This is the final entry point instance in which you specify the default parameterization and other details about the normalizer. In the reference you can see all of the available [configuration options for a `NormalizerEntryPoint`](../../reference/plugins.md#normalizerentrypoint).


The entry point instance should then be added to the `[project.entry-points.'nomad.plugin']` table in `pyproject.toml` in order for the normalizer to be automatically detected:

```toml
[project.entry-points.'nomad.plugin']
mynormalizer = "nomad_example.normalizers:mynormalizer"
```

## `Normalizer` class

The resource returned by a normalizer entry point must be an instance of a `nomad.normalizing.Normalizer` class. This normalizer definition should be contained in a separate file (e.g. `*/normalizer/mynormalizer.py`) and could look like this:

```python
from typing import Dict

from nomad.datamodel import EntryArchive
from nomad.normalizing import Normalizer


class MyNormalizer(Normalizer):
    def normalize(
        self,
        archive: EntryArchive,
        logger=None,
    ) -> None:
        logger.info('MyNormalizer called')
```

The minimal requirement is that your class has a `normalize` function, which as input takes:

 - `archive`: The [`EntryArchive` object](../../reference/glossary.md#archive) in which the normalization results will be stored
 - `logger`: Logger that you can use to log normalization events into

## `SystemBasedNormalizer` class

`SystemBasedNormalizer` is a special base class for normalizing systems that allows to run the normalization on all (or only the resulting) `representative` systems:

```python
from nomad.normalizing import SystemBasedNormalizer
from nomad.atomutils import get_volume

class UnitCellVolumeNormalizer(SystemBasedNormalizer):
    def _normalize_system(self, system, is_representative):
        system.unit_cell_volume = get_volume(system.lattice_vectors.magnitude)
```

For `SystemBasedNormalizer`, we implement the `_normalize_system` method. The parameter `is_representative` will be true for the `representative` systems. The representative system refers to the system that corresponds to the calculation result. It is determined by scanning the archive sections starting with `workflow2` until the system fitting the criterion is found. For example, it refers to the final step in a geometry optimization or other workflow.

Of course, if you add new information to the archive, this also needs to be defined in the schema (see [How-to extend the schema](schema_packages.md#extending-existing-sections)). For example you could extend the section system with a special system definition that extends the existing section system definition:

```python
import numpy as np
from nomad.datamodel.metainfo import runschema
from nomad.metainfo import Section, Quantity

class UnitCellVolumeSystem(runschema.system.System):
    m_def = Section(extends_base_section=True)
    unit_cell_volume = Quantity(np.dtype(np.float64), unit='m^3')
```

Here, we used the schema definition for the `run` section defined in this [plugin](schema_packages.md#schema-packages-developed-by-fairmat).

## Control normalizer execution order

`NormalizerEntryPoints` have an attribute `level`, which you can use to control their execution order. Normalizers are executed in order from lowest level to highest level. The default level for normalizers is `0`, but this can be changed per installation using `nomad.yaml`:

```yaml
plugins:
  entry_points:
    options:
      "nomad_example.normalizers:mynormalizer1":
        level: 1
      "nomad_example.normalizers:mynormalizer2":
        level: 2
```

## Running the normalizer

If you have the plugin package and `nomad-lab` installed in your Python environment, you can run the normalization as a part of the parsing process using the NOMAD CLI:

```shell
nomad parse <filepath> --show-archive
```

The output will return the final archive in JSON format.

Normalization can also be run within a python script (or Jupyter notebook), e.g., to facilate debugging, with the following code:

```python
from nomad.datamodel import EntryArchive
from nomad_example.normalizers.mynormalizer import MyNormalizer
import logging

p = MyNormalizer()
a = EntryArchive()
p.normalize(a, logger=logging.getLogger())

print(a.m_to_dict())
```

## Normalizers developed by FAIRmat

The following is a list of plugins containing normalizers developed by FAIRmat:

| Normalizer class             | Path/Project url                                                            |
| ---------------------------- | --------------------------------------------------------------------------- |
| SimulationWorkflowNormalizer | <https://github.com/nomad-coe/nomad-schema-plugin-simulation-workflow.git>  |
| SystemNormalizer             | <https://github.com/nomad-coe/nomad-normalizer-plugin-system.git>           |
| SoapNormalizer               | <https://github.com/nomad-coe/nomad-normalizer-plugin-soap.git>             |
| SpectraNormalizer            | <https://github.com/nomad-coe/nomad-normalizer-plugin-spectra.git>          |
| DosNormalizer                | <https://github.com/nomad-coe/nomad-normalizer-plugin-dos.git>              |
| BandStructureNormalizer      | <https://github.com/nomad-coe/nomad-normalizer-plugin-bandstructure.git>    |

To refine an existing normalizer, you should install it via the `nomad-lab` package:

```shell
pip install nomad-lab
```

Clone the normalizer project:

```shell
git clone <normalizer-project-url>
cd <normalizer-dir>
```

Either remove the installed normalizer and `pip install` the cloned version:

```shell
rm -rf <path-to-your-python-env>/lib/python3.9/site-packages/<normalizer-module-name>
pip install -e .
```

Or set `PYTHONPATH` so that the cloned code takes precedence over the installed code:

```shell
PYTHONPATH=. nomad parse <path-to-example-file>
```

Alternatively, you can also do a full [developer setup](../develop/setup.md) of the NOMAD infrastructure and
enhance the normalizer there.