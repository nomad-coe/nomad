# Archive and Metainfo

## Introduction

NOMAD stores all processed data in a *well defined*, *structured*, and *machine readable*
format. Well defined means that each element is supported by a formal definition that provides
a name, description, location, shape, type, and possible unit for that data. It has a
hierarchical structure that logically organizes data in sections and subsections and allows
cross-references between pieces of data. Formal definitions and corresponding
data structures enable the machine processing of NOMAD data.

![archive example](assets/archive-example.png)
#### The Metainfo is the schema for Archive data.
The Archive stores descriptive and structured information about materials-science
data. Each entry in NOMAD is associated with one Archive that contains all the processed
information of that entry. What information can possibly exist in an archive, how this
information is structured, and how this information is to be interpreted is governed
by the Metainfo.

#### On schemas and definitions
Each piece of Archive data has a formal definition in the Metainfo. These definitions
provide data types with names, descriptions, categories, and further information that
applies to all incarnations of a certain data type.

Consider a simulation `Run`. Each
simulation run in NOMAD is characterized by a *section*, that is called *run*. It can contain
*calculation* results, simulated *systems*, applied *methods*, the used *program*, etc.
What constitutes a simulation run is *defined* in the metainfo with a *section definition*.
All other elements in the Archive (e.g. *calculation*, *system*, ...) have similar definitions.

Definitions follow a formal model. Depending on the definition type, each definition
has to provide certain information: *name*, *description*, *shape*, *units*, *type*, etc.

#### Types of definitions

- *Sections* are the building block for hierarchical data. A section can contain other
  sections (via *subsections*) and data (via *quantities*).
- *Subsections* define a containment relationship between sections.
- *Quantities* define a piece of data in a section.
- *References* are special quantities that allow to define references from a section to
  another section or quantity.
- *Categories* allow to categorize definitions.
- *Packages* are used to organize definitions.

#### Interfaces
The Archive format and Metainfo schema is abstract and not not bound to any
specific storage format. Archive and Metainfo can be represented in various ways.
For example, NOMAD internally stores archives in a binary format, but serves them via
API in json. Users can upload archive files (as `.archive.json` or `.archive.yaml`) files.
Metainfo schema can be programmed with Python classes, but can also be uploaded as
archive files (the Metainfo itself is just a specific Archive schema). The following
chart provides a sense of various ways that data can be entered into NOMAD:

![nomad data flow](assets/data-flow.png)

There are various interface to provide or retrieve Archive data and Metainfo schemas.
The following documentation sections will explain a few of them.

## Archive JSON API

The [API section](api.md#access-archives) demonstrates how to access an Archive. The
API will give you JSON data likes this:

```json title="https://nomad-lab.eu/prod/v1/api/v1/entries/--dLZstNvL_x05wDg2djQmlU_oKn/archive"
{
    "run": [
        {
            "program": {...},
            "method": [...],
            "system": [
                {...},
                {...},
                {...},
                {...},
                {
                    "type": "bulk",
                    "configuration_raw_gid": "-ZnDK8gT9P3_xtArfKlCrDOt9gba",
                    "is_representative": true,
                    "chemical_composition": "KKKGaGaGaGaGaGaGaGaGa",
                    "chemical_composition_hill": "Ga9K3",
                    "chemical_composition_reduced": "K3Ga9",
                    "atoms": {...},
                    "springer_material": [...],
                    "symmetry": [...]
                }
            ]
            "calculation": [...],
        }
    ],
    "workflow": [...],
    "metadata": {...},
    "results":{
        "material": {...},
        "method": {...},
        "properties": {...},
    }
}
```

This will show you the Archive as a hierarchy of JSON objects (each object is a section),
where each key is a property (e.g. a quantity or subsection). Of course you can use
this data in this JSON form. You can expect that the same keys (each item has a formal
definition) always provides the same type of data. However, not all keys are present in
every archive, and not all lists might have the same number of objects. This depends on the
data. For example, some *runs* contain many systems (e.g. geometry optimizations), others
don't; typically *bulk* systems will have *symmetry* data, non bulk systems might not.
To learn what each key means, you need to look up its definition in the Metainfo.

{{ metainfo_data() }}


## Archive Python interface

In Python, JSON data is typically represented as nested combinations of dictionaries
and lists. Of course, you could work with this right away. To make it easier for Python
programmers the [NOMAD Python package](pythonlib.md) allows you to use this
JSON data with a higher level interface, which provides the following advantages:

- code completion in dynamic coding environments like Jupyter notebooks
- a cleaner syntax that uses attributes instead of dictionary access
- all higher dimensional numerical data is represented as numpy arrays
- allows to navigate through references
- numerical data has a Pint unit attached to it

For each section the Python package contains a Python class that corresponds to its
definition in the metainfo. You can use these classes to access `json_data` downloaded
via API:
```python
from nomad.datamodel import EntryArchive

archive = EntryArchive.m_from_dict(json_data)
calc = archive.run[0].calculation[-1]
total_energy_in_ev = calc.energy.total.value.to(units.eV).m
formula = calc.system_ref.chemical_formula_reduced
```

Archive data can also be serialized into JSON again:
```python
import json

print(json.dumps(calc.m_to_dict(), indent=2))
```

## Metainfo (and Archive) .yaml interface

Calling it `.yaml` interface is a bit wrong. Yaml and Json are simple formats for
structured data. In many respects they are the same. After these files got loaded
(e.g into Python dictionaries and array or Javascript objects) all differences disappear.
Therefore, you can apply everything about the `archive.yaml` "format" also to `archive.json`.

Above you saw how the JSON API is representing Archive data in JSON. The same can be used
to upload data to NOMAD. NOMAD will interpret all files called `*.archive.yaml` or
`*.archive.json` accordingly. This means the data in those wills has to match the existing
NOMAD schemas.

Here is a basic example for a `archive.yaml` file that contains data and a matching schema:

```yaml
--8<-- "examples/data/custom-schema/intra-entry.archive.yaml"
```

## Metainfo Python interface

To learn more about the Python interface, look at the [Metainfo documentation](metainfo.md)
that explains how the underlying Python classes work, and how you can extend the
metainfo by providing your own classes.


## Custom metainfo schemas (e.g. for ELNs)

With custom metainfo schemas, you can extend the existing NOMAD schemas. This will allow
you to provide .json or .yaml data even if NOMAD does not support it out of the box. It
also allows you to define editable data (e.g. to create an ELN).


The idea is that a schema defines all possible data structures. NOMAD uses this information
to parse data or provide editing functions with schema specific forms.
This is what a custom schema can look like:

```yaml
--8<-- "examples/data/custom-schema/simple-schema.archive.yaml"
```

### Annotations (for ELNs)

Schema elements can have annotations. These annotations provide additional information
that NOMAD can use to alter its behavior around these definitions. A reference for
these annotations will be added here soon.

Many annotations control the representation of data in the GUI. This can be for plots
or data entry/editing capabilities. As part of the GUI, you'll find an overview about all
ELN edit annotations and components [here]({{ nomad_url() }}/../gui/dev/editquantity).

The ELN example upload, that you can create from the upload page, contains a schema that
uses and demonstrates all our current annotations:

```yaml
--8<-- "examples/data/eln/schema.archive.yaml"
```

### Custom normalizers

For custom schemas, you might want to add custom normalizers. All files are parsed
and normalized when they are uploaded or changed. The NOMAD metainfo Python interface
allows you to add functions that are called when your data is normalized.

Here is an example:

```python
--8<-- "examples/archive/custom_schema.py"
```

To add a `normalize` function, your section has to inherit from `ArchiveSection` which
provides the base for this functionality. Now you can overwrite the `normalize` function
and add you own behavior. Make sure to call the `super` implementation properly to
support schemas with multiple inheritance.

If we parse an archive like this:

```yaml
--8<-- "examples/archive/custom_data.archive.yaml"
```

we will get a final normalized archive that contains our data like this:

```json
{
  "data": {
    "m_def": "examples.archive.custom_schema.SampleDatabase",
    "samples": [
      {
        "added_date": "2022-06-18T00:00:00+00:00",
        "formula": "NaCl",
        "sample_id": "2022-06-18 00:00:00+00:00--NaCl"
      }
    ]
  }
}
```