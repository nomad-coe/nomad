# How to write a schema plugin

## Getting started

Fork and clone the [schema example project](https://github.com/nomad-coe/nomad-schema-plugin-example){:target="_blank"} as described in [How-to mount a plugin](../oasis/plugins_install.md).

## Writing schemas in Python compared to YAML schemas

In this [guide](basics.md), we explain how to write and upload schemas in the `.archive.yaml` format. Writing and uploading such YAML schemas is a good way for NOMAD users to add schemas. But it has limitations. As a NOMAD developer or Oasis administrator you can add Python schemas to NOMAD. All built-in NOMAD schemas (e.g. for electronic structure code data) are written in Python and are part of the NOMAD sources (`nomad.datamodel.metainfo.*`).

There is a 1-1 translation between Python schemas (written in classes) and YAML (or JSON) schemas (written in objects). Both use the same fundamental concepts, like *section*, *quantity*, or *subsection*, introduced in [YAML schemas](basics.md).

## Starting example

```python
from nomad.metainfo import MSection, Quantity, SubSection, Units

class System(MSection):
    '''
    A system section includes all quantities that describe a single simulated
    system (a.k.a. geometry).
    '''

    n_atoms = Quantity(
        type=int, description='''
        A Defines the number of atoms in the system.
        ''')

    atom_labels = Quantity(
        type=MEnum(ase.data.chemical_symbols), shape['n_atoms'])
    atom_positions = Quantity(type=float, shape=['n_atoms', 3], unit=Units.m)
    simulation_cell = Quantity(type=float, shape=[3, 3], unit=Units.m)
    pbc = Quantity(type=bool, shape=[3])

class Run(MSection):
    section_system = SubSection(sub_section=System, repeats=True)
```

We define a simple metainfo schema with two *sections* called `System` and `Run`.
Each section can have two types of properties: *quantities* and *subsections*. Sections and their properties are defined with
Python classes and their attributes. Each *quantity* defines a piece of data. Basic quantity attributes are `type`, `shape`,
`unit`, and `description`.

*Subsections* allow the placement of sections within each other, forming containment
hierarchies. Basic subsection attributes are
`sub_section`&mdash;a reference to the section definition of the subsection&mdash;and `repeats`&mdash;determines whether a subsection can be included once or multiple times.

To use the above-defined schema and create actual data, we have to
instantiate the classes:

```python
run = Run()
system = run.m_create(System)
system.n_atoms = 3
system.atom_labels = ['H', 'H', 'O']

print(system.atom_labels)
print(n_atoms = 3)
```

Section *instances* can be used like regular Python objects: quantities and subsections
can be set and accessed like any other Python attribute. Special metainfo methods, starting
with `m_` allow us to realize more complex semantics. For example `m_create` will
instantiate a subsection and add it to the *parent* section in one step.
<!-- ? m_create is deprecated? -->

Another example for an `m_`-method is:

```python
run.m_to_json(indent=2)
```

This will convert the data into JSON:

```json
{
    "m_def" = "Run",
    "systems": [
        {
            "n_atoms" = 3,
            "atom_labels" = [
                "H",
                "H",
                "O"
            ]
        }
    ]
}
```

## Definitions

The following describes in detail the schema language for the NOMAD Metainfo and how it is expressed in Python.


### Common attributes of Metainfo Definitions

In the example, you have already seen the basic Python interface to the Metainfo. *Sections* are
represented in Python as objects. To define a section, you write a Python class that inherits
from `MSection`. To define subsections and quantities you use Python properties. The
definitions themselves are also objects derived from classes. For subsections and
quantities, you directly instantiate `:class:SubSection` and `:class:Quantity`. For sections
there is a generated object derived from `:class:Section` and available via
`m_def` from each *section class* and *section instance*.
<!-- TODO Either fix all cross references with :: syntax here and throughout or remove them -->

These Python classes, used to represent metainfo definitions, form an inheritance
hierarchy to share common properties

- `name`: each definition has a name. This is typically defined by the corresponding
Python property. For example, a section class name becomes the section name; a quantity gets the name
from the variable name used in its Python definition, etc.
- `description`: each definition should have one. Either set it directly or use *doc strings*
- `links`: a list of useful internet references.
- `more`: a dictionary of custom information. Any additional `kwargs` set when creating a definition
    are added to `more`.

### Sections

Sections are defined with Python classes that extend `MSection` (or other section classes).

- `base_sections`: automatically taken from the base classes of the Python class.
- `extends_base_section`: a boolean that determines the inheritance. If this is `False`,
normal Python inheritance implies and this section will inherit all properties (subsections,
quantities) from all base classes. If `True`, all definitions in this section
will be added to the properties of the base class section. This allows the extension of existing
sections with additional properties.

### Quantities

Quantity definitions are the main building block of metainfo schemas. Each quantity
represents a single piece of data. Quantities can be defined with the following attributes:

- `type`: can be a primitive Python type (`str`, `int`, `bool`), a numpy
data type (`np.dtype('float64')`), an `MEnum('item1', ..., 'itemN')`, a predefined
metainfo type (`Datetime`, `JSON`, `File`, ...), or another section or quantity to define
a reference type.
- `shape`: defines the dimensionality of the quantity. Examples are: `[]` (number),
`['*']` (list), `[3, 3]` (3 by 3 matrix), `['n_elements']` (a vector of length defined by
another quantity `n_elements`).
- `unit`: a physical unit. We use [Pint](https://pint.readthedocs.io/en/stable/){:target="_blank"} here. You can
use unit strings that are parsed by Pint, e.g. `meter`, `m`, `m/s^2`. As a convention the
NOMAD Metainfo uses only SI units.

### SubSection

A subsection defines a named property of a section that refers to another section. It
allows to define that a section that contains another section.

- `sub_section`: (aliases `section_def`, `sub_section_def`) defines the section that can
be contained.
- `repeats`: a boolean that determines whether the subsection relationship allows multiple sections
or only one.

### References and Proxies

Besides creating hierarchies with subsections (e.g. tree structures), the metainfo
also allows one to create a reference within a section that points to either another section or a quantity
value:

```python
class Calculation(MSection):
    system = Quantity(type=System.m_def)
    atom_labels = Quantity(type=System.atom_labels)

calc = Calculation()
calc.system = run.systems[-1]
calc.atom_labels = run.systems[-1]
```

To define a reference, define a normal quantity and simply use the section or quantity
you want to refer to as type. Then you can assign respective section instances
as values.

In Python memory, quantity values that reference other sections simply contain a
Python reference to the respective *section instance*. However, upon serializing/storing
metainfo data, these references have to be represented differently.

Value references work a little differently. When you read a value reference, it behaves like
the reference value. Internally, we do not store the values, but instead a reference to the
section that holds the referenced quantity is stored. Therefore, when you want to
assign a value reference, use the section with the quantity and not the value itself.
<!-- TODO Add a simply example here -->

References are serialized as URLs. There are different types of reference URLs:

- `#/run/0/calculation/1`: a reference in the same Archive
- `/run/0/calculation/1`: a reference in the same archive (legacy version)
- `../upload/archive/mainfile/{mainfile}#/run/0`: a reference into an Archive of the same upload
- `/entries/{entry_id}/archive#/run/0/calculation/1`: a reference into the Archive of a different entry on the same NOMAD installation
- `/uploads/{upload_id}/archive/{entry_id}#/run/0/calculation/1`: similar to the previous one but based on uploads
- `https://myoasis.de/api/v1/uploads/{upload_id}/archive/{entry_id}#/run/0/calculation/1`: a global reference towards a different NOMAD installation (Oasis)

The host and path parts of URLs correspond with the [NOMAD API](../programmatic/api.md). The anchors are paths from the root section of an Archive, over its subsections, to the referenced section or quantity value. Each path segment is the name of the subsection or an index in a repeatable subsection: `/system/0` or `/system/0/atom_labels`.

References are automatically serialized by `:py:meth:MSection.m_to_dict`. When de-serializing
data with `:py:meth:MSection.m_from_dict` these references are not resolved right away,
because the reference section might not yet be available. Instead references are stored
as `:class:MProxy` instances. These objects are automatically replaced by the referenced
object when a respective quantity is accessed.

If you want to define references, it might not be possible to define the referenced
section or quantity beforehand, due to the way Python definitions and imports work. In these
cases, you can use a proxy to reference the reference type. There is a special proxy
implementation for sections:

```python
class Calculation(MSection):
    system = Quantity(type=SectionProxy('System')
```

The strings given to `SectionProxy` are paths within the available definitions.
The above example works, if `System` is eventually defined in the same package.

### Categories

In the old metainfo this was known as *abstract types*.

Categories are defined with Python classes that have `:class:MCategory` as base class.
Their name and description are taken from the name and docstring of the class. An example
category looks like this:

``` python
class CategoryName(MCategory):
    ''' Category description '''
    m_def = Category(links=['http://further.explanation.eu'], categories=[ParentCategory])
```

### Packages

Metainfo packages correspond to Python packages. Typically your metainfo Python files should follow this pattern:
```python
from nomad.metainfo import Package

m_package = Package()

# Your section classes and categories

m_package.__init_metainfo__()
```

## Adding Python schemas to NOMAD

The following describes how to integrate new schema modules into the existing code according
to best practices.

### Schema super structure

You should follow the basic [developer's getting started](../develop/setup.md) to setup a development environment. This will give you all the necessary libraries and allows you
to place your modules into the NOMAD code.

The `EntryArchive` section definition sets the root of the archive for each entry in
NOMAD. It therefore defines the top level sections:

- `metadata`: all "administrative" metadata (ids, permissions, publish state, uploads, user metadata, etc.)
- `results`: a summary with copies and references to data from method specific sections. This also
presents the [searchable metadata](../develop/search.md).
- `workflows`: all workflow metadata
- Method-specific subsections: e.g. `run`. This is were all parsers are supposed to
add the parsed data.
<!-- TODO Update!!! -->

The main NOMAD Python project includes Metainfo definitions in the following modules:

- `nomad.metainfo`: defines the Metainfo itself. This includes a self-referencing schema. E.g. there is a section `Section`, etc.
- `nomad.datamodel`: defines the section `metadata` that contains all "administrative"
metadata. It also contains the root section `EntryArchive`.
- `nomad.datamodel.metainfo`: defines all the central, method specific (but not parser specific) definitions.
For example the section `run` with all the simulation definitions (computational material science definitions)
that are shared among the respective parsers.

### Extending existing sections

Parsers can provide their own definitions. By convention, these are placed into a
`metainfo` sub-module of the parser Python module. The definitions here can add properties
to existing sections (e.g. from `nomad.datamodel.metainfo`). By convention, use a `x_mycode_`
prefix. This is done with the
`extends_base_section` [Section property](#sections). Here is an example:
<!-- ? Do we want to encourage this as best practice in the future? -->

```py
from nomad.metainfo import Section
from nomad.datamodel.metainfo.workflow import Workflow

class MyCodeRun(Workflow)
    m_def = Section(extends_base_section=True)
    x_mycode_execution_mode = Quantity(
        type=MEnum('hpc', 'parallel', 'single'), description='...')
```

### Schema conventions

- Use lower snake case for section properties; use upper camel case for section definitions.
- Use a `_ref` suffix for references.
- Use subsections rather than inheritance to add specific quantities to a general section.
E.g. the section `workflow` contains a section `geometry_optimization` for all geometry optimization specific
workflow quantities.
- Prefix parser-specific and user-defined definitions with `x_name_`, where `name` is the
short handle of a code name or other special method prefix.
<!-- TODO add case examples to the reference pages and add corresponding links here and throughout -->

## Use Python schemas to work with data

### Access structured data via API

The [API section](../programmatic/api.md#access-archives) demonstrates how to access an Archive, i.e.
retrieve the processed data from a NOAMD entry. This API will give you JSON data likes this:

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


### Wrap data with Python schema classes

In Python, JSON data is typically represented as nested combinations of dictionaries
and lists. Of course, you could work with this right away. To make it easier for Python
programmers, the [NOMAD Python package](../programmatic/pythonlib.md) allows you to use this
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

### Access structured data via the NOMAD Python package

The NOMAD Python package provides utilities to [query large amounts of
archive data](../programmatic/archive_query.md). This uses the built-in Python schema classes as
an interface to the data.

## Custom normalizers

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
support schemas with multiple inheritance. In order to control the order by which the
normalizations are executed, one can define `normalizer_level` which is set to 0 by
default. The normalize functions are always called for any sub section before the parent section.
However, the order for any sections on the same level will be from low values of `normalizer_level`
to high.

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
## Schema plugin metadata
{{pydantic_model('nomad.config.models.plugins.Schema')}}

## Pre-defined schemas in NOMAD
Several schemas are currently being developed in NOMAD for various metadata.
The following lists these projects:

| Description                  | Project url                                                                |
| ---------------------------- | -------------------------------------------------------------------------- |
| simulation run               | <https://github.com/nomad-coe/nomad-schema-plugin-run.git>                 |
| simulation data              | <https://github.com/nomad-coe/nomad-schema-plugin-simulation-data.git>     |
| simulation workflow          | <https://github.com/nomad-coe/nomad-schema-plugin-simulation-workflow.git> |
| NEXUS                        | <https://github.com/FAIRmat-NFDI/pynxtools.git>                            |
| synthesis                    | <https://github.com/FAIRmat-NFDI/AreaA-data_modeling_and_schemas.git>      |
| material processing          | <https://github.com/FAIRmat-NFDI/nomad-material-processing.git>            |
| measurements                 | <https://github.com/FAIRmat-NFDI/nomad-measurements.git>                   |
