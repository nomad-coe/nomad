# Extending the Archive and Metainfo

In [Using the Archive and Metainfo](archive.md), we learned what the Archive and Metainfo
are. It also demonstrated the Python interface and how to use it on Archive data. The
Metainfo is written down in Python code as a bunch of classes that define *sections*
and their properties. Here, we will look at how the Metainfo classes work and how
the metainfo can be extended with new definitions.

## Starting example

```python
from nomad.metainfo import MSection, Quantity, SubSection, Units

class System(MSection):
    '''
    A system section includes all quantities that describe a single a simulated
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

We define simple metainfo schema with two *sections* called `System` and `Run`. Sections
allow to organize related data into, well, *sections*. Each section can have two types of
properties: *quantities* and *sub-sections*. Sections and their properties are defined with
Python classes and their attributes.

Each *quantity* defines a piece of data. Basic quantity attributes are its `type`, `shape`,
`unit`, and `description`.

*Sub-sections* allow to place section into each other and therefore allow to form containment
hierarchies or sections and the respective data in them. Basic sub-section attributes are
`sub_section`(i.e. a reference to the section definition of the sub-section) and `repeats`
(determines if a sub-section can be contained once or multiple times).

The above simply defines a schema, to use the schema and create actual data, we have to
instantiate the above classes:

```python
run = Run()
system = run.m_create(System)
system.n_atoms = 3
system.atom_labels = ['H', 'H', 'O']

print(system.atom_labels)
print(n_atoms = 3)
```

Section *instances* can be used like regular Python objects: quantities and sub-sections
can be set and access like any other Python attribute. Special meta-info methods, starting
with `m_` allow us to realize more complex semantics. For example `m_create` will
instantiate a sub-section and add it to the *parent* section in one step.

Another example for an `m_`-method is:

```python
run.m_to_json(indent=2)
```

This will serialize the data into JSON:

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

## Definitions and Instances

As you already saw in the example, we first need to define what data can look like (schema),
before we can actually program with data. Because schema and data are often discussed in the
same context, it is paramount to clearly distingish between both. For example, if we
just say "system", it is unclear what we refer to. We could mean the idea of a system, i.e. all
possible systems, a data structure that comprises a lattice, atoms with their elements and
positions in the lattice. Or we mean a specific system of a specific calculation, with
a concrete set of atoms, real numbers for lattice vectors and atoms positions as concrete data.

The NOMAD Metainfo is just a collection of definition that describe what materials science
data could be (a schema). The NOMAD Archive is all the data that we extract from all data
provided to NOMAD. The data in the NOMAD Archive follows the definitions of the NOMAD metainfo.

Similarely, we need to distingish between the NOMAD Metainfo as a collection of definitions
and the Metainfo system that defines how to define a section or a quantity. In this sense,
we have a three layout model, were the Archive (data) is an instance of the Metainfo (schema)
and the Metainfo is an instance of the Metainfo system (schema of the schema).

This documentation describes the Metainfo by explaining the means of how to write down definitions
in Python. Conceptually we map the Metainfo system to Python language constructs, e.g.
a section definition is a Python class, a quantity a Python property, etc. If you are
familiar with databases, this is similar to what an object relational mapping (ORM) would do.


### Common attributes of Metainfo Definitions

In the example, you already saw the basic Python interface to the Metainfo. *Sections* are
represented in Python as objects. To define a section, you write a Python classes that inherits
from `MSection`. To define sub-sections and quantities you use Python properties. The
definitions themselves are also objects derived from classes. For sub-sections and
quantities, you directly instantiate :class`SubSection` and :class`Quantity`. For sections
there is a generated object derived from :class:`Section` that is available via
`m_def` from each *section class* and *section instance*.

These Python classes that are used to represent metainfo definitions form an inheritance
hierarchy to share common properties

- `name`, each definition has a name. This is typically defined from the corresponding
Python. E.g. a sections class name, becomes the section name; a quantity gets the name
from its Python property, etc.
- `description`, each definition should have one. Either set it directly or sue *doc strings*
- `links`, a list of useful internet references.
- `more`, a dictionary of custom information. All additional `kwargs` that are set when
creating a definition are added to `more`.

### Sections

Sections are defined with Python classes that extend `MSection` (or other section classes)

- `base_sections` are automatically taken from the Python class's base classes.
- `extends_base_section` a boolean that determines the inheritance. If this is `False`,
normal Python inheritance implies and this section will inherit all properties (sub-sections,
quantities) from all base classes. If this is `True`, all definition in this section
will be added to the properties of the base class section. This allows to extend existing
sections with additional properties.

### Quantities
Quantity definitions are the main building block of meta-info schemas. Each quantity
represents a single piece of data. Quantities can define:

- A `type`, this can be a primitive Python type (`str`, `int`, `bool`), a numpy
data type (`np.dtype('float64')`), an `MEnum('item1', ..., 'itemN')`, a predefined
metainfo type (`Datetime`, `JSON`, `File`, ...), or another section or quantity to define
a reference type.
- A `shape` that defines the dimensionality of the quantity. Examples are: `[]` (number),
`['*']` (list), `[3, 3]` (3 by 3 matrix), `['n_elements']` (a vector of length defined by
another quantity `n_elements`).
- A physics `unit`. We use [Pint](https://pint.readthedocs.io/en/stable/) here. You can
use unit strings that are parsed by Pint, e.g. `meter`, `m`, `m/s^2`. As a convention the
metainfo uses only SI units.

### Sub-Section
A sub-section defines a named property of a section that refers to another section. It
allows to define that a section can contain another section.

- `sub_section` (aliases `section_def`, `sub_section_def`) defines the section that can
be contained.
- `repeats` is a boolean that determines if the sub-section relationship allows for multiple
or only one section.

### References and Proxies

Beside creating hierarchies (e.g. tree structures) with sub sections, the metainfo
also allows to create cross references between sections and other sections or quantity
values:

```python
class Calculation(MSection):
    system = Quantity(type=System.m_def)
    atom_labels = Quantity(type=System.atom_labels)

calc = Calculation()
calc.system = run.systems[-1]
calc.atom_labels = run.systems[-1]
```

To define a reference, you define a normal quantity and simply use the section or quantity
that you want to reference as type. Then you can assign respective section instances
as values.

When in Python memory, quantity values that reference other sections simply contain a
Python reference to the respective *section instance*. However, upon serializing/storing
metainfo data, these references have to be represented differently.

Value references are a little different. When you read a value references it behaves like
the references value. Internally, we do not store the values, but a reference to the
section that holds the referenced quantity. Therefore, when you want to
assign a value reference, you use the section with the quantity and not the value itself.

References are serialized as URLs. There are different types of reference URLs:

- `#/run/0/calculation/1`, a reference in the same Archive
- `/run/0/calculation/1`, a reference in the same archive (legacy version)
- `../upload/archive/mainfile/{mainfile}#/run/0`, a reference into an Archive of the same upload
- `/entries/{entry_id}/archive#/run/0/calculation/1`, a reference into the Archive of a different entry on the same NOMAD installation
- `/uploads/{upload_id}/archive/{entry_id}#/run/0/calculation/1`, similar but based on uploads
- `https://myoasis.de/api/v1/uploads/{upload_id}/archive/{entry_id}#/run/0/calculation/1`, a global reference towards a different NOMAD installation (Oasis)

The host and path parts of URLs correspond with the NOMAD API. The anchors are paths from the root section of an Archive, over its sub-sections, to the referenced section or quantity value. Each path segment is the name of the sub-section or an index in a repeatable sub-section: `/system/0` or `/system/0/atom_labels`.

References are automatically serialized by :py:meth:`MSection.m_to_dict`. When de-serializing
data with :py:meth:`MSection.m_from_dict` these references are not resolved right away,
because the references section might not yet be available. Instead references are stored
as :class:`MProxy` instances. These objects are automatically replaced by the referenced
object when a respective quantity is accessed.

If you want to defined references, it might not be possible to define the referenced
section or quantity before hand, due to how Python definitions and imports work. In these
cases, you can use a proxy to reference the reference type. There is a special proxy
implementation for sections:

```python
class Calculation(MSection):
    system = Quantity(type=SectionProxy('System')
```

The strings given to `SectionProxy` are paths within the available definitions.
The above example works, if `System` is eventually defined in the same package.

### Categories

In the old meta-info this was known as *abstract types*.

Categories are defined with Python classes that have :class:`MCategory` as base class.
Their name and description is taken from the class's name and docstring. An example
category looks like this:

``` python
class CategoryName(MCategory):
    ''' Category description '''
    m_def = Category(links=['http://further.explanation.eu'], categories=[ParentCategory])
```

### Packages

Metainfo packages correspond with Python packages. Typically you want your metainfo
Python files follow this pattern:
```python
from nomad.metainfo import Package

m_package = Package()

# Your section classes and categories

m_package.__init_metainfo__()
```

## Adding definition to the existing metainfo schema

Now you know how to define new sections and quantities, but how should your additions
be integrated in the existing schema and what conventions need to be followed?

### Metainfo schema super structure

The `EntryArchive` section definition set the root of the archive for each entry in
NOMAD. It therefore defines the top level sections:

- `metadata`, all "administrative" metadata (ids, permissions, publish state, uploads, user metadata, etc.)
- `results`, a summary with copies and references to data from method specific sections. This also
presents the [searchable metadata](search.md).
- `workflows`, all workflow metadata
- Method specific sub-sections, e.g. `run`. This is were all parsers are supposed to
add the parsed data.

The main NOMAD Python project include Metainfo definitions in the following modules:

- `nomad.metainfo` Defines the Metainfo itself. This includes a self-referencing schema
of itself. E.g. there is a section `Section`, etc.
- `nomad.datamodel` Mostly defines the section `metadata` that contains all "administrative"
metadata. It also contains the root section `EntryArchive`.
- `nomad.datamodel.metainfo` Defines all the central, method specific (but not parser specific) definitions.
For example the section `run` with all the simulation (computational material science definitions)
definition that are shared among the respective parsers.

### Extending existing sections

Parsers can provide their own definitions. By conventions those are places into a
`metainfo` sub-module of the parser Python module. The definitions here can add properties
to existing sections (e.g. from `nomad.datamodel.metainfo`). By convention us a `x_mycode_`
prefix. This is done with the
`extends_base_section` [Section property](#sections). Here is an example:

```py
from nomad.metainfo import Section
from nomad.datamodel.metainfo.simulation import Method

class MyCodeRun(Method)
    m_def = Section(extends_base_section=True)
    x_mycode_execution_mode = Quantity(
        type=MEnum('hpc', 'parallel', 'single'), description='...')
```

### Metainfo schema conventions

- Use lower snake case for section properties; use upper camel case for section definitions.
- Use a `_ref` suffix for references.
- Use sub-sections rather than inheritance to add specific quantities to a general section.
E.g. section `workflow` contains a section `geometry_optimization` for all geometry optimization specific
workflow quantities.
- Prefix parser specific and custom definitions with `x_name_`. Where `name` is the
short handle of a code name or other special method prefix.
