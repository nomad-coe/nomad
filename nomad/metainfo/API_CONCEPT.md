**! This is not yet aligned with the ideas of CONCEPT.md !**

# A new metainfo schema, interface, and *file formats* support

This is a design document (later documentation) for a re-implementation of nomad's old
meta-info *system*. Here *system* refers to the set of software components necessary to
define (meta-)data structures and types, CRUD data according to schema via an abstract
interface, and manage data in various resources
(diff. file formats, databases, web resources, etc)

## Example use cases

We'll try to explain/design the system through a serious of use cases first. The
respective examples are centered around a hypothetical python library that works as
the meta-info systems interface.

### Scientists using nomad parsers

Imagine a student that uses VASP to simulate a system and wants to plot the density of states.
She should be able to access the respective energies with:

```python
import nomad
nomad.parse('TiO3.xml').run.single_configuration_calculation.dos.dos_energies.values
```

If she does not know what exactely `dos_energies` refers to:
```python
my_calc = nomad.parse('TiO3.xml')
my_calc.run.single_configuration_calculation.dos.dos_energies.definition
```

It should give you something like:
```json
{
    "description": "Array containing the set of discrete energy values for the density ...",
    "share": ["number_of_dos_values"],
    "type": "float",
    "units": "J"
}
```

But the units should actually not be fixed:

```python
my_calc.run.system.atom_positions.convert(nomad.units.angstrom).values
```

Values can be regular python lists or np arrays:
```python
my_calc.run.system.atom_positions.np_array
```

In the cases above, `system` is a list of systems. Therefore, the lines should return
a list of actual values (e.g. a list of position matrices/np_arrays). To access
a particular system:
```python
my_calc.run.system[0].atom_positions.np_array
```

To create more complex dict structures:
```python
from nomad import x
my_calc.run.system(atom_positions=x.np_array, atom_labels=x.values, lattice_vector=x.np_array)
```
Should return a list of dictionaries with `atom_position`, `atom_labels`, `lattice_vector` as keys.
The `x` acts as a surrogate for the retrived meta-info objects and everything accessed
on `x`, will be accessed on the actual values.

You can also create recursive *GraphQL* like queries:
```python
my_calc.run.system(atom_labels=x.values, symmetry=x(spacegroup=x.value)))
```

Or if another syntax is prefered:
```python
my_calc.run.system({
    atom_labels: x.values,
    symmetry: {
        spacegroup: x.value,
        hall_number: x.value
    }
})
```

### Working with uploaded data

There needs to be support for various resources, e.g. resources on the web, like the
nomad archive:
```python
nomad.archive(upload_id='hh41jh4l1e91821').run.system.atom_labels.values
```

This can also be used to extend queries for calculations with queries for certain
data points:
```python
nomad.archive(user='me@email.org', upload_name='last_upload').run.system(
    atom_labels=x.values, atom_positions=x.convert(units.angstrom).np_array)
```
In this case, it will return a generator that hides any API pagination. This this is now
not a single run, it would generate a list (runs) or lists (systems) of dicts
(with labels, and positions).


### A parser creating data
The CoE NOMAD parsing infrastructure used the concept of a *backend* that acted as an
abstract interface to enter new data to the nomad archive. We like to think in terms of
*resources*, where a resource can represent various storage media
(e.g. in-memory, hdf5 file, something remote, etc.). Lets assume that the parser gets
such a *resource* from the infrastructure to populate it with new data. Lets call the
*resource* `backend` for old times sake:

```python
run = backend.run()
system = run.system()
system.atom_labels.values = ['Ti', 'O', 'O']
system.atom_positions.values = [a1, a2, a3]
system.close()
system = run.system()
# ...
resource.close()
```

This basically describes the write interface. Of course parsers should also allow to read
back all the properties they entered so far:
```python
system.atom_labels.values = ['Ti', 'O', 'O']
system.atom_labels.values
```

The old backend allowed to build arrays piece by piece. This could also be possible here:
```python
positions = system.atom_positions.create()
positions.add(a1)
positions.add(a2)
```

## Core concepts


### Conventions

When mapping the following concepts to python implementations, we use the prefix `m_` on
all methods and attributes that might conflict with user given names for meta info
definitions.

### Resources

A *resource* refers to anything that can be used to *hold* data. This can be basic
in python memory, a JSON file, an HDF5 file, a search index, a mongodb, or a remote
resource that is accessible via REST API. Respective *resource* classes and their objects
are used to parameterize access to the data. For example, to access a JSON file a
file path is required, to store something in mongodb a connection to mongo and a collection
is necessary, to read from an API an endpoint and possible parameters have to be configured.

Beyond parameterized data access, all resources offer the same interface to navigate,
enter, or modify data. The only exception are readonly resources that do not allow
to add or modify data.

The creation of resources could be allowed directly from a `nomad` package to create a
very simple interface:
```python
nomad.parse('vasp_out.xml')
nomad.open('archive.json')
nomad.open('archive.hdf5')
nomad.connect('mongodb://localhost/db_name', 'calculations')
nomad.connect('https://nomad.fairdi.eu/archive/3892r478323/23892347')
```

The various resource implementations should offer the same interface. Necessary methods are
- `close(save: bool = True)` to close a open file or connection
- `save()` to save changes

### Data Objects

When navigating the contents of a resource (e.g. via `nomad.open('archive.json').run.system.atom_labels`),
we start from a resource object (`nomad.open('archive.json')`) and pass various *data objects* (`.run.system.atom_labels`).
There are obviously different types of data objects, i.e. *sections* (like `run`, `system`) and
*properties* (like `atom_labels`). Sections and properties have to offer different interfaces.
Sections need to allow to access and create subsections and properties. Properties have to allow to
access, set, or modify the stored data.

Independent of the object type, all data objects should allow to navigate to the definition (e.g. `run.m_definition`, `run.system.atom_labels.definition`).

Navigation uses the user defined names of meta-info definitions for sections and properties.

Section object need to support:
- access to subsection via subsection name
- access of properties via property name
- array access for repeatable sections
- navigation to its containing section: `.m_def`
- allow to create/(re-)open subsections via calling the subsection name as a method: `.system()`
- close a section so that the underlying resource implementation can potentially remove the section from memory and write it to a database/.hdf5 file
- the *GraphQL* like access methods with dictionary to specify multiple sub-sections

Property objects
- access to values, depending on the shape and desired representation: `.value`, `.values`, `.np_array`
- set values: `.values = ['Ti', 'O']`
- create arrays and add values: `.atom_positions().add(...)`
- convert units: `.convert(nomad.units.angstrom).values` or `.convert(nomad.units.angstrom).values = ...`

### References and Resource Sets

The meta-info allows to define references as a special property type. In this case the
values of a property are references to other sections. This can be sections in the same
resource or even other resources (across files). From an interface perspective, references
can be set like any other property:

```python
system = run.system()
calc = run.single_configuration_calculation()
calc.system.value = system
```

Within the underlying resource, references are represented as URLs and these can be
`local` or `remote`. URL paths define the position of a section in a resource. Local
URLs only contain the path, remote URLs also contain a part that identifies the resource.

The local reference in the previous example would be `/run/system/0`, referring to the first
system in the run.

If we create two different resource:
```python
systems = nomad.open('systems.json')
calcs = nomad.open('calc.json')
system_1 systems.system()
system_2 systems.system()
calc = calcs.run().single_configuration_calculation()
calc.system = system_2
systems.close()
calcs.close()
```

the reference in the `calc.json` would be `file://systems.json/system/1`.

When accessing a resource, other resources will be accessed on demand, when ever a reference
needs to be resolved. The library can keep track of all accessed (and therefore open)
resources through a global resource set. In more involved use-cases, it might be desirable
that users can control resource sets via python.

### Schemas and Definitions

#### Introduction to schemas and definitions

A schema comprises a set of definitions that define rules which govern what data does
adhere to the schema and what not. We also say that we validate data against a schema to
check if the data follows all the rules. In this sense, a schema defines an unlimited
set of possible data that can be expressed in this schema.

The definitions a schema can possibly contain is also govern by rules and these rules are also
defined in a schema and this schema would be the schema of the schema. To be even
more confusing, a schema can be the schema of itself. Meaning we can use the same set of
definitions to formally define the definitions themselves. But lets start with an
informal definition of schema elements.

The following *elements* or kinds of definitions are used in metainfo schemas.

#### Elements
We call all kinds of definitions and their parts elements. Element is an abstract concept,
in the sense that you cannot define elements directly. The abstract definition for elements
barely defines a set of properties that is then shared by other elements.

These properties include:
- the elements *name* in the sense of a python compatible identifier (only a-bA-Z0-9_ and conventionally in camel case)
- a human readable description, potentially in markdown format
- a list of tags
- a reference to a *section* definition. This denotes that instances of this element definition can be contained
in instances of the referenced (parent) section definition.

#### Sections
A section definition is a special element. Sections will be used to create
hierarchical structures of data. A section can contain other section and a set of properties.
A section definition has the following properties: name, description, (parent) section (as all element definitions have), plus
- a boolean *abstract* that determine if this section definition can be instantiated
- a boolean *repeats* that determines if instances of this section can appear multiple times in their respective parent section
- a references to another section definition called *extends* that denotes that instance of this definition can
contain instances of all the element definitions that state the extended section definition as their (parent) section.

#### Propertys
A property definition is a special element. A property can be contained in a
section and a property has a value. The type of the property values is determined by its property
definition. Property value can be scalar, vectors, matrices, or higher-dimensional matrices.
Each possible dimension can contain *primitive values*. Possible primitive value types are
numerical types (int, float, int32, ...), bool, str, or references to other sections. Types
of references are defined by referencing the respective section definition.

A property definition has all element properties (name, description, (parent) section, tags) and has the following properties
- *type* the data type that determine the possible values that the various dimensions can contain
- a *shape* that determines how many primitive values each dimension of the properties value can contain. This can be a fix integer,
a *wildcard* like 1..n or 0..n, or a references to another property with one dimension with a single int.
- *units* is a list of strings that determines the physical unit of the various dimensions.

#### Packages
Packages are special elements. Packages can contain section definitions
and property definitions. They have an additional property:
- *dependencies* that list references to those packages that have definitions which complete the definitions
in this package.

#### annotation
- name, value


If we use the python interface described by example above, we can use it to define the
metainfo schema schema in itself. Be aware that this might break your brain:
```python
from nomad import metainfo as mi
from nomad.metainfo import types, dimensions

resource = mi.open('metainfo.json')
metainfo = resource.package(
    name='metainfo',
    description='The nomad meta-info schema schema')

element = metainfo.section(name='element', abstract=True)
metainfo.property(name='name', type=str, section=element)
metainfo.property(name='description', type=str, section=element)
super_section = metainfo.property(name='section', section=element)

package = metainfo.section(name='package', extends=element)
metainfo.property(name='dependency', section=package, types=[types.reference(package)], shape=[dimensions.0_n])

element.section = package

section = metainfo.section(name='section', repeats=True, extends=element, section=package)

super_section.types = [nomad.types.reference(section)]

metainfo.property(name='abstract', type=bool, section=section)
metainfo.property(name='repeats', type=bool, section=section)
metainfo.property(name='extends', type=types.reference(section), shape=[dimensions.0_n], section=section)

property = metainfo.section(name='property', section=package, repeats=True)
metainfo.property(name=types, type=types.type], section=property)
metainfo.property(name=shape, type=types.union(types.reference(property), int, types.dimension), shape=[dimensions.0_n], section=property)
metainfo.property(name=units, type=str, shape=[dimensions.0_n])

annotation = metainfo.section(name='annotation', section=element, repeats=True)
metainfo.property(name='key', type=str, section=annotations)
metainfo.property(name='value', type=str, section=annotations)
```

In a similar manner, you can use the metainfo python interface to create metainfo
schema definitions:
```python
resource = mi.open('common.json')
common = resource.package(name='common')

run = common.section(name='run')
system = common.section(name='system', section='run')
number_of_atoms = common.property(name='number_of_atoms', section=system)
common.property(name='atom_labels', types=[str], shape=[number_of_atoms], section=system)
common.property(name='atom_positions', types=[])

```
