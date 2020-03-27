# NOMAD MetaInfo concept

## History

The NOMAD MetaInfo was devised within the first NOMAD CoE; over 2000 quantities have
been defined in this *old MetaInfo*. The experience with this system revealed the following drawbacks:

- The Python libraries that allow to use the MetaInfo are non pythonic and incomplete.
- The MetaInfo is only used for the archive, not for the encyclopedia and repository data.
- There is no direct support to map MetaInfo definitions to DI technologies (databases, search indices, APIs).
- There is no support for namespaces. MetaInfo names are cumbersome. This will not scale to expected levels of FAIRmat metadata.
- MetaInfo packages are not version controlled. They are part of the same git and do not belong to the independently evolving parsers. This does not allow for "external" parser development and makes it hard to keep versions consistent.
- The MetaInfo is defined in JSON. The syntax is inadequate, checks are not immediate.
- Attempts to revise the MetaInfo have failed in the past.

## Goals

### Common language to define physics (meta-)data quantities and their relationships

The *physics quantities* part includes
- each quantity MAY have a physics *unit*
- each quantity MUST have a *shape* that precisely define vectors, matrices, tensors, and their dimensions
- each quantity MAY have a numpy dtype that allows to map physics data to numpy arrays

The *relationship* parts entails:
- hierarchies for quantity *values* (e.g. *sections*)
- hierarchies for quantity *definition* (e.g. *categories*, former *abstract types*)
- *derived* quantities that can be computed from other quantities
- *synonyms* as a special trivial case for derived quantities
- *shapes* might also define a type of relationship through one quantity being the dimension of another
- *references* between sections, via quantities that have a section definition as type

In addition there are the *typical* data-type definition (schema, ontology, ...) features:
- names/namespaces
- modularization (i.e. Metainfo packages)
- extensions: section inheritance, sections that add to other sections after definition
- documentation
- basic primitive types (int, string, bool)
- simple compound types (lists, dictionaries, unions)
- MAYBE an event mechanism

### Complex, evolving, extendable packages of quantities

There are a lot of quantities, and they need to be organized. There are three mechanisms
to organize quantities:
- *Packages* (a.k.a modules) allow to modularize large sets of quantities, e.g. one package per code
- *Sections* allow to organize quantity values into containment (a.k.a whole-part, parent-child) hierarchies, e.g. `system` *contains* all quantity values that describe the simulated system.
- *Categories* allow to organize quantity definitions via generalization (a.k.a specialization, inheritance) relationships, e.g. `atom_labels` and `formula_hill` (*special*) both express `chemical_composition` (*general*)

Quantities and their relationships change over time. This requires (at least) a versioning mechanism to track changes and reason whether a pieces of data adheres to a certain version of the MetaInfo or not.

The MetaInfo needs to be extendable. It must be possible to add *packages*, quantities in new *packages* must be addable to existing sections and categories. Existing sections must be extendable. It must be possible to develop and version packages independently.

### Mappings to DI technologies

The core of the MetaInfo is about defining data and their physics. But in the end, the data needs to be managed with DI components, such as file formats, databases, search indices, onotology tools, APIs, GUIs, programming languages, etc. While all these tools come with their own ways of defining data, it can be cumbersome to manually map the MetaInfo to the corresponding DI technology. Furthermore, this usually comprises both mapping definitions and transforming values.

The MetaInfo will allow for quantity *annotations*. Annotations allow to add additional
information to quantity definitions that carry the necessary information to automatically map/transform definitions and their values to underlying DI components. Annotations can be easily stripped/filtered to present the MetaInfo either clean or under technology specific lenses.

### Intuitive programming interface to create, access, and use (meta-)data defined with the NOMAD MetaInfo

While MetaInfo definitions and MetaInfo values should have a *native* serialization format (JSON), the primary interface to deal with definitions and data should be made from programming language (Python) primitives. By the way, both things are basically just mappings from the logical MetaInfo into concrete technologies (i.e. JSON, Python).

As a programming language, Python has a far richer set of syntax to define and use data
than JSON has. We should use this. It was not used for definitions in the NOMAD CoE, and
the *backend*s for data were designed for creating data only and not very *pythonic*.

## Concepts for a new NOMAD MetaInfo

A schema comprises a set of definitions that define rules which govern what data does
adhere to the schema and what not. We also say that we validate data against a schema to
check if the data follows all the rules. In this sense, a schema defines an unlimited
set of possible data that can be expressed in this schema.

The definitions a schema can possibly contain is also govern by rules and these rules are also
defined in a schema and this schema would be the schema of the schema. To be even
more confusing, a schema can be the schema of itself. Meaning we can use the same set of
definitions to formally define the definitions themselves. But lets start with an
informal definition of schema elements.

### Conventions

When mapping the following concepts to python implementations, we use the prefix `m_` on
all methods and attributes that might conflict with user given names for meta info
definitions.

### Definition
We call all kinds of definitions and their parts definitions. Definition is an abstract concept,
in the sense that you cannot define definitions directly. The abstract definition for elements
barely defines a set of properties that is then shared by other elements.

These properties include:
- the elements *name* in the sense of a python compatible identifier (only a-bA-Z0-9_ and conventionally in camel case)
- a human readable description, potentially in markdown format
- a list of categories
- annotations that attach possible metainfo extensions (db, search support, etc.) to the definition

`Definition` is the abstract base for all definitions in the MetaInfo.

- `name`, a string
- `description`, a string
- `links`, a list of URLs
- `categories`, a list of references to category definitions
- `annotations`, a list of `Annotations`

- *derived*: `qualified_name`


### Property

`Property` is a special `Definition` and an abstract base for section properties.
Properties define what data a section instance can hold. Properties are mapped to Python
*descriptors*.

- `section` specialized `parant` relation with the containing `Section`


#### SubSections

`SubSection` is a special `Property` that defines that a section instance can **contain**
the instances of a sub section.

- `sub_section` reference to the `Section` definition for the children
- `repeats` is a boolean that determines if this sub section can be contain only once of multiple times

- *constraint*: sub sections are not circular


### Quantities (incl. dimensions, incl. references)

A Quantity definition is a special definition. A quantity can be contained in a
section and a quantity has a value. The type of the quantity values is determined by its quantity
definition. Quantity value can be scalar, vectors, matrices, or higher-dimensional matrices.
Each possible dimension can contain *primitive values*. Possible primitive value types are
numerical types (int, float, int32, ...), bool, str, or references to other sections. Types
of references are defined by referencing the respective section definition.

A quantity definition has all definition quantities (name, description, ...) and has the following properties
- *type* the data type that determine the possible values that the various dimensions can contain
- a *shape* that determines how many primitive values each dimension of the quantities value can contain. This can be a fix integer,
a *wildcard* like 1..n or 0..n, or a references to another property with one dimension with a single int.
- *units* is a list of strings that determines the physical unit of the various dimensions.

A `Quantity` definition is a special and concrete `Property` definition:

- `shape`, a list of either `int`, references to a dimension (quantity definition), or limits definitions (e.g. `'1..n'`, `'0..n'`.)
- `type`, a primitive or MEnum type
- `unit`, a (computed) units, e.g. `units.F * units.m`
- `derived_from`, a list of references to other quantity definitions
- `synonym`, a reference to another quantity definition

*Dimensions* are quantity definitions with empty shape and int type.

- *constraint*: `synonym`, `derived_from`, and dimensions come from the same section


### Sections (incl. references)

A section definition is a special definition. Sections will be used to create
hierarchical structures of data. A section can contain other section and a set of properties.
A section definition has the following properties: name, description, (parent) section (as all element definitions have), plus
- a boolean *abstract* that determine if this section definition can be instantiated
- a boolean *repeats* that determines if instances of this section can appear multiple times in their respective parent section
- a references to another section definition called *extends* that denotes that instance of this definition can
contain instances of all the element definitions that state the extended section definition as their (parent) section.

A `Section` is a special and concrete `Definition`.

- `adds_to`, a reference to another section definition. All quantities of this *pseudo* section are added to the given section. (Might not be necessary)
- `repeats`, a boolean
- `extends`, list of reference to other section definitions. This section automatically inherits all quantities of the other sections. (Might not be necessary)

- *derived*: `all_sub_sections`, all sub sections, included added and inherited ones, by name
- *derived*: `all_quantities`, all quantities, included added and inherited ones, by name
- *derived*: `all_properties`, all properties, included added and inherited ones, by name

- *constraint*: `extends` is not circular
- *constraint*: `adds_to` is not circular
- *constraint*: all quantities that have *this* (or an `extends`) section as `section` have unique names

`Section`s are mapped to Python classes/objects. `extends` is mapped to Python inheritance.


### Categories

A `Category` is a special `Definition`.

- *constraint:* `Category` definition and its `categories` attribute do not form circles


### Packages
Packages are special definitions. Packages can contain section definitions
and category definitions.

A `Package` is a special `Definition` that contains definitions. `Packages` are mapped
to Python modules.

- *derived*: `definitions`, all definitions in this package
- *derived*: `sections`, all sections in this package
- *derived*: `categories`, all categories in this package


### Annotations

Arbitrary objects that can be attached to all definitions and contain additional information.


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

### MSection

`MSection` is a Python base-class for all sections and provides additional reflection.

- `m_def`: Python variable with the definition of this section
- `m_data`: container for all the section data
- `m_parent`: Python variable with the parent section instance
- `m_parent_index`: Python variable with the index in the parent's repeatable sub section
- `m_contents()`: all sub section instances
- `m_all_contents()`: traverse all sub and sub sub section instances
- `m_to_dict()`: serializable dict form
- `m_to_json()`


## Examples (of the Python interface)

### Definitions

This could be code, from a python module that represents the NOMAD *common* package `nomad.metainfo.common`:
```python
class System(MSection):
    '''
    The system is ...
    '''

    n_atoms = Quantity(type=int, derived_from='atom_labels')

    atom_labels = Quantity(
        shape=['n_atoms'],
        type=MEnum(ase.data.chemical_symbols),
        annotations=[ElasticSearchQuantity('keyword')])
    '''
    Atom labels are ...
    '''

    formula_hill = Quantity(type=str, derived_from=['atom_labels'])

    atom_species = Quantity(shape=['n_atoms'], type=int, derived_from='atom_labels')

    atom_positions = Quantity(shape=['n_atoms', 3], type=float, unit=units.m)

    cell = Quantity(shape=[3, 3], type=float, unit=units.m)
    lattice_vectors = Quantity(synonym='cell')

    pbc = Quantity(shape=[3], type=bool)

    # Not sure if this should be part of the definition. It will not serialize to
    # JSON. It might get complex for more involved cases. In many cases, we would
    # need both directions anyways. On the other hand, it allows to formally define
    # the derive semantics.
    def m_derive_atom_species(self) -> List[int]:
        return [ase.data.atomic_numbers[label] for label in self.atom_labels]

    def m_derive_n_atoms(self) -> int:
        return len(self.atom_labels)


class Run(MSection):

    systems = SubSection(System, repeats=True)
```

This could be part of the VASP source code:
```python
class Method(MSection):
    m_definition = Section(adds_to=nomad.metainfo.common.Method)

    incar_nbands = Quantity(
        type=int, links=['https://cms.mpi.univie.ac.at/wiki/index.php/NBANDS'])
```

### (Meta-)data

```python
from nomad.metainfo.common import Run, System

run = Run()

system = run.systems.create(atom_labels=['H', 'H', 'O'])
system.atom_positions = [[0, 0, 0], [1, 0, 0], [0.5, 0.5, 0]]
system.cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
system.pbc = [False, False, False]

print(system.atom_species)  # [1, 1, 96]
print(system.lattice_vectors)
print(system.n_atoms)

print(run.m_to_json(indent=2))
```

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

# Glossary

A list of words with very specific and precise meaning. This meaning might not yet be
fully expressed, but its there.

- annotation
- category
- derived quantity
- dimension
- new MetaInfo
- old MetaInfo
- package
- pythonic
- quantity
- reference
- section
- shape
- synonym
- unit
