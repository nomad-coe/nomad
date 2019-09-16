# NOMAD MetaInfo

## History

The NOMAD MetaInfo was devised within the first NOMAD CoE; over 2000 quantities have
been defined in this *old MetaInfo*. The experience with this system revealed the following drawbacks:

- The Python libraries that allow to use the MetaInfo are non pythonic and incomplete.
- The MetaInfo is only used for the archive, not for the encyclopedia and repository data.
- There is direct support to map MetaInfo definitions to DI technologies (databases, search indices, APIs).
- There is no support for namespaces. MetaInfo names are cumbersome. This will not scale
to expected levels of FAIRmat metadata.
- MetaInfo packages are not version controlled. They are part of the same git and do not belong to the independently evolving parsers. This does not allow for "external" parser development and makes it hard to keep versions consistent.
- The MetaInfo is defined in JSON. The syntax is inadequate, checks are not immediate.

Attempts to revise the MetaInfo have failed in the past.

## Goals

### Common language to define physics (meta-)data quantities and their relationships

The *physics quantities* part includes
- each quantity has a physics *unit*
- *shapes* that precisely define vectors, matrices, tensors, and their dimensions

The *relationship* parts entails:
- hierarchies for quantity *values* (e.g. *sections*)
- hierarchies for quantity *definition* (e.g. *categories*, former *abstract types*)
- *derived* quantities that can be computed from other quantities
- *synonyms* as a special trivial case for derived quantities
- *shapes* might also define a type of relationship through one quantity being the dimension of another

In addition there are the *typical* data-type definition (schema, ontology, ...) features:
- names/namespaces
- modularization (i.e. Metainfo packages)
- documentation
- basic primitive types (int, string, bool)
- simple compound types (lists, dictionaries, unions)
- references between data objects

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

### Definition

`Definition` is the abstract base for all definitions in the MetaInfo.

### Attributes:
- `name`, a string
- `description`, a string
- `links`, a list of URLs
- `annotations`, a list of `Annotations`

### Property

`Property` is a special `Definition` and an abstract base for section properties.

#### Attributes
- `section`, a reference to a section definition

### Quantities (incl. dimensions)

A `Quantity` definition is a special and concrete `Property` definition:

#### Attributes
- `shape`, a list of either `int`, references to a dimension (quantity definition), or limits definitions (e.g. `'1..n'`, `'0..n'`.)
- `type`, a primitive or Enum type
- `categories`, a list of references to category definitions
- `section`, a reference to the parent section definition
- `unit`, a (computed) units, e.g. `units.F * units.m`
- `derived_from`, a list of references to other quantity definitions
- `synonym`, a reference to another quantity definition

A `Quantity`s are mapped to Python *descriptors*. *Dimensions* are quantity definitions
with empty shape and int type.

#### Contrains
- `synonym`, `derived_from`, and dimensions come from the same section

### Sections (incl. references)

A `Section` is a special and concrete `Definition`.

#### Attributes
- `adds_to`, a reference to another section definition. All quantities of this *pseudo* section are added to the given section.
- `parent_section`, a reference to another section definition
- `repeats`, a boolean
- `extends`, list of reference to other section definitions. This section automatically inherits all quantities of the other sections. (Might not be necessary)

#### Contrains
- `parent_section` is not circular
- `extends` is not circular
- `adds_to` is not circular
- all quantities that have *this* (or an `extends`) section as `section` have unique names

`Section`s are mapped to Python classes/objects.

### Categories

A `Category` is a special `Definition`.

#### Attributes
- `super_categories`, a list of references to other category definitions

#### Contrains
- `super_categories` is not circular

### Packages

A `Package` is a special `Definition` that contains definitions. `Packages` are mapped
to Python modules.

### References

A `Reference` is a special `Property`.

#### Attributes
- `referenced_section`, reference to a section definition

### Annotations

Arbitrary serializable objects that can contain additional information.


## Examples (of the Python interface)

### Definitions

This could be code, from a python module that represents the NOMAD *common* package `nomad.metainfo.common`:
```python
class System(MetainfoObject):
    """
    The system is ...
    """
    m_definition = Section(parent_section=Run, repeats=True)

    n_atoms = Quantity(type=int, derived_from='atom_labels')

    atom_labels = Quantity(
        shape=['n_atoms'],
        type=Enum(ase.data.chemical_symbols),
        annotations=[ElasticSearchQuantity('keyword')])
    """
    Atom labels are ...
    """

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
```

This could be part of the VASP source code:
```python
class Method(MetainfoObject):
    m_definition = Section(adds_to=nomad.metainfo.common.Method)

    incar_nbands = Quantity(
        type=int, links=['https://cms.mpi.univie.ac.at/wiki/index.php/NBANDS'])
```

### (Meta-)data

```python
from nomad.metainfo.common import Run, System

run = Run()

system = run.m_create(System)
system.atom_labels = ['H', 'H', 'O']
system.atom_positions = [[0, 0, 0], [1, 0, 0], [0.5, 0.5, 0]]
system.cell = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
system.pbc = [False, False, False]

print(system.atom_species)  # [1, 1, 96]
print(system.lattice_vectors)
print(system.n_atoms)
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
