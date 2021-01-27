#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

'''
The NOMAD meta-info allows to define schemas for physics data independent of the used
storage format. It allows to define physics quantities with types, complex shapes
(vetors, matrices, etc.), units, links, and descriptions. It allows to organize large
amounts of these quantities in containment hierarchies of extendable sections, references
between sections, and additional quantity categories.

NOMAD uses the meta-info to define all archive data, repository meta-data, (and encyclopedia
data). The meta-info provides a convenient Python interface to create,
manipulate, and access data. We also use it to map data to various storage formats,
including JSON, (HDF5), mongodb, and elastic search.

Starting example
----------------

.. code-block:: python

    from nomad.metainfo import MSection, Quantity, SubSection, Units

    class System(MSection):
        \'\'\'
        A system section includes all quantities that describe a single a simulated
        system (a.k.a. geometry).
        \'\'\'

        n_atoms = Quantity(
            type=int, description=\'\'\'
            A Defines the number of atoms in the system.
            \'\'\')

        atom_labels = Quantity(type=MEnum(ase.data.chemical_symbols), shape['n_atoms'])
        atom_positions = Quantity(type=float, shape=['n_atoms', 3], unit=Units.m)
        simulation_cell = Quantity(type=float, shape=[3, 3], unit=Units.m)
        pbc = Quantity(type=bool, shape=[3])

    class Run(MSection):
        section_system = SubSection(sub_section=System, repeats=True)


We define simple metainfo schema with two `sections` called ``System`` and ``Run``. Sections
allow to organize related data into, well, `sections`. Each section can have two types of
properties: `quantities` and `sub-sections`. Sections and their properties are defined with
Python classes and their attributes.

Each `quantity` defines a piece of data. Basic quantity attributes are its `type`, `shape`,
`unit`, and `description`.

`Sub-sections` allow to place section into each other and therefore allow to form containment
hierarchies or sections and the respective data in them. Basic sub-section attributes are
`sub_section`(i.e. a reference to the section definition of the sub-section) and `repeats`
(determines if a sub-section can be contained once or multiple times).

The above simply defines a schema, to use the schema and create actual data, we have to
instantiate the above classes:

.. code-block:: python

    run = Run()
    system = run.m_create(System)
    system.n_atoms = 3
    system.atom_labels = ['H', 'H', 'O']

    print(system.atom_labels)
    print(n_atoms = 3)


Section `instances` can be used like regular Python objects: quantities and sub-sections
can be set and access like any other Python attribute. Special meta-info methods, starting
with ``m_`` allow us to realize more complex semantics. For example ``m_create`` will
instantiate a sub-section and add it to the `parent` section in one step.

Another example for an ``m_``-method is:

.. code-block:: python

    run.m_to_json(indent=2)

This will serialize the data into JSON:

.. code-block:: JSON

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

Definitions and Instances
-------------------------

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


Common attributes of Metainfo Definitions
-----------------------------------------

In the example, you already saw the basic Python interface to the Metainfo. *Sections* are
represented in Python as objects. To define a section, you write a Python classes that inherits
from ``MSection``. To define sub-sections and quantities you use Python properties. The
definitions themselves are also objects derived from classes. For sub-sections and
quantities, you directly instantiate :class:`SubSection` and :class`Quantity`. For sections
there is a generated object derived from :class:`Section` that is available via
`m_def` from each `section class` and `section instance`.

These Python classes that are used to represent metainfo definitions form an inheritance
hierarchy to share common properties

.. autoclass:: Definition


Quantities
----------

Quantity definitions are the main building block of meta-info schemas. Each quantity
represents a single piece of data.

.. autoclass:: Quantity


Sections and Sub-Sections
-------------------------

.. _metainfo-sections:


The NOMAD Metainfo allows to create hierarchical (meta-)data structures. A hierarchy
is basically a tree, where *sections* make up the root and inner nodes of the tree,
and *quantities* are the leaves of the tree. In this sense a section can hold further
sub-sections (more branches of the tree) and quantities (leaves). We say a section can
have two types of *properties*: *sub-sections* and *quantities*.

There is a clear distinction between *section* and *sub-section*. The term *section*
refers to an object that holds data (*properties*) and the term *sub-section* refers to a
relation between *sections*, i.e. between the containing section and the contained
sections of the same type. Furthermore, we have to distinguish between *sections* and
*section definitions*, as well as *sub-sections* and *sub-section definitions*. A *section
definition* defines the possible properties of its instances, and an instantiating *section*
can store quantities and sub-sections according to its definition. A *sub-section*
definition defines a principle containment relationship between two *section definitions*
and a *sub-section* is a set of *sections* (contents) contained in another *section* (container).
*Sub-section definitions* are parts (*properties*) of the definition of containing sections.

.. autoclass:: Section

.. autoclass:: SubSection

.. _metainfo-categories:


References and Proxies
----------------------


Beside creating hierarchies (e.g. tree structures) with ``SubSection``, the metainfo
also allows to create cross references between sections and other sections or quantity
values:

.. code-block:: python

    class Calculation(MSection):
        system = Quantity(type=System.m_def)
        atom_labels = Quantity(type=System.atom_labels)

    calc = Calculation()
    calc.system = run.systems[-1]
    calc.atom_labels = run.systems[-1]


To define a reference, you define a normal quantity and simply use the section or quantity
that you want to reference as type. Then you can assign respective section instances
as values.

When in Python memory, quantity values that reference other sections simply contain a
Python reference to the respective `section instance`. However, upon serializing/storing
metainfo data, these references have to be represented differently.

Value references are a little different. When you read a value references it behaves like
the references value. Internally, we do not store the values, but a reference to the
section that holds the referenced quantity. Therefore, when you want to
assign a value reference, you use the section with the quantity and not the value itself.

Currently this metainfo implementation only supports references within a single
section hierarchy (e.g. the same JSON file). References are stored as paths from the
root section, over sub-sections, to the referenced section or quantity value. Each path segment is
the name of the sub-section or an index in a repeatable sub-section:
``/system/0`` or ``/system/0/atom_labels``.

References are automatically serialized by :py:meth:`MSection.m_to_dict`. When de-serializing
data with :py:meth:`MSection.m_from_dict` these references are not resolved right away,
because the references section might not yet be available. Instead references are stored
as :class:`MProxy` instances. These objects are automatically replaced by the referenced
object when a respective quantity is accessed.

.. autoclass:: MProxy

If you want to defined references, it might not be possible to define the referenced
section or quantity before hand, due to how Python definitions and imports work. In these
cases, you can use a proxy to reference the reference type:

.. code-block:: python

    class Calculation(MSection):
        system = Quantity(type=MProxy('System')
        atom_labels = Quantity(type=MProxy('System/atom_labels')

The strings given to ``MProxy`` are paths within the available definitions. The above example
works, if ``System`` and ``System/atom_labels`` are eventually defined in the same package.


Categories
----------

In the old meta-info this was known as `abstract types`.

Categories are defined with Python classes that have :class:`MCategory` as base class.
Their name and description is taken from the class's name and docstring. An example
category looks like this:

.. code-block:: python

    class CategoryName(MCategory):
        \'\'\' Category description \'\'\'
        m_def = Category(links=['http://further.explanation.eu'], categories=[ParentCategory])

Packages
--------

.. autoclass:: Package

.. _metainfo-custom-types:


Environments
------------

.. autoclass:: Environment


Custom data types
-----------------

.. autoclass:: DataType
    :members:

.. autoclass:: MEnum


.. _metainfo-reflection:

Reflection and custom data storage
----------------------------------

When manipulating metainfo data in Python, all data is represented as Python objects, where
objects correspond to `section instance` and their attributes to `quantity values` or
`section instances` of sub-sections. By defining sections with `section classes` each
of these Python objects already has an interface that allows to get/set quantities and
sub-sections. But often this interface is too limited, or the specific section and
quantity definitions are unknown when writing code.

.. autoclass:: MSection
    :members:

.. autoclass:: MetainfoError
.. autoclass:: DeriveError
.. autoclass:: MetainfoReferenceError

.. _metainfo-urls:

Resources
---------

.. autoclass:: MResource

A more complex example
----------------------

.. literalinclude:: ../nomad/metainfo/example.py
    :language: python

'''


from .metainfo import (
    MSectionBound,
    MSection,
    MCategory,
    Definition,
    Property,
    Quantity,
    SubSection,
    Section,
    Category,
    Package,
    Environment,
    MEnum,
    Datetime,
    Capitalized,
    MProxy,
    MetainfoError,
    DeriveError,
    MetainfoReferenceError,
    DataType,
    Reference,
    QuantityReference,
    Datetime,
    Unit,
    JSON,
    Dimension,
    MResource,
    m_package,
    Annotation,
    DefinitionAnnotation,
    SectionAnnotation,
    SectionProxy,
    derived,
    constraint,
    units)
