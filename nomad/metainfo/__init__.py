# Copyright 2018 Markus Scheidgen
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an"AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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
        systems = SubSection(sub_section=System, repeats=True)


We define simple metainfo schema with two `sections` called ``System`` and ``Run``. Sections
allow to organize related data into,  well, `sections`. Each section can have two types of
properties: `quantities` and `sub-sections`. Sections and their properties are defined with
Python classes and their attributes.

Each `quantity` defines a piece of data. Basic quantity attributes are its `type`, `shape`,
`unit`, and `description`.

`Sub-sections` allow to place section into each other and there allow to form containment
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

Definitions
-----------

.. autoclass:: Definition


Quantities
----------

.. autoclass:: Quantity

.. _metainfo-sections:

Sections
--------

With sections it is paramount to always be clear what is talked about. The lose
term `section` can reference one of the following three:

* `section definition`
    Which is a Python object that represents the definition of
    a section, its sub-sections and quantities. `Section definitions` should not be not
    written directly. `Section definitions` are objects of :class:`Section`.
* `secton class`
    Which is a Python class and :class:`MSection` decendant that is
    used to express a `section defintion` in Python. Each `section class` is tightly
    associated with its `section definition`. The `section definition` can be access
    with the class attribute ``m_def``. The `section definition` is automatically created
    from the `section class` upon defining the class through metaclass vodoo.
* `section instance`
    The instance (object) of a `section class`, it `follows` the
    definition associated with the instantiated `section class`. The followed
    section definition can be accessed with the object attribute ``m_def``.

A `section class` looks like this:

.. code-block:: python

    class SectionName(BaseSection):
        \'\'\' Section description \'\'\'
        m_def = Section(**section_attributes)

        quantity_name = Quantity(**quantity_attributes)
        sub_section_name = SubSection(**sub_section_attributes)

The various Python elements of this class are mapped to a respective *section definition*.
The ``SectionName`` becomes the *name*. The
``BaseSection`` is either :class:`MSection` or another *section class*. The
``section_attributes`` become additional attributes of the `section definition`. The
various ``Quantity`` and ``SubSection`` become the *quantities* and *sub sections*.

Each *section class* has to directly or indirectly extend :class:`MSection`. This will
provided certain class and object features to all *section classes* and all *section instances*.
Read :ref:metainfo-reflection to learn more.

.. autoclass:: Section

Sub-Sections
------------

.. autoclass:: SubSection

.. _metainfo-categories:

Categories
----------

.. autoclass:: Quantity

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

References and metainfo URLs
----------------------------

When in Python memory, quantity values that reference other sections simply contain a
Python reference to the respective `section instance`. However, upon serializing/storing
metainfo data, these references have to be represented differently.

Currently this metainfo implementation only supports references within a single
section hierarchy (e.g. the same JSON file). References are stored as paths from the
root section, over sub-sections, to the references section. Each path segment is
the name of the sub-section or an index in a repeatable sub-section:
``/system/0/symmetry``.

References are automatically serialized by :py:meth:`MSection.m_to_dict`. When de-serializing
data with :py:meth:`MSection.m_from_dict` these references are not resolved right away,
because the references section might not yet be available. Instead references are stored
as :class:`MProxy` instances. These objects are automatically replaced by the referenced
object when a respective quantity is accessed.

.. autoclass:: MProxy

Resources
_________

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
    MProxy,
    MetainfoError,
    DeriveError,
    MetainfoReferenceError,
    DataType,
    Reference,
    MResource,
    m_package,
    Annotation,
    DefinitionAnnotation,
    SectionAnnotation,
    SectionProxy,
    derived,
    constraint,
    units)
