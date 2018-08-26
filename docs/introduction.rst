Introduction
============

Nomad xt (and NOMAD coe) is about software for storing and processing (computational)
material science data based on *big-data* and *disributed computing*.

Architecture
------------

The following depicts the *nomad-xt* architecture with respect to software compenent
in terms of python modules, gui components, and 3rd party services (e.g. databases,
search enginines, etc.)

.. figure:: components.png
   :alt: nomad xt components

   The main components of nomad xt and their dependencies.

Principles
----------

* simple first, complicated only when necessary
* adopting generic established 3rd party solutions before implementing specific solutions
* only uni directional dependencies between components/modules, no circles
* only one language: Python (except, GUI of course)