# How to work with units

Units are a very important part of any scientific work. Units are also a common
source of problems when multiple people are working with the same data. Sometimes this has far-reaching consequences as demonstrated by the [Mars Climate Orbiter
incident](https://en.wikipedia.org/wiki/Mars_Climate_Orbiter#Cause_of_failure){:target="_blank"}. This document explains the possibilities and best practices when working with units within the NOMAD infrastructure.

## Available unit names

The available unit names are controlled by the `nomad/units/defaults_en.txt`
file. This is a plain text file that is internally read by the
[Pint](https://pint.readthedocs.io/en/stable/){:target="_blank"} library that we
use for unit transformations in the NOMAD Python backend. These definitions and the
associated conversion factors are then used always when transforming between two
units in the NOMAD platform. Our graphical user interface also performs unit
conversions in Javascript, and the same definitions and conversion factors are
translated to the frontend through an environment config file (`env.js`).

It is possible to add new units in this file to make them available both in the
Python and in the Javascript environment. You can read more on how units are
defined
[here](https://pint.readthedocs.io/en/stable/advanced/defining.html){:target="_blank"}.

All units support the use of SI-prefixes. This means that if the unit `meter`
has been defined, you are able to then automatically use `kilometers`,
`centimeters`, `millimeters`, etc.

## Defining units for storing data

From a practical point of view, there needs to be single choice for the unit in
which the data is stored on hard disks, databases or in search engines. This
choice is controlled by the `unit` attribute of a
`Quantity`:

```python
from nomad.metainfo import Quantity

my_energies = Quantity(
  dtype=float,
  unit='eV'
)
```

The data will always be stored in this unit and will be returned in this unit
when using e.g. the API.

## Using units in Python

When creating custom schemas or parsers in the NOMAD framework, it is important
to use the same unit definitions and conversion factors throughout. **You should
never define custom unit conversion routines of factors, but instead use the
`nomad.units` package**. This is important to ensure the interoperability of
data.

Here is an example of how you could work with units in Python:

```python
import numpy as np
from nomad.units import ureg  # Always import from here, never use another registry!

from nomad.metainfo import MSection, Quantity


# Here is a section with a quantity definition
class MySection(MSection):
    my_energies = Quantity(
        type=np.float64,
        shape=[2],
        unit='eV'
    )


my_section = MySection()

# If we assign a plain number array to a quantity, it is assumed to be given in
# the units defined in the Quantity
energies = np.array([1, 1])
my_section.my_energies = energies
print(my_section.my_energies)  # [1 1] electron_volt

# We can make a plain number array into a Pint Quantity by multiplying with
# a unit
energies_with_unit = energies * ureg('hartree')

# All numpy operations still work as normal, but the unit information is always
# stored alongside
energies_with_unit *= 10

# If you now assign this data into a NOMAD Quantity, the unit conversion is done
# automatically
my_section.my_energies = energies_with_unit
print(my_section.my_energies)  # [272.11386245988473 272.11386245988473] electron_volt
```

## Defining units for displaying data

When data is being displayed by the GUI, the unit can be choosen independently
from the unit used in storing the data. There are several reasons for doing
this:

 - Maybe the unit is stored in SI units for consistency, but when viewing the
 data you want to view it in some more field-specific units
 - Maybe the unit is stored in some field specific units, but when demonstrating
 your work you wish to use more standardized units.
 - Maybe due to your background, you are more familiar with a specific unit, and
 viewing the data in this unit helps you to understand it better. E.g. a
 physicist might be familiar with working with electron volts, whereas a chemist
 might prefer kilocalorie per mole.

Currently the display unit is controlled through the [ELN annotation](../../reference/annotations.md#eln-annotations), like this:

```python
distance = Quantity(
  dtype=float,
  unit='meter',
  a_eln=dict(defaultDisplayUnit='millimeter')
)
```