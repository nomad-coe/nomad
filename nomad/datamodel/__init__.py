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

"""
This module contains classes that allow to represent the core
nomad data entities :class:`Upload` and :class:`Calc` on a high level of abstraction
independent from their representation in the different modules
:py:mod:`nomad.processing`, :py:mod:`nomad.coe_repo`, :py:mod:`nomad.parsing`,
:py:mod:`nomad.search`, :py:mod:`nomad.api`, :py:mod:`nomad.migration`.

It is not about representing every detail, but those parts that are directly involved in
api, processing, migration, mirroring, or other 'infrastructure' operations.

Transformations between different implementations of the same entity can be build
and used. To ease the number of necessary transformations the classes
:class:`UploadWithMetadata` and :class:`CalcWithMetadata` can act as intermediate
representations. Therefore, implement only transformation from and to these
classes. These are the implemented transformations:

.. image:: datamodel_transformations.png

.. autoclass:: nomad.datamodel.UploadWithMetadata
    :members:
.. autoclass:: nomad.datamodel.CalcWithMetadata
    :members:

The class :class:`CalcWithMetadata` only defines non domain specific metadata quantities
about ids, user metadata, etc. To define domain specific quantities the classes
:class:`Domain` and :class:`DomainQuantity` must be used.

.. autoclass:: nomad.datamodel.Domain
    :members:
.. autoclass:: nomad.datamodel.DomainQuantity
    :members:
"""

import sys

from nomad.datamodel.base import UploadWithMetadata, CalcWithMetadata, Domain
from nomad.datamodel import ems, dft
from nomad.datamodel.dft import DFTCalcWithMetadata
from nomad.datamodel.ems import EMSEntryWithMetadata

# Override the CalcWithMetadata with the domain specific decendant
setattr(sys.modules['nomad.datamodel'], 'CalcWithMetadata', Domain.instance.domain_entry_class)
