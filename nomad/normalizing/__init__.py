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
After parsing calculations have to be normalized with a set of *normalizers*.
In NOMAD-coe those were programmed in python (we'll reuse) and scala (we'll rewrite).

Currently the normalizers are:
- system.py (contains aspects of format stats, system, system type, and symmetry normalizer)
- optimade.py
- fhiaims.py
- dos.py
- encyclopedia.py (used to create the data in NOMAD-coe Encyclopedia)

The normalizers are available via

.. autodata:: nomad.normalizing.normalizers

There is one ABC for all normalizer:

.. autoclass::nomad.normalizing.normalizer.Normalizer
    :members:
'''

from typing import List, Any, Iterable, Type

from .system import SystemNormalizer
from .optimade import OptimadeNormalizer
from .fhiaims import FhiAimsBaseNormalizer
from .dos import DosNormalizer
from .normalizer import Normalizer
from .band_structure import BandStructureNormalizer
from .encyclopedia.encyclopedia import EncyclopediaNormalizer
from .workflow import WorkflowNormalizer

normalizers: Iterable[Type[Normalizer]] = [
    SystemNormalizer,
    OptimadeNormalizer,
    # FhiAimsBaseNormalizer,
    DosNormalizer,
    BandStructureNormalizer,
    EncyclopediaNormalizer,
    WorkflowNormalizer,
]
