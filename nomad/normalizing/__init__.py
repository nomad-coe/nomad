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
from typing import List, Any

from .normalizer import Normalizer
from .system import SystemNormalizer
from .symmetry import SymmetryAndTypeNormalizer
from .systemtype import SystemTypeNormalizer
from .fhiaims import FhiAimsBaseNormalizer

"""
After parsing calculations have to be normalized with a set of *normalizers*.
In NOMAD-coe those were programmed in python (we'll reuse) and scala (we'll rewrite).

Currently the normalizers are:
- system.py
- symmetry.py
- fhiaims.py
- systemtype.py

The normalizers are available via

.. autodata:: nomad.normalizing.normalizers

There is one ABC for all normalizer:

.. autoclass::nomad.normalizing.normalizer.Normalizer
    :members:
"""

# The loose explicit type is necessary to avoid a ABC class as item type that causes
# further errors on instantiating the normalizers. A solution would be to use objects
# instead of classes.
normalizers: List[Any] = [
    SystemNormalizer,
    FhiAimsBaseNormalizer,
]
