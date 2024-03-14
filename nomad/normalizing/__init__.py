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

"""
After parsing entries have to be normalized with a set of *normalizers*.
In NOMAD-coe those were programmed in python (we'll reuse) and scala (we'll rewrite).

Currently the normalizers are:
- system.py (contains aspects of format stats, system, system type, and symmetry normalizer)
- optimade.py
- fhiaims.py
- dos.py

The normalizers are available via

.. autodata:: nomad.normalizing.normalizers

There is one ABC for all normalizer:

.. autoclass::nomad.normalizing.normalizer.Normalizer
    :members:
"""

import importlib
from typing import Any, Iterator
from collections import UserList
from .normalizer import Normalizer

from nomad import config


def import_normalizer(class_name: str):
    try:
        package, cls = class_name.rsplit('.', 1)
        return getattr(importlib.import_module(package), cls)
    except Exception as e:
        raise ImportError(f'Cannot import normalizer {class_name}', e)


class SortedNormalizers(UserList):
    def __iter__(self) -> Iterator:
        self.sort(key=lambda x: x.normalizer_level)
        return super().__iter__()


class NormalizerInterface:
    def __init__(self, path: str) -> None:
        self._path = path
        self._cls = None

    @property
    def normalizer_class(self):
        if self._cls is None:
            self._cls = import_normalizer(self._path)
        return self._cls

    def normalize(self, logger=None):
        self.normalizer_class.normalize(logger)

    def __call__(self, *args: Any) -> Any:
        return self.normalizer_class(*args)

    def __getattr__(self, name: str):
        return getattr(self.normalizer_class, name, None)


normalizers = SortedNormalizers([])


for plugin_name, plugin in config.plugins.options.items():
    if isinstance(plugin, config.Normalizer) and config.plugins.filter(plugin_name):
        normalizers.append(NormalizerInterface(plugin.normalizer_class_name))

for normalizer in config.normalize.normalizers.filtered_values():
    normalizers.append(NormalizerInterface(normalizer))
