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
'''

from typing import List, Type
import importlib

from nomad import config

from .normalizer import Normalizer


normalizers: List[Type[Normalizer]] = []

def add_normalizer(class_name: str):
    try:
        package, cls = class_name.rsplit('.', 1)
        normalizer = getattr(importlib.import_module(package), cls)
        normalizers.append(normalizer)
    except Exception as e:
        raise ImportError(f'Cannot import normalizer {class_name}', e)


for plugin_name, plugin in config.plugins.options.items():
    if isinstance(plugin, config.Normalizer) and config.plugins.filter(plugin_name):
        add_normalizer(plugin.normalizer_class_name)

for normalizer in config.normalize.normalizers.filtered_values():
    add_normalizer(normalizer)
