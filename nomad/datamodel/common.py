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

import numpy as np

from nomad import config


def get_optional_backend_value(backend, key, section, unavailable_value=None, logger=None):
    # Section is section_system, section_symmetry, etc...
    val = None  # Initialize to None, so we can compare section values.
    # Loop over the sections with the name section in the backend.
    for section_index in backend.get_sections(section):
        if section == 'section_system':
            try:
                if not backend.get_value('is_representative', section_index):
                    continue
            except (KeyError, IndexError):
                continue

        try:
            new_val = backend.get_value(key, section_index)
        except (KeyError, IndexError):
            new_val = None

        # Compare values from iterations.
        if val is not None and new_val is not None:
            if val.__repr__() != new_val.__repr__() and logger:
                logger.warning(
                    'The values for %s differ between different %s: %s vs %s' %
                    (key, section, str(val), str(new_val)))

        val = new_val if new_val is not None else val

    if val is None and logger:
        logger.warning(
            'The values for %s where not available in any %s' % (key, section))
        return unavailable_value if unavailable_value is not None else config.services.unavailable_value
    else:
        if isinstance(val, np.generic):
            return val.item()

        return val
