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
from typing import Tuple, Generator, cast

from nomad.metainfo.metainfo import Quantity, Reference
from nomad.metainfo.search_extension import Search
from nomad.search import search_quantities


def provider_specific_fields() -> Generator[Tuple[str, Search], None, None]:
    for search_quantity in search_quantities.values():
        quantity = cast(Quantity, search_quantity.definition)
        if isinstance(quantity.type, Reference):
            # we can't yet support those
            continue

        nmd_name = search_quantity.qualified_name
        nmd_name_split = nmd_name.split('.')

        if len(nmd_name_split) == 1:
            # plain metadata
            pass
        elif not nmd_name_split[0] in ['dft', 'encyclopedia']:
            # other domains fields that do not make sense in the optimade context
            continue
        elif len(nmd_name_split) > 2 and nmd_name_split[1] == 'optimade':
            # these are already in optimade
            continue

        opt_name = nmd_name.replace('.', '_')
        yield opt_name, search_quantity
