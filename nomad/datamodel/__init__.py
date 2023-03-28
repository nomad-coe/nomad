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

import sys

from nomad.metainfo import Environment
from .data import (User, Author, user_reference, author_reference)
from .datamodel import (
    Dataset, EditableUserMetadata, AuthLevel,
    MongoUploadMetadata, MongoEntryMetadata, MongoSystemMetadata,
    EntryMetadata, EntryArchive)
from .optimade import OptimadeEntry, Species
from .metainfo import m_env
from .results import Results
from .data import EntryData, ArchiveSection
from .context import Context, ClientContext, ServerContext

m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.datamodel'].m_package)  # type: ignore
m_env.m_add_sub_section(Environment.packages, sys.modules['nomad.datamodel.optimade'].m_package)  # type: ignore
