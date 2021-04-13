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
The archive storage is made from two tiers. First the whole archive is stored in
files, secondly parts of the archive are stored in mongodb documents.

The file storage is done in msg-pack files. Each file contains the archive of many
entries (all entries of an upload). These msg-pack files contain the JSON serialized
version of the metainfo archive (see module:`nomad.metainfo`). In addition msg-pack
contains TOC information for quicker access of individual sections. See :func:`write_archive`
and :func:`read_archvive`. In addition there query functionality to partially read
specified sections from an archive: func:`query_archive`.

The mongo storage uses mongodb's native bson to store JSON serialized metainfo archive
data. Each document in mongodb holds the partial archive of single entry. Which parts
of an archive are stored in mongo is determined by the metainfo and
section annotations/categories.
'''

from .storage import (
    write_archive, read_archive, ArchiveError, ArchiveReader, ArchiveWriter,
    ArchiveObject, ArchiveList, ArchiveItem)
from .query import query_archive, filter_archive, ArchiveQueryError
from .partial import (
    read_partial_archive_from_mongo, read_partial_archives_from_mongo,
    write_partial_archive_to_mongo, delete_partial_archives_from_mongo,
    create_partial_archive, compute_required_with_referenced)
from .required import RequiredReader, RequiredValidationError
