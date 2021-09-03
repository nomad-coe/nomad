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

from nomad.metainfo import MCategory, Category


class FastAccess(MCategory):
    '''
    Used to mark archive objects that need to be stored in a fast 2nd-tier storage medium,
    because they are frequently accessed via archive API.

    If applied to a sub_section, the section will be added to the fast storage. Currently
    this only works for *root* sections that are sub_sections of `EntryArchive`.

    If applied to a reference types quantity, the referenced section will also be added to
    the fast storage, regardless if the referenced section has the category or not.
    '''

    m_def = Category()
