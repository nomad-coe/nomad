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

import re

from nomad import utils

# ../entries/<entry_id>/archive#<path>
# /entries/<entry_id>/archive#<path>
_regex_form_a = re.compile(r'^(?:\.\.)?/entries/([^?]+)/(archive|raw)#([^?]+?)$')

# ../upload/<upload_id>/archive/<entry_id>#<path>
# /uploads/<upload_id>/archive/<entry_id>#<path>
# <installation>/uploads/<upload_id>/archive/<entry_id>#<path>
_regex_form_b = re.compile(r'^([^?]+?)?/uploads?/([\w=-]*)/?(archive|raw)/([^?]+?)#([^?]+?)$')


def parse_path(url: str, upload_id: str = None):
    '''
    Parse a reference path.

    The upload_id of current upload is taken as the input to account for that the relative reference has no
    information about the upload_id, and it may also contain no entry_id. Has to know the upload_id when only
    path to mainfile is given.

    Returns:
        (installation, upload_id, entry_id_or_mainfile, kind, path): successfully parsed path
        None: fail to parse path
    '''

    url_match = _regex_form_b.match(url)
    if not url_match:
        # try another form
        url_match = _regex_form_a.match(url)
        if not url_match:
            # not valid
            return None

        entry_id_or_mainfile = url_match.group(1)
        kind = url_match.group(2)  # archive or raw
        path = url_match.group(3)
        if not path.startswith('/'):
            path = '/' + path

        return None, upload_id, entry_id_or_mainfile, kind, path

    installation = url_match.group(1)
    if installation == '':
        installation = None
    elif installation == '..':
        installation = None

    # if empty, it is a local reference to the same upload, use the current upload_id
    other_upload_id = upload_id if url_match.group(2) == '' else url_match.group(2)

    kind = url_match.group(3)  # archive or raw
    entry_id_or_mainfile = url_match.group(4)
    path = url_match.group(5)
    if not path.startswith('/'):
        path = '/' + path

    if kind == 'archive':
        if entry_id_or_mainfile.startswith('mainfile/'):
            entry_id_or_mainfile = utils.generate_entry_id(
                other_upload_id, entry_id_or_mainfile.replace('mainfile/', ''))
        elif '/' in entry_id_or_mainfile:  # should not contain '/' in entry_id
            return None

    return installation, other_upload_id, entry_id_or_mainfile, kind, path
