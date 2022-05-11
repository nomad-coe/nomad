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
Definitions that are used in the documentation via mkdocs-macro-plugin.
'''

from nomad.app.v1.models import (
    query_documentation,
    owner_documentation)
from nomad.app.v1.routers.entries import archive_required_documentation
from nomad import utils

doc_snippets = {
    'query': query_documentation,
    'owner': owner_documentation,
    'archive-required': archive_required_documentation
}


def define_env(env):
    @env.macro
    def nomad_url():  # pylint: disable=unused-variable
        # TODO Fix the configuration during build time.
        return 'https://nomad-lab.eu/prod/v1/staging/api'
        # return config.api_url()

    @env.macro
    def doc_snippet(key):  # pylint: disable=unused-variable
        return doc_snippets[key]

    @env.macro
    def metainfo_data():  # pylint: disable=unused-variable
        return utils.strip('''
            You can browse the [NOMAD metainfo schema](../gui/analyze/metainfo)
            or the archive of each entry (e.g. [a VASP example](../gui/search/entries/entry/id/d5OYC0SJTDevHMPk7YHd4A/-7j8ojKkna2NLXdytv_OjV4zsBXw/archive))
            in the web-interface.''')
