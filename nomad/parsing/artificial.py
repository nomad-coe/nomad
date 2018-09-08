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

"""
Parser for creating artificial test, brenchmark, and demonstration data.
"""

from typing import Callable, IO, Any
import json
import os.path
import numpy as np

from nomadcore.local_meta_info import loadJsonFile, InfoKindEl

from nomad.parsing.backend import LocalBackend
from nomad.parsing.parser import Parser


class TemplateParser(Parser):
    """
    A parser that generates data based on a template given via the
    mainfile. The template is basically some archive json. Only
    """
    def __init__(self):
        super().__init__()
        # use vasp metainfo, not to really use it, but because it works
        file_dir = os.path.dirname(os.path.abspath(__file__))
        relative_metainfo_path = "../../.dependencies/nomad-meta-info/meta_info/nomad_meta_info/vasp.nomadmetainfo.json"
        meta_info_path = os.path.normpath(os.path.join(file_dir, relative_metainfo_path))
        self.meta_info_env, _ = loadJsonFile(filePath=meta_info_path, dependencyLoader=None, extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS, uri=None)
        self.name = 'parsers/template'
        self.backend = None

    def is_mainfile(self, filename: str, open: Callable[[str], IO[Any]]) -> bool:
        return filename.endswith('template.json')

    def add_section(self, section):
        if not isinstance(section, dict):
            print(section)

        name = section['_name']
        index = self.backend.openSection(name)

        for key, value in section.items():
            if key.startswith('x_') or key.startswith('_'):
                continue

            if key.startswith('section_'):
                values = value if isinstance(value, list) else [value]
                for value in values:
                    self.add_section(value)
            else:
                if isinstance(value, list):
                    shape = self.meta_info_env[key].get('shape')
                    if shape is None or len(shape) == 0:
                        for single_value in value:
                            self.backend.addValue(key, single_value, index)
                    else:
                        self.backend.addArrayValues(key, np.asarray(value), index)
                else:
                    self.backend.addValue(key, value, index)

        self.backend.closeSection(name, index)

    def run(self, mainfile: str) -> LocalBackend:
        self.backend = LocalBackend(metaInfoEnv=self.meta_info_env, debug=False)
        template_json = json.load(open(mainfile, 'r'))
        section = template_json['section_run'][0]
        self.add_section(section)
        self.backend.finishedParsingSession('ParseSuccess', [])
        return self.backend
