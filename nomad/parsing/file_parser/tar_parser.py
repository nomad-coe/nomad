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

import tarfile

from .file_parser import FileParser


class TarParser(FileParser):
    def __init__(self, mainfile=None, logger=None):
        super().__init__(mainfile, logger, tarfile.open)
        self._names_map = None

    @property
    def names_map(self):
        if self._names_map is None:
            if self.mainfile_obj is not None:
                self._names_map = {f.lower(): f for f in self.mainfile_obj.getnames()}
        return self._names_map

    def parse(self, key):
        if self._results is None:
            self._results = dict()

        if self.mainfile_obj is None:
            return

        name = self.names_map.get(key, key)
        try:
            val = self.mainfile_obj.extractfile(name)
        except Exception:
            val = None

        self._results[key] = val
