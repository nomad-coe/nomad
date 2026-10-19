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

from pydantic.v1.error_wrappers import ErrorWrapper


class CustomErrorWrapper(ErrorWrapper):
    """ErrorWrapper that supports dict conversion via dict(instance)."""

    def __iter__(self):
        """Yields (key, value) pairs for dict conversion: error and location."""
        yield 'msg', str(self.exc)
        yield 'loc', self.loc_tuple()

    def __getitem__(self, key):
        """Enables dictionary-style access: error['loc'] and error['msg']."""
        if key == 'loc':
            return self.loc_tuple()
        elif key == 'msg':
            return str(self.exc)
        raise KeyError(key)
