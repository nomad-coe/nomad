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

from nomad.metainfo import MSection, Package, Section, Quantity

m_package = Package(name='material_library')


class Chemical(MSection):
    pass


class Sample(MSection):
    pass


class Instrument(MSection):
    pass


class Process(MSection):
    pass


class Measurement(MSection):
    pass


m_package.__init_metainfo__()
