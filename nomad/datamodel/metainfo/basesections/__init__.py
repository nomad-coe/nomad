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

# --- For Backward Compatibility ---
# Make all v1 contents available at the top level of the `basesections` package.
# Users can continue to do: from nomad.datamodel.metainfo.basesections import SomeV1Class
from nomad.datamodel.metainfo.basesections.v1 import *

# --- For Forward Compatibility and Clarity ---
# Import the v2 module itself, so it's available as a namespace.
# Users can explicitly access v2 content like this:
# from nomad.datamodel.metainfo.basesections import v2
from nomad.datamodel.metainfo.basesections import v2
