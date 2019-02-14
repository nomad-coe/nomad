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


import matid.symmetry.symmetryanalyzer
import matid.utils.segfault_protect


# A patch for the segfault protection of systax (internally uses protection for spglib calls.)
# We basically disable the protection. The multiprocessing based original protection.
# somehow interfers with the celery work infrastructure and leads to a deadlock. Its a TODO.
# It also seems to deadlock without celery .. just not working consistently.
def segfault_protect_patch(f, *args, **kwargs):
    return f(*args, **kwargs)


matid.symmetry.symmetryanalyzer.segfault_protect = segfault_protect_patch
matid.utils.segfault_protect.segfault_protect = segfault_protect_patch
