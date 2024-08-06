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

import os
import pkgutil


def get_package_path(package_name: str) -> str:
    """Given a python package name, returns the filepath of the package root folder."""
    package_path = None
    try:
        # We try to deduce the package path from the top-level package
        package_path_segments = package_name.split('.')
        root_package = package_path_segments[0]
        package_dirs = package_path_segments[1:]
        package_path = os.path.join(
            os.path.dirname(
                pkgutil.get_loader(root_package).get_filename()  # type: ignore
            ),
            *package_dirs,
        )
        if not os.path.isdir(package_path):
            # We could not find it this way. Let's try to official way
            package_path = os.path.dirname(
                pkgutil.get_loader(package_name).get_filename()  # type: ignore
            )
    except Exception as e:
        raise ValueError(f'The python package {package_name} cannot be loaded.', e)

    return package_path
