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
A simple script that applies the output of `pip freeze` to the packages mentioned
in a requirements.txt while leaving the rest untouched.
'''

import re


freeze_path = 'freeze.txt'
requirements_path = 'requirements.txt'


with open(freeze_path, 'rt') as f:
    pkgs = {}
    for line in f.readlines():
        if '==' in line:
            key, value = line.split('==')
            pkgs[key] = value


with open(requirements_path, 'rt') as f:
    for line in f.readlines():
        match = re.search(r'^([a-zA-Z0-9_\-]+)(\[[^\]]+\])?', line)
        if match:
            pkg = match.group(1)
            if pkg in pkgs:
                print(f'{match.group(0)}=={pkgs[pkg]}', end='')
                continue

        print(line, end='')
