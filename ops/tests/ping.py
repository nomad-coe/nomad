#
# Copyright (c) 2018-2020 The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
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

import requests
import sys
import time
from ase.data import chemical_symbols
import random
from datetime import datetime

base_url = 'https://nomad-lab.eu/prod/rae/api'
if len(sys.argv) > 1:
    base_url = sys.argv[1]

while True:
    try:
        start = time.time()
        atoms = '&atoms=%s&atoms%s' % (random.choice(chemical_symbols), random.choice(chemical_symbols))
        response = requests.get('%s%s%s%s' % (
            base_url,
            '/repo/',
            '?page=1&per_page=10&order_by=upload_time&order=-1&domain=dft&owner=public&statistics=atoms&exclude=atoms,only_atoms,dft.files,dft.quantities,dft.optimade,dft.labels,dft.geometries',
            atoms))
        end = time.time()
        print('PING – %s – %f - %s' % (response.status_code, end - start, datetime.now()))
        time.sleep(1)
    except Exception as e:
        print('ERROR – %s' % e)
