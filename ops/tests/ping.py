import requests
import time
from ase.data import chemical_symbols
import random
from datetime import datetime

base_url = 'https://nomad-lab.eu/prod/rae/api'

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
