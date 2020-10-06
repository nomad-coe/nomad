from locust import HttpUser, task, between
from ase.data import chemical_symbols
import random


# These are the API requests from the search UI with various tabs and statistics
query_params = [
    'page=1&per_page=10&order_by=upload_time&order=-1&domain=dft&owner=public&atoms=Co&statistics=atoms&exclude=atoms,only_atoms,dft.files,dft.quantities,dft.optimade,dft.labels,dft.geometries',
    'page=1&per_page=10&order_by=upload_time&order=-1&domain=dft&owner=public&statistics=dft.labels_springer_compound_class&statistics=dft.system&statistics=dft.crystal_system&statistics=dft.compound_type&exclude=atoms,only_atoms,dft.files,dft.quantities,dft.optimade,dft.labels,dft.geometries',
    'page=1&per_page=10&order_by=upload_time&order=-1&domain=dft&owner=public&statistics=dft.code_name&statistics=dft.basis_set&statistics=dft.xc_functional&exclude=atoms,only_atoms,dft.files,dft.quantities,dft.optimade,dft.labels,dft.geometries',
    'page=1&per_page=10&order_by=upload_time&order=-1&domain=dft&owner=public&statistics=dft.searchable_quantities&statistics=dft.labels_springer_classification&statistics=dft.workflow.workflow_type&exclude=atoms,only_atoms,dft.files,dft.quantities,dft.optimade,dft.labels,dft.geometries',
    'page=1&per_page=10&order_by=upload_time&order=-1&domain=dft&owner=public&statistics=dft.searchable_quantities&statistics=dft.labels_springer_classification&statistics=dft.workflow.workflow_type&exclude=atoms,only_atoms,dft.files,dft.quantities,dft.optimade,dft.labels,dft.geometries&datasets_grouped=true',
    'page=1&per_page=10&order_by=upload_time&order=-1&domain=dft&owner=public&metrics=dft.calculations&statistics=atoms&exclude=atoms,only_atoms,dft.files,dft.quantities,dft.optimade,dft.labels,dft.geometries&datasets_grouped=true'
]


class QuickstartUser(HttpUser):
    wait_time = between(1, 2)

    @task
    def empty_search(self):
        self.client.get('/prod/rae/beta/api/repo/?%s' % query_params[0])

    @task(3)
    def elements_search(self):
        self.client.get("/prod/rae/beta/api/repo/?%s&atoms=%s" % (
            random.choice(query_params),
            random.choice(chemical_symbols[1:])))

    def on_start(self):
        pass
