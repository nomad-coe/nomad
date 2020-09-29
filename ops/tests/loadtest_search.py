from locust import HttpUser, task, between
from ase.data import chemical_symbols
import random


class QuickstartUser(HttpUser):
    wait_time = between(1, 2)

    @task
    def empty_search(self):
        self.client.get("/prod/rae/beta/api/repo/")

    @task
    def elements_search(self):
        self.client.get("/prod/rae/beta/api/repo/?atoms=%s" % random.choice(chemical_symbols[1:]))

    def on_start(self):
        pass
