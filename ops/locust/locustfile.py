from locust import HttpUser, task
from ase.data import chemical_symbols
import random


class HelloWorldUser(HttpUser):
    @task
    def prod(self):
        self.client.get('http://cloud.nomad-lab.eu/prod/v1/alive')

    @task
    def homepage(self):
        self.client.get('http://cloud.nomad-lab.eu/nomad-lab/index.html')

    @task
    def search(self):
        url = 'http://cloud.nomad-lab.eu/prod/v1/api/v1/entries?owner=public&include=entry_id'
        element = random.choice(chemical_symbols)
        url += f'&q=results.material.elements__{element}'
        self.client.request_name = (
            '/prod/v1/api/v1/entries?q=results.material.elements__X'
        )
        self.client.get(url)

    @task
    def search_sort(self):
        url = 'http://cloud.nomad-lab.eu/prod/v1/api/v1/entries?owner=public&include=entry_id&order_by=upload_create_time&order=desc'
        element = random.choice(chemical_symbols)
        url += f'&q=results.material.elements__{element}'
        self.client.request_name = '/prod/v1/api/v1/entries?q=results.material.elements__X&order_by=upload_create_time'
        self.client.get(url)

    @task
    def file(self):
        url = 'https://cloud.nomad-lab.eu/prod/v1/api/v1/uploads/Hy3Otg19QAesKS_Muv9GLw/raw/archive/S-edge/RIs_0.375ML/S-edge-CO%2BOH%2B3H.vasprun.xml?length=16384&decompress=true&ignore_mime_type=true'
        self.client.request_name = '/prod/v1/api/v1/uploads/<id>/raw/archive/<path>'
        self.client.get(url)

    @task
    def archive(self):
        url = 'https://cloud.nomad-lab.eu/prod/v1/api/v1/entries/xQiOo1umiJlIe_ZlTzoVnL-yFt3F/archive'
        self.client.request_name = '/prod/v1/api/v1/entries/<id>/archive'
        self.client.get(url)
