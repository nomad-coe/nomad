
import os
import sys
import importlib

# patch optimade python tools config and log handling
os.environ['OPTIMADE_CONFIG_FILE'] = os.path.join(os.path.dirname(__file__), 'optimade_config.json')
sys.modules['optimade.server.logger'] = importlib.import_module('nomad.app_fastapi.optimade_logger')

# patch optimade base path
from nomad import config  # nopep8
from optimade.server.config import CONFIG  # nopep8
CONFIG.root_path = "%s/optimade" % config.services.api_base_path

from optimade.server import main as optimade  # nopep8
from optimade.server.main import app as optimade_app  # nopep8
from optimade.server.routers import structures  # nopep8

# remove all the test data
from optimade.server.routers import ENTRY_COLLECTIONS  # nopep8

for name, collection in ENTRY_COLLECTIONS.items():
    if name == 'links':
        collection.collection.drop()
        collection.collection.insert_one({
            "id": "index",
            "type": "links",
            "name": "Index meta-database",
            "description": "Index for NOMAD databases",
            "base_url": "http://providers.optimade.org/index-metadbs/nmd",
            "homepage": "https://nomad-lab.eu",
            "link_type": "root"
        })
    else:
        collection.collection.drop()

# patch the structure collection with out elasticsearch implementation
from .elasticsearch import ElasticsearchStructureCollection  # nopep8
from .filterparser import parse_filter  # nopep8

structures.structures_coll = ElasticsearchStructureCollection()
optimade.add_major_version_base_url(optimade.app)
