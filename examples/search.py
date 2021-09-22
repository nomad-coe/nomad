from nomad import config, infrastructure
from nomad.search.v1 import search
from nomad.app.v1.models import MetadataPagination


config.elastic.host = 'localhost'
config.elastic.port = 19202
config.elastic.index_name = 'fairdi_nomad_prod_v0_8'

infrastructure.setup_elastic()

for entry in search(pagination=MetadataPagination(page_size=1000), query=dict(authors='Emre Ahmetcik')).data:
    print('entry_id', entry['entry_id'])
