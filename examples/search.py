from nomad import infrastructure
from nomad.config import config
from nomad.search import search
from nomad.app.v1.models import MetadataPagination


config.elastic.host = 'localhost'
config.elastic.port = 19202
config.elastic.entries_index = 'fairdi_nomad_prod_v0_8'

infrastructure.setup_elastic()

for entry in search(
    pagination=MetadataPagination(page_size=1000), query=dict(authors='Emre Ahmetcik')
).data:
    print('entry_id', entry['entry_id'])
