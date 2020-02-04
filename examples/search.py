from nomad import config, infrastructure, search


config.elastic.host = 'localhost'
config.elastic.port = 19202
config.elastic.index_name = 'fairdi_nomad_prod_v0_7'

infrastructure.setup_logging()
infrastructure.setup_elastic()


req = search.SearchRequest()
req.search_parameter('authors', 'Emre Ahmetcik')
upload_id = None

i = 0
for entry in req.execute_scan(order_by='upload_id', size=1000):
    i += 1
    if entry['upload_id'] != upload_id:
        upload_id = entry['upload_id']
        print(upload_id, i)
