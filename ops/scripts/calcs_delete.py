# Simple python script template for deleting calculations from mongo and elastic
# based on elastic query

import elasticsearch_dsl as es

from nomad import infrastructure, config, processing

infrastructure.setup_logging()
infrastructure.setup_mongo()
infrastructure.setup_elastic()

query = es.Q('match', code_name='Octopus') \
    & ~es.Q('exists', field='pid') \
    & es.Q('wildcard', mainfile='*/inp') \
    & es.Q('match', published=True)

search = es.Search(index=config.elastic.index_name).query(query)

calc_ids = [hit.calc_id for hit in search.scan()]

input('Will delete %d calcs. Press any key to continue ...' % len(calc_ids))


def chunks(l, n):
    for i in range(0, len(l), n):
        print(i)
        yield l[i:i + n]


for chunk in chunks(calc_ids, 1000):
    processing.Calc.objects(calc_id__in=chunk).delete()

search.delete()
