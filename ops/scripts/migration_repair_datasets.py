from nomad import infrastructure
from dateutil.parser import parse
import datetime

infrastructure.setup_logging()
mongo = infrastructure.setup_mongo()
calcs = mongo.fairdi_nomad_prod_v0_7.calc
datasets = mongo.fairdi_nomad_prod_v0_7.dataset


def doit(dataset_id):
    example = calcs.find_one({'metadata.datasets': dataset_id})
    if example is None:
        print('no example for %s' % dataset_id)
        return

    user_id = example['metadata']['uploader']

    if 'upload_time' not in example['metadata']:
        print('no upload time in %s' % dataset_id)
        upload_time = datetime.datetime.now()
    else:
        upload_time = example['metadata']['upload_time']

    update = [ 
        {'_id': dataset_id}, 
        {'$set': {
            'user_id': str(user_id),
            'created': parse(str(upload_time))
        }}
    ]
    
    datasets.update_one(*update)


for item in datasets.distinct('_id'):
    if 'user_id' not in item and 'created' not in item:
        doit(item)
