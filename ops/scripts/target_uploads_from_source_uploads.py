import sys
from nomad import infrastructure
from nomad import config
from nomad.migration import Package
from pymongo import MongoClient
from nomad.processing import Upload

infrastructure.setup_logging()
infrastructure.setup_mongo()

sources = sys.argv[1:]
names = list(package.package_id for package in Package.objects(upload_id__in=sources))

config.mongo.db_name = 'fairdi_nomad_migration'
client = infrastructure.setup_mongo()

client = MongoClient()
db = client['fairdi_nomad_migration']
collection = db['upload']

targets = collection.find({'name': {'$in': names}})
for target in targets:
    print(target['_id'])
