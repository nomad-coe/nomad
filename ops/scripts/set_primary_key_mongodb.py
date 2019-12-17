import pymongo

from nomad import config

client = pymongo.MongoClient(config.mongo.host, int(config.mongo.port))
db = client[config.mongo.db_name]
collection = db.dataset

new_collection = []
for obj in collection.find({}):
    obj['_id'] = obj.pop('dataset_id')
    new_collection.append(obj)


collection.insert_many(new_collection)
collection.delete_many({'dataset_id': {'$exists': True}})
