'''
This examplefies how to send raw queries to elasticsearch.
'''
from typing import Iterable, Tuple, Set, Any
from concurrent.futures import ThreadPoolExecutor
from threading import Event, Lock
from queue import Queue, Empty
import json

from nomad import search, infrastructure, files, config

infrastructure.setup_files()
infrastructure.setup_elastic()


results = []
after = None
while True:
    request = {
        "query": {
            "bool": {
                "must": [
                    {
                        "match": {
                            "dft.quantities": "section_dos_fingerprint"
                        },
                    },
                    {
                        "match": {
                            "published": True
                        },
                    },
                    {
                        "match": {
                            "with_embargo": False
                        }
                    }
                ]
            }
        },
        "size": 0,
        "aggs": {
            "results": {
                "composite": {
                    "sources": {
                        "materials": {
                            "terms": {
                                "field": "encyclopedia.material.material_id"
                            }
                        }
                    },
                    "size": 10000
                },
                "aggs": {
                    "calcs": {
                        "top_hits": {
                            "sort": {
                                "_script": {
                                    "type": "number",
                                    "script": {
                                        "lang": "painless",
                                        "source": '''
                                            int result = 0;
                                            String code = doc['dft.code_name'].value;
                                            String functional = doc['dft.xc_functional'].value;
                                            if (functional == 'GGA') result += 100;
                                            if (code == 'VASP')
                                                result += 1;
                                            else if (code == 'FHI-aims')
                                                result += 2;
                                            return result;
                                        '''
                                    },
                                    "order": "asc"
                                },
                            },
                            "_source": {
                                "includes": ['upload_id', 'calc_id', 'dft.code_name', 'dft.xc_functional']
                            },
                            "size": 1
                        }
                    }
                }
            }
        }
    }

    if after is not None:
        request['aggs']['results']['composite']['after'] = after
    res = infrastructure.elastic_client.search(index=config.elastic.index_name, body=request)

    if len(res['aggregations']['results']['buckets']) == 0:
        break

    after = res['aggregations']['results']['after_key']
    for material_bucket in res['aggregations']['results']['buckets']:
        material_id = material_bucket['key']['materials']
        entry = material_bucket['calcs']['hits']['hits'][0]['_source']
        upload_id = entry['upload_id']
        calc_id = entry['calc_id']
        results.append(dict(material_id=material_id, upload_id=upload_id, calc_id=calc_id))

    print('.')

results.sort(key=lambda item: item['upload_id'])
print(json.dumps(results, indent=2))
