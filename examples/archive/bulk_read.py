'''
In this example, we select calculations/upload by query and iterate through all
uploads in parallel, reading the uploads in full to extract some information for each
calculation.

The motivation behind this is that the selective access of sections in an archvie might
be slow. If something is read from almost all calculations, it might be faster to
sequentially read all of the upload's archive file.

This is not an API example, but directly accesses archive files.
'''
from typing import Iterable, Tuple, Set, Any
from concurrent.futures import ThreadPoolExecutor
from threading import Event, Lock
from queue import Queue, Empty
import json

from nomad import search, infrastructure, files, config

infrastructure.setup_files()
infrastructure.setup_elastic()


def query() -> Iterable[Tuple[str, Set[str]]]:
    after = None
    while True:
        res = infrastructure.elastic_client.search(index=config.elastic.index_name, body={
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
                        "sources": [
                            {
                                "uploads": {
                                    "terms": {
                                        "field": "upload_id"
                                    },
                                }
                            },
                            {
                                "materials": {
                                    "terms": {
                                        "field": "encyclopedia.material.material_id"
                                    }
                                }
                            }
                        ],
                        "size": 100
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
        })
        print(json.dumps(res, indent=2))
        raise
        # searchRequest = search.SearchRequest()
        # searchRequest.quantity(
        #     name='encyclopedia.material.material_id',
        #     examples=2,
        #     examples_source=['upload_id', 'calc_id', 'dft.code_name', 'dft.xc_functional'],
        #     order_by='upload_id',
        #     after=after)

        # result = searchRequest.execute()['quantities']['encylcopedia.material.material_id']
        # after = result['after']
        # if len(result) == 0:
        #     break

        # for material_id, calcs in result.items():
        # print(json.dumps(result, indent=2))
        # raise
        # calc_ids: Set[str] = set()
        # upload_id = None
        # for entry in searchRequest.execute_scan(order_by='upload_id'):
        #     entry_upload_id, entry_calc_id = entry['upload_id'], entry['calc_id']
        #     if upload_id is not None and upload_id != entry_upload_id:
        #         yield upload_id, calc_ids
        #         upload_id = entry_calc_id
        #         calc_ids = set()

        #     upload_id = entry_upload_id
        #     calc_ids.add(entry_calc_id)

        # if upload_id is not None:
        #     yield upload_id, calc_ids


for _ in query():
    pass

def read_archive(upload_id, calc_ids):
    try:
        upload_files = files.UploadFiles.get(upload_id, lambda *args: True)
        for calc_id in calc_ids:
            with upload_files.read_archive(calc_id) as archive:
                material_id = archive[calc_id]['section_metadata']['encyclopedia']['material']['material_id']
                for run in archive[calc_id].get('section_run', []):
                    for calc in run.get('section_single_configuration_calculation', []):
                        for dos in calc.get('section_dos', []):
                            fingerprint = dos.get('section_dos_fingerprint')
                            if fingerprint:
                                yield {
                                    'upload_id': upload_id,
                                    'calc_id': calc_id,
                                    'material_id': material_id,
                                    'fingerprint': fingerprint}
    except Exception:
        import traceback
        traceback.print_exc()


nworker = 10
nended_worker = 1
upload_queue: Any = Queue(maxsize=100)
result_queue: Any = Queue(maxsize=100)
producer_end = Event()
result_end = Event()
ended_worker_lock = Lock()


def worker():
    global nended_worker
    while not (producer_end.is_set() and upload_queue.empty()):
        try:
            upload_id, calc_ids = upload_queue.get(block=True, timeout=0.1)
        except Empty:
            continue

        for result in read_archive(upload_id, calc_ids):
            result_queue.put(result)

    with ended_worker_lock:
        nended_worker += 1
        if nended_worker == nworker:
            print('result end')
            result_end.set()

    print('end worker')


def writer():
    while not (result_end.is_set() and result_queue.empty()):
        try:
            result = result_queue.get(block=True, timeout=0.1)
        except Empty:
            continue

        print(json.dumps(result, indent=2))

    print('end writer')


def producer():
    for upload_id, calc_ids in query():
        upload_queue.put((upload_id, calc_ids), block=True)

    producer_end.set()
    print('end producer')


with ThreadPoolExecutor(max_workers=nworker + 2) as executor:
    for _ in range(nworker):
        executor.submit(worker)
    executor.submit(producer)
    executor.submit(writer)
