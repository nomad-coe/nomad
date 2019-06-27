import pymongo
import zipfile
import sys

from nomad import parsing

client = pymongo.MongoClient()
packages = client['coe_migration']['package']

def check(upload_id, mainfile):
    content = None
    for package in packages.find(dict(upload_id=upload_id)):
        package_path = package['package_path']
        with zipfile.ZipFile(package_path, 'r') as zf:
            try:
                with zf.open(mainfile, 'r') as f:
                    content = f.read(5000)

            except KeyError:
                pass

    if content is None:
        print('mainfile does not exist')
        sys.exit(1)

    match = None
    for parser in parsing.parsers:
        if parser.is_mainfile(mainfile, 'text/plain', content, None):
            match = parser

    if match is None:
        try:
            print(content.decode('utf-8'))
        except Exception:
            print('not unicode decodable, probably binary file')

    return match is not None

import json
with open('local/missing_calcs_data.json') as f:
    data = json.load(f)

for cause in data['others'] + data['no_calcs']:
    if 'investigated_cause' not in cause and 'phonopy' not in cause['example_mainfile'] and not check(cause['source_upload_id'], cause['example_mainfile']):
        input(cause)
