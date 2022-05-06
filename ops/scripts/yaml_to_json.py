import yaml
import json
from typing import Dict

os_list: Dict[str, str] = {}

with open("examples/data/ExampleUploads.yml") as infile:
    os_list = yaml.load(infile, Loader=yaml.FullLoader)

with open("gui/src/ExampleUploads.json", 'w') as outfile:
    json.dump(os_list, outfile, indent=4)
