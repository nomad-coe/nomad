import csv
import sys
from inspect import cleandoc as strip


if __name__ == '__main__':
    with open(sys.argv[1], 'rt') as f:
        data = [row for row in csv.reader(f)]

    for item in data:
        print(strip(f'''
- [ ] {item[1]} entries ({item[4]} uploads) x
[**{item[6]}**: *{item[3]}/{item[4]}* (v1)](https://nomad-lab.eu/prod/v1/gui/search/entries/entry/id/{item[3]}/{item[2]}),
[raw files (v0)](https://nomad-lab.eu/prod/rae/gui/entry/id/{item[3]}/{item[2]}/raw)
```
{item[7]}
```'''))
