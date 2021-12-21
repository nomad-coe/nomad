# Using the Python library

NOMAD provides a Python package called `nomad-lab`.
## Install

The package is hosted on [pypi](https://pypi.org/project/nomad-lab/)
and you can install it with *pip* (or conda).

To install the latest stable pypi release, simply use pip:
```sh
pip install nomad-lab
```

!!! attention
    Since NOMAD v1 is still in beta, there is no official pypi release. `pip install` will
    still give you the Python library for the NOMAD v0.10.x version.

To install the latest release developer release (e.g. v1) from our servers use:
```sh
pip install nomad-lab --extra-index-url https://gitlab.mpcdf.mpg.de/api/v4/projects/2187/packages/pypi/simple
```

There are different layers of dependencies that you have to install, in order to use certain functions of NOMAD.
The base install above, will only install the necessary packages for
accessing the NOMAD Archive and use the NOMAD metainfo (see access the archive).

Other functions, e.g. using the NOMAD parsers to parse your code output, require additional dependencies.
You can use the [extra] notation to install these extra requirements:

```
pip install nomad-lab[parsing]
pip install nomad-lab[infrastructure]
pip install nomad-lab[dev]
pip install nomad-lab[all]
```
The various extras have the following meaning:

- *parsing*, everything necessary to run the parsers

- *infrastructure*, everything to run NOMAD services

- *dev*, additional tools that are necessary to develop NOMAD

- *all*, all of the above

## Access parsed NOMAD data with `ArchiveQuery`

The `ArchiveQuery` allows you to search for entries and access their parsed *archive* data
at the same time. Furthermore, all data is accessible through a convenient Python interface
based on the [NOMAD metainfo](metainfo.md) rather than plain JSON.

Here is an example:
```py
query = ArchiveQuery(
    query={
        'results.method.simulation.program_name': 'VASP',
        'results.material.elements': ['Ti', 'O'],
        'results.method.simulation.geometry_optimization': {
            'convergence_tolerance_energy_difference:lt': 1e-22
        }
    },
    required={
        'workflow': {
            'calculation_result_ref': {
                'energy': '*',
                'system_ref': {
                    'chemical_composition_reduced': '*'
                }
            }
        }
    },
    parallel=10,
    max=100)
```

This instantiates an `ArchiveQuery`. You can print some details about the query:

```py
print(query)
```

This gives you a general overview about the query. For example what search is used on
the NOMAD API. How many entries were found. What was already downloaded, etc.
```
Query: {
  "and": [
    {
      "results.method.simulation.program_name": "VASP",
      "results.material.elements": [
        "Ti",
        "O"
      ],
      "results.method.simulation.geometry_optimization": {
        "convergence_tolerance_energy_difference:lt": 1e-22
      }
    },
    {
      "quantities": [
        "run.system.chemical_composition_reduced",
        "run.calculation.system_ref",
        "run.calculation.energy",
        "workflow",
        "workflow.calculation_result_ref"
      ]
    }
  ]
}
Total number of entries that fulfil the query: 252
Number queried entries: 252
Number of entries loaded in the last api call: 70
Bytes loaded in the last api call: 53388
Bytes loaded from this query: 53388
Number of downloaded entries: 70
Number of made api calls: 1
```

This `ArchiveQuery` is not downloaded all archive data immediately. More and more data
will be downloaded as you iterate through the query:
```py
for result in query:
    calc = result.workflow[0].calculation_result_ref
    formula = calc.system_ref.chemical_composition_reduced
    total_energy = calc.energy.total.value.to(units.eV)
    print(f'{formula}: {total_energy}')
```

The resulting output can look like this:
```
O10K2Ti3La2: -136.76387842 electron_volt
Li2O10Ti3La2: -139.15455203 electron_volt
O8Ti4: -107.30373862 electron_volt
O8Ca2Ti4: -116.52240913000001 electron_volt
...
```

Let's discuss the used `ArchiveQuery` parameters:

- `query`, this is an arbitrary API query as discussed in the under [Queries in the API section](api.md#queries).
- `required`, this optional parameter allows you to specify what parts of an archive you need. This is also
described in under [Access archives in API section](api.md#access-archives).
- `per_page`, with this optional parameter you can determine, how many results are downloaded at a time. For bulk downloading many results, we recommend ~100. If you are just interested in the first results a lower number might increase performance.
- `max`, with this optional parameter, we limit the maximum amount of entries that are downloaded, just to avoid accidentally iterating through a result set of unknown and potentially large size.
- `owner` and `auth`, allows you to access private data or specify you only want to
query your data. See also [owner](api.md#owner) and [auth](api.md#authentication) in the API section. Her is an example with authentication:
```py
from nomad.client import ArchiveQuery, Auth

query = ArchiveQuery(
    owner='user',
    required={
        'run': {
            'system[-1]': '*'
        }
    },
    authentication=Auth(user='yourusername', password='yourpassword'))
```

The archive query object can be treated as a Python list-like. You use indices and ranges to select results. Each result is a Python object. The attributes of these objects are
determined by NOMAD's schema, [the metainfo and it's Python interface](metainfo).
This energy value is a number with an attached unit (Joule), which can be converted to something else (e.g. eV). {{ metainfo_data() }}

The create query object keeps all results in memory. Keep this in mind, when you are accessing a large amount of query results.

## Use NOMAD parser locally

If you install `nomad-lab[parsers]`, you can use the NOMAD parsers locally on your computer.
To use the NOMAD parsers from the command line, you can use the parse CLI command. The parse command will automatically match the right parser to your code output file and run the parser. There are two output formats, `--show-metadata` (a JSON representation of the basic metadata), `--show-archive` (a JSON representation of the full parse results).

```sh
nomad parser --show-archive <path-to-your-mainfile-code-output-file>
```

You can also use the NOMAD parsers from within Python. This will give you the parse results as metainfo objects to conveniently analyse the results in Python. See metainfo for more details on how to use the metainfo in Python.

```python
import sys
from nomad.client import parse, normalize_all

# match and run the parser
archive = parse(sys.argv[1])
# run all normalizers
normalize_all(archive)

# get the 'main section' section_run as a metainfo object
section_run = archive.run[0]

# get the same data as JSON serializable Python dict
python_dict = section_run.m_to_dict()
```


You can also clone a parser project and use this to debug or fix a parser:
```sh
git clone https://github.com/nomad-coe/nomad-parser-vasp.git
cd nomad-parser-vasp
git checkout metainfo-refactor
python -m nomad.cli nomad parser --show-archive <path-to-your-vasp-file>
```

Our parsers are hosted in github. They are in the [nomad-coe](https://github.com/nomad-coe) organization. They are typically named `nomad-parser-<code-name>`. The parser version
that fits the NOMAD v1 metainfo schema is typically in the `metainfo-refactor` branch.
Run the CLI with `python -m nomad.cli` to automatically include the current working directory
in the Python path. This will use the cloned parser code over the installed parser code.