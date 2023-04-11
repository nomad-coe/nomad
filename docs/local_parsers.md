# Run NOMAD parser locally

If you install `nomad-lab[parsing]`, you can use the NOMAD parsers locally on your computer.
To use the NOMAD parsers from the command line, you can use the parse CLI command. The parse command will automatically
match the right parser to your code output file and run the parser. There are two output formats, `--show-metadata` (a
JSON representation of the basic metadata) and `--show-archive` (a JSON representation of the full parse results).

```sh
nomad parse --show-archive <path-to-your-mainfile-code-output-file>
```

You can also use the NOMAD parsers within Python, as shown below.
This will give you the parse results as metainfo objects to conveniently analyze the results in Python.
See metainfo for more details on how to use the metainfo in Python.

```python
# requires: nomad-lab[parsing,infrastructure]
import sys
from nomad.client import parse, normalize_all

# match and run the parser
archives = parse(sys.argv[1])
# run all normalizers
for archive in archives:
    normalize_all(archive)

    # get the 'main section' section_run as a metainfo object
    section_run = archive.run[0]

    # get the same data as JSON serializable Python dict
    python_dict = section_run.m_to_dict()
```

You can also clone a parser project to debug or fix a parser:

```sh
git clone https://github.com/nomad-coe/nomad-parser-vasp.git
cd nomad-parser-vasp
git checkout metainfo-refactor
python -m nomad.cli nomad parse --show-archive <path-to-your-vasp-file>
```

Our parsers are hosted in GitHub.
They are in the [nomad-coe](https://github.com/nomad-coe) organization.
They are typically named `nomad-parser-<code-name>`.
The parser version that fits the NOMAD v1 metainfo schema is typically in the `metainfo-refactor` branch.
Run the CLI with `python -m nomad.cli` to automatically include the current working directory in the Python path.
This will use the cloned parser code over the installed parser code.