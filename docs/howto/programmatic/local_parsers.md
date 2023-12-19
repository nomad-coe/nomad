# How to run a parser

You can find a [list of all parsers](../../reference/parsers.md) and supported files in the reference.

First you need to have the `nomad-lab` pypi package installed. You find more detailed
instructions [here](pythonlib.md):

```
pip install nomad-lab
```

## From the command line

You can run NOMAD parsers from the [command line interface](../../reference/cli.md) (CLI).
The parse command will automatically match the right parser to your file and run the parser.
There are two output formats:

- `--show-metadata` a JSON representation of the basic metadata
-  `--show-archive` a JSON representation of the full parse results

```
nomad parse --show-archive <path-to-your-mainfile-code-output-file>
```

!!! note

    If you run into missing dependency errors, you might want to install additional
    dependencies via `pip install nomad-lab[parsing]`. Only a view parsers require
    extra dependencies. Please refer to the parser projects for more details.

To skip the parser matching, i.e. the process that determined which parser fits to
the given file, and state the parser directly, you can use the `--parser` argument
to provide a [parser name](../../reference/parsers.md).

```
nomad parser --parser parsers/vasp <path-to-your-mainfile-code-output-file>
```

To skip the potentially error-prone and depending on your use-case unnecessary
normalization, you can use the `--skip-normalizers` argument:

```
nomad parser --skip-normalizers <path-to-your-mainfile-code-output-file>
```

## From a python program

You can also use the NOMAD parsers within Python, as shown below.
This will give you the parse results as metainfo objects to conveniently analyze the results in Python.
See metainfo for more details on how to use the metainfo in Python.

```python
# requires: nomad-lab
import sys
from nomad.client import parse, normalize_all

# match and run the parser
archives = parse('path/to/you/file')
# run all normalizers
for archive in archives:
    normalize_all(archive)

    # get the 'main section' section_run as a metainfo object
    section_run = archive.run[0]

    # get the same data as JSON serializable Python dict
    python_dict = section_run.m_to_dict()
```
