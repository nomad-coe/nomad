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
The base install above, only installs the necessary packages for
accessing the NOMAD Archive and using the NOMAD metainfo (see access the archive).

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
based on the [NOMAD metainfo](archive.md) rather than plain JSON.

Please check [this page](async_archive.md) for details.

## Use NOMAD parser locally

If you install `nomad-lab[parsers]`, you can use the NOMAD parsers locally on your computer.
To use the NOMAD parsers from the command line, you can use the parse CLI command. The parse command will automatically match the right parser to your code output file and run the parser. There are two output formats, `--show-metadata` (a JSON representation of the basic metadata) and `--show-archive` (a JSON representation of the full parse results).

```sh
nomad parser --show-archive <path-to-your-mainfile-code-output-file>
```

You can also use the NOMAD parsers within Python, as shown below. This will give you the parse results as metainfo objects to conveniently analyze the results in Python. See metainfo for more details on how to use the metainfo in Python.

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


You can also clone a parser project to debug or fix a parser:
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