# How to write a parser

NOMAD uses parsers to convert raw code input and output files into NOMAD's common archive
format. This is the documentation on how to develop such a parser.

## Getting started

We have prepared an example parser in a github repository to learn how to write parsers. To explore the example, you can clone the project with the command:

```shell
git clone https://github.com/nomad-coe/nomad-parser-example.git --branch hello-world
```

Alternatively, fork the example project on Github and create your own parser.

Once you clone the example project, the file structure is:

```txt
example
   ├── exampleparser
   │   ├── __init__.py
   │   ├── __main__.py
   │   ├── metainfo.py
   │   └── parser.py
   ├── LICENSE.txt
   ├── README.md
   └── pyproject.toml
```

Create a virtual environment (**make sure to use Python 3.9**) and install the new parser:

```shell
cd nomad-parser-example/
python3.9 -m venv .pyenv
source .pyenv/bin/activate
pip install --upgrade pip
pip install -e .
```

The last command will install both the `nomad-lab` Python package and the new parser.
The `-e` parameter installs the parser in *development* mode, which allows you to change the sources without having to reinstall.

The main parser class is found in `exampleparser/parser.py`:

```python
class ExampleParser:
    def parse(self, mainfile: str, archive: EntryArchive, logger):
        # Log a hello world, just to get us started. TODO remove from an actual parser.
        logger.info('Hello World')

        archive.workflow2 = Workflow(name='EXAMPLE')
```

A parser is a simple Python module containing a single class. The convention
for the class name is `<Name>Parser` where `Name` is the file type or code name
e.g. `VASPParser`. It has a main function, `parse` which takes the path to the mainfile
and an empty [`EntryArchive` object](../../reference/glossary.md#archive) as input to be
populated with the parsed quantities. The development of parsers is up to each user, and
will heavily depend on what the user wants to parse. In the simple example above, we created
a logger info entry and populated the archive with a root section called [`Workflow`](../../explanation/data.md#archive-files-a-shared-entry-structure). We then set the workflow name to `EXAMPLE`.

You can run the parser (see the included `__main__.py`) with the path to the file to be
parsed as an argument:

```shell
python -m exampleparser tests/data/example.out
```

The output show the log entry and the minimal archive with a `workflow2` section with the
quantity `name` as in the following:

```json
{
  "workflow2": {
    "name": "EXAMPLE"
  }
}
```

Parsing can also be run within a python script (or Jupyter notebook), e.g., to facilate debugging, with the following code:

```python
from exampleparser import ExampleParser
from nomad.datamodel import EntryArchive
import logging

p = ExampleParser()
a = EntryArchive()
p.parse('tests/data/example.out', a, logger=logging.getLogger())

print(a.m_to_dict())
{'workflow2': {'name': 'EXAMPLE'}}
```
<!-- TODO Add some tips for working with archives in python somewhere -->

## Parsing text files

ASCII text files are amongst the most common files used. Here, we show you how to parse the text by matching specific [regular expressions](https://realpython.com/regex-python/) in these files. For the following example, we will use the project file `tests/data/example.out`:

<!-- TODO can I get rid of this? -->
Check out the `master` branch of the `exampleparser` project,

```shell
git checkout master
```

and examine the file to be parsed in `tests/data/example.out`:

```text
2020/05/15
               *** super_code v2 ***

system 1
--------
sites: H(1.23, 0, 0), H(-1.23, 0, 0), O(0, 0.33, 0)
latice: (0, 0, 0), (1, 0, 0), (1, 1, 0)
energy: 1.29372

*** This was done with magic source                                ***
***                                x°42                            ***


system 2
--------
sites: H(1.23, 0, 0), H(-1.23, 0, 0), O(0, 0.33, 0)
cell: (0, 0, 0), (1, 0, 0), (1, 1, 0)
energy: 1.29372
```

At the top there is some general information such as date, name of the code (`super_code`)
and its version (`v2`). Then is information for two systems (`system 1` and `system 2`),
separated with a string containing a code-specific value `magic source`. Both system sections contain the quantities `sites` and `energy`, but each have a unique quantity as well, `latice` and `cell`, respectively.

In order to convert the information from this file into the NOMAD archive, we first have to
parse the necessary quantities. The `nomad-lab` Python package provides a `text_parser`
module for declarative (i.e., semi-automated) parsing of text files. You can define text file parsers as follows:

```python
def str_to_sites(string):
    sym, pos = string.split('(')
    pos = np.array(pos.split(')')[0].split(',')[:3], dtype=float)
    return sym, pos


calculation_parser = TextParser(
    quantities=[
        Quantity(
            'sites',
            r'([A-Z]\([\d\.\, \-]+\))',
            str_operation=str_to_sites,
            repeats=True,
        ),
        Quantity(
            Model.lattice,
            r'(?:latice|cell): \((\d)\, (\d), (\d)\)\,?\s*\((\d)\, (\d), (\d)\)\,?\s*\((\d)\, (\d), (\d)\)\,?\s*',
            repeats=False,
        ),
        Quantity('energy', r'energy: (\d\.\d+)'),
        Quantity(
            'magic_source',
            r'done with magic source\s*\*{3}\s*\*{3}\s*[^\d]*(\d+)',
            repeats=False,
        ),
    ]
)

mainfile_parser = TextParser(
    quantities=[
        Quantity('date', r'(\d\d\d\d\/\d\d\/\d\d)', repeats=False),
        Quantity('program_version', r'super\_code\s*v(\d+)\s*', repeats=False),
        Quantity(
            'calculation',
            r'\s*system \d+([\s\S]+?energy: [\d\.]+)([\s\S]+\*\*\*)*',
            sub_parser=calculation_parser,
            repeats=True,
        ),
    ]
)
```

The quantities to be parsed can be specified as a list of `Quantity` objects in `TextParser`.
Each quantity should have a name and a *regular expression (re)* pattern to match the value.
The matched value should be enclosed in a group(s) denoted by `(...)`. In addition, we can
specify the following arguments:

- `findall (default=True)`: Switches simultaneous matching of all quantities using `re.findall`.
In this case, overlap between matches is not tolerated, i.e. two quantities cannot share the same
block in the file. If this cannot be avoided, set `findall=False` switching to`re.finditer`.
This will perform matching one quantity at a time which is slower but with the benefit that
matching is done independently of other quantities.
- `repeats (default=False)`: Switches finding multiple matches for a quantity. By default,
only the first match is returned.
- `str_operation (default=None)`: An external function to be applied on the matched value
to perform more specific string operations. In the above example, we defined `str_to_sites` to
convert the parsed value of the atomic sites.
- `sub_parser (default=None)`: A nested parser to be applied on the matched block. This can
also be a `TextParser` object with a list of quantities to be parsed or [other `FileParser` objects](#other-fileparser-classes).
- `dtype (default=None)`: The data type of the parsed value.
- `shape (default=None)`: The shape of the parsed data.
- `unit (default=None)`: The pint unit of the parsed data.
- `flatten (default=True)`: Switches splitting the parsed string into a flat list.
- `convert (default=True)`: Switches automatic conversion of parsed value.
- `comment (default=None)`: String preceding a line to ignore.

A `metainfo.Quantity` object can also be passed as first argument in place of name in order
to define the data type, shape, and unit for the quantity. `TextParser` returns a dictionary
of key-value pairs, where the key is defined by the name of the quantities and the value is
based on the matched re pattern.

To parse a file, simply do:
To parse a file, specify the path to such file and call the `parse()` function of `TextParser`:

```python
mainfile_parser.mainfile = mainfile
mainfile_parser.parse()
```

This will populate the `mainfile_parser` object with parsed data and it can be accessed
like a Python dict with quantity names as keys or directly as attributes:

```python
mainfile_parser.get('date')
'2020/05/15'

mainfile_parser.calculation
[TextParser(example.out) --> 4 parsed quantities (sites, lattice_vectors, energy, magic_source), TextParser(example.out) --> 3 parsed quantities (sites, lattice_vectors, energy)]
```

The next step is to write the parsed data into the NOMAD archive. We can use one of the
[predefined schemas plugins](schemas.md#pre-defined-schemas-in-nomad) in NOMAD.
However, to better illustrate the connection between a parser and a schema we will define our own schema in this example (See [How to write a schema in python](./schemas.md#writing-schemas-in-python-compared-to-yaml-schemas) for additional information on this topic). We define a root section called `Simulation` containing two subsections, `Model` and `Output`. The definitions are found in `exampleparser/metainfo/example.py`:

```python
class Model(ArchiveSection):
    m_def = Section()

    n_atoms = Quantity(
        type=np.int32, description="""Number of atoms in the model system."""
    )

    labels = Quantity(
        type=str, shape=['n_atoms'], description="""Labels of the atoms."""
    )

    positions = Quantity(
        type=np.float64, shape=['n_atoms'], description="""Positions of the atoms."""
    )

    lattice = Quantity(
        type=np.float64,
        shape=[3, 3],
        description="""Lattice vectors of the model system.""",
    )

class Output(ArchiveSection):
    m_def = Section()

    model = Quantity(
        type=Reference(Model), description="""Reference to the model system."""
    )

    energy = Quantity(
        type=np.float64,
        unit='eV',
        description="""Value of the total energy of the system.""",
    )


class Simulation(ArchiveSection):
    m_def = Section()

    code_name = Quantity(
        type=str, description="""Name of the code used for the simulation."""
    )

    code_version = Quantity(type=str, description="""Version of the code.""")

    date = Quantity(type=Datetime, description="""Execution date of the simulation.""")

    model = SubSection(sub_section=Model, repeats=True)

    output = SubSection(sub_section=Output, repeats=True)
```
Each of the classes innherit from the base class `ArchiveSection`. This is the abstract class used in NOMAD to define sections and subsections in a schema. The `Model` section is used to store the `sites` and `lattice/cell` information, while the `Output` section is used to store the `energy` quantity.
Each of the classes that we defined is a sub-class of `ArchiveSection`. This is required in order to assign these sections to the `data` section
of the NOMAD archive.

The following is the implementation of the `parse` function of `ExampleParser` to write the parsed quantities from our mainfile parser into the archive:

```python
def parse(self, mainfile: str, archive: EntryArchive, logger):
    simulation = Simulation(
        code_name='super_code', code_version=mainfile_parser.get('program_version')
    )
    date = datetime.datetime.strptime(mainfile_parser.date, '%Y/%m/%d')
    simulation.date = date

    for calculation in mainfile_parser.get('calculation', []):
        model = Model()
        model.lattice = calculation.get('lattice_vectors')
        sites = calculation.get('sites')
        model.labels = [site[0] for site in sites]
        model.positions = [site[1] for site in sites]
        simulation.model.append(model)

        output = Output()
        output.model = model
        output.energy = calculation.get('energy') * units.eV
        magic_source = calculation.get('magic_source')
        if magic_source is not None:
            archive.workflow2 = Workflow(x_example_magic_value=magic_source)
        simulation.output.append(output)
    # put the simulation section into archive data
    archive.data = simulation
```
We first assign the code name and version as well as the date that the simulation was performed. For each of the parsed
calculations, we create a model and an output section to which we write the corresponding
parsed quantities. Finally, we assign the simulation section to the archive data subsection.

Now, [run the parser again](#getting-started) and check that the new archive stores the intended quantities from `tests/data/example.out`.

Additionally, the standard [normalizers](../../explanation/processing.md#normalizing) will be applied as well. This is run automatically during parsing, one can skip these by passing the
argument `skip-normalizers`.

## Extending the Metainfo
There are several built-in schemas NOMAD (`nomad.datamodel.metainfo`).
<!-- ? What about restructuring this part into the idea of "NOMAD has some predefined section and quantities ... please, check HERE... and HERE... for more information and details"? -->
In the example above, we have made use of the base section for workflow and extended
it to include a code-specific quantity `x_example_magic_value`.
```python
# We extend the existing common definition of section Workflow
class ExampleWorkflow(Workflow):
    # We alter the default base class behavior to add all definitions to the existing
    # base class instead of inheriting from the base class
    m_def = Section(extends_base_section=True)

    # We define an additional example quantity. Use the prefix x_<parsername>_ to denote
    # non common quantities.
    x_example_magic_value = Quantity(
        type=int, description='The magic value from a magic source.'
    )
```
<!-- TODO remove x_ notation in the future -->
This is the approach for domain-specific schemas such as for [simulation workflows](https://github.com/nomad-coe/nomad-schema-plugin-simulation-workflow.git). Refer to [how to extend schemas](schemas.md#extending-existing-sections).


## Testing a parser

Good software development practice involves adding tests during parser development in order to catch bugs and extend the maintainaibility of the parser in the future. For this purpose, we use the Python unit test framework `pytest`. A typical test would take one example file,
parse it, and check assertions about the output:

<!-- TODO suggest a more compartmentalized testing framework -->
```python
def test_example():
    parser = ExampleParser()
    archive = EntryArchive()
    parser.parse('tests/data/example.out', archive, logging)

    sim = archive.data
    assert len(sim.model) == 2
    assert len(sim.output) == 2
    assert archive.workflow2.x_example_magic_value == 42
```

Run all the tests in the `tests/` directory with:

```shell
python -m pytest -svx tests
```

You should follow good [python testing best practices](https://realpython.com/python-testing/).
<!-- TODO add some specific guidance here!  -->

## Other FileParser classes
Aside from `TextParser`, other `FileParser` classes are also defined. These include:

- `DataTextParser`: in addition to matching strings as in `TextParser`, this parser uses the `numpy.loadtxt` function to load structured data files. The loaded `numpy.array` data can then be accessed from the property data.

- `XMLParser`: uses the ElementTree module to parse an XML file. The `parse` method of
the parser takes in an XPath-style key to access individual quantities. By default,
automatic data type conversion is performed, which can be switched off by setting
`convert=False`.

## Adding the parser to NOMAD
NOMAD has to manage multiple parsers and must decide during processing which parsers to run
on which files. To accomplish this, specific parser attributes are matched to a
file. These are specified by interfacing the parser with `MatchingParser`. This can be achieved
by either 1. adding it as a plugin (`nomad.config.__init__.py::plugins`) or 2. directly adding it to the list of parsers (`nomad.parsing.parsers.py::parsers`),
the former being the preferred route. See [How to add a plugin to your NOMAD](plugins.md#add-a-plugin-to-your-nomad)
to learn more.

```python
MatchingParserInterface(
    'parsers/example',
    mainfile_contents_re=(r'^\s*#\s*This is example output'),
    mainfile_mime_re=r'(application/.*)|(text/.*)',
    supported_compressions=["gz", "bz2", "xz"],
    mainfile_alternative=False,
    mainfile_contents_dict={'program': {'version': '1', 'name': 'EXAMPLE'}})
```

- `mainfile_mime_re`: A regular expression on the MIME type of files. The parser is run
  only on files with matching MIME type. The MIME type is *guessed* with libmagic.

- `mainfile_contents_re`: A regular expression that is applied to the first 4k characters of a file.
  The parser is run only on files where this matches.

- `mainfile_name_re`: A regular expression that can be used to match against the name and
  path of the file.

- `supported compressions`: A list of [`gz`, `bz2`] if the parser supports compressed
  files.

- `mainfile_alternative`: If `True`, a file is `mainfile` unless another file in the same
  directory matches `mainfile_name_re`.

- `mainfile_contents_dict`: A dictionary to match the contents of the file. If provided,
  it will load the file and match the value of the key(s) provided. One can also specify
  the keys that should be present by using the tags `__has_key`, `__has_all_keys`
  and `__has_only_keys`. For example, one can have
  `{'program': {'__has_all_keys': ['version', 'name']}}` to specify that `version` and `name`
  must be present in the file to be matched.

The NOMAD infrastructure keeps a list of [Supported Parsers](../../reference/parsers.md#supported-parsers) in
`nomad/parsing/parsers.py::parsers`. These parsers are considered in the order they
appear in the list. The first matching parser is used to parse a given file.

Once the parser is successfully installed and added, it will also become available through the NOMAD [Command Line Interface (CLI)](../../reference/cli.md):

```shell
nomad parse tests/data/example.out
```

## Developing an existing parser

A number of parsers are constantly being developed in NOMAD.

| Description                  | Project url                                            |
| ---------------------------- | ------------------------------------------------------ |
| electronic structure codes   | <https://github.com/nomad-coe/electronic-parsers.git>  |
| atomistic codes              | <https://github.com/nomad-coe/atomistic-parsers.git>   |
| workflow engines             | <https://github.com/nomad-coe/workflow-parsers.git>    |
| databases                    | <https://github.com/nomad-coe/database-parsers.git>    |

To refine an existing parser, you should install the parser via the `nomad-lab` package:

```shell
pip install nomad-lab
```

Clone the parser project:

```shell
git clone <parser-project-url>
cd <parser-dir>
```

Either remove the installed parser and `pip install` the cloned version:

```shell
rm -rf <path-to-your-python-env>/lib/python3.9/site-packages/<parser-module-name>
pip install -e .
```

Or set `PYTHONPATH` so that the cloned code takes precedence over the installed code:

```shell
PYTHONPATH=. nomad parse <path-to-example-file>
```

Alternatively, you can also do a full [developer setup](../develop/setup.md) of the NOMAD infrastructure and
enhance the parser there.

## Developing a Parser plugin
In order to make a parser into a [plugin](plugins.md#add-a-plugin-to-your-nomad), one needs to provide the parser metadata. We provide a separate project for the example parser as plugin.
Fork and clone the [parser example project](https://github.com/nomad-coe/nomad-parser-plugin-example){:target="_blank"} as described in [How to install a plugin](plugins.md).

## Parser plugin metadata
{{pydantic_model('nomad.config.models.plugins.Parser', hide=['code_name','code_category','code_homepage','metadata'])}}

