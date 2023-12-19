# How to write a parser

NOMAD uses parsers to convert raw code input and output files into NOMAD's common archive
format. This is the documentation on how to develop such a parser.

## Getting started

Let's assume we need to write a new parser from scratch.

First we need to install the `nomad-lab` Python package to get the necessary libraries.
See [how to get started](../develop/setup.md) for the installation guide:

```shell
pip install nomad-lab
```

We prepared an example parser in a github repository to learn how to write parsers. Clone it by:

```shell
git clone https://github.com/nomad-coe/nomad-parser-example.git --branch hello-world
```

Alternatively, you can fork the example project on GitHub to create your own parser. Clone
your fork accordingly.

The project structure should be:

```txt
├── example
│   ├── exampleparser
│   │   ├── __init__.py
│   │   ├── __main__.py
│   │   ├── metainfo.py
│   │   ├── parser.py
│   ├── LICENSE.txt
│   ├── README.md
│   ├── setup.py
```

Next, you should install your new parser with pip. The `-e` parameter installs the parser
in *development*. This means you can change the sources without having to reinstall:

```shell
cd example
pip install -e .
```

The main parser class is found in `exampleparser/parser.py`:

```python
class ExampleParser:
    def parse(self, mainfile: str, archive: EntryArchive, logger):
        # Log a hello world, just to get us started. TODO remove from an actual parser.
        logger.info('Hello World')

        run = archive.m_create(Run)
        run.program = Program(name='EXAMPLE')
```

A parser is a simple Python module containing a single class. For simulation, the convention
for the class name is `<CodeName>Parser` e.g. `VASPParser`. It has a main function, `parse`
which takes the path to the mainfile and an empty `EntryArchive` object as input to be
populated with the parsed quantities. The development of parsers is up to each user, and
will heavily depend on what the user wants to parse. In the simple example above, we created
a logger info entry and populated the archive with a root section called `Run`. We then
created the program section and set the program name to `Example`.

You can run the parser (see the included `__main__.py`) with the path to the file to be
parsed as argument:

```shell
python -m exampleparser tests/data/example.out
```

The output show the log entry and the minimal archive with a `run` section and the
respective `program.name` as in the following:

```json
{
  "run": [
    {
      "program": {
        "name": "EXAMPLE"
      }
    }
  ]
}
```

## Match your raw file

!!! warning "Attention"
    This part of the documentation should be more substantiated. There will be a section
    about this topic soon.

## Parsing test files

We will now show you how to parse ASCII files containing some structure information, a
typical output of simulation codes.

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
and its version (`v2`). Then there are two systems sections (`system 1` and `system 2`)
separated with a string containing a code-specific value `magic source`. Each system section
contains data about the atom positions (`sites`), the lattice information (`latice`),
and a variable `energy`.

In order to convert the information from this file into the archive, we first have to
parse the necessary quantities: the date, system, energy, etc. The `nomad-lab` Python
package provides a `text_parser` module for declarative parsing of text files. You can
define text file parsers as in the following:

```python
def str_to_sites(string):
    sym, pos = string.split('(')
    pos = np.array(pos.split(')')[0].split(',')[:3], dtype=float)
    return sym, pos


calculation_parser = TextParser(quantities=[
    Quantity('sites', r'([A-Z]\([\d\.\, \-]+\))', str_operation=str_to_sites),
    Quantity(
        Atoms.lattice_vectors,
        r'(?:latice|cell): \((\d)\, (\d), (\d)\)\,?\s*\((\d)\, (\d), (\d)\)\,?\s*\((\d)\, (\d), (\d)\)\,?\s*',
        repeats=False),
    Quantity('energy', r'energy: (\d\.\d+)'),
    Quantity('magic_source', r'done with magic source\s*\*{3}\s*\*{3}\s*[^\d]*(\d+)', repeats=False)])

mainfile_parser = TextParser(quantities=[
    Quantity('date', r'(\d\d\d\d\/\d\d\/\d\d)', repeats=False),
    Quantity('program_version', r'super\_code\s*v(\d+)\s*', repeats=False),
    Quantity(
        'calculation', r'\s*system \d+([\s\S]+?energy: [\d\.]+)([\s\S]+\*\*\*)*',
        sub_parser=calculation_parser,
        repeats=True)
])
```

The quantities to be parsed can be specified as a list of `Quantity` objects in `TextParser`.
Each quantity should have a name and a *regular expression (re)* pattern to match the value.
The matched value should be enclosed in a group(s) denoted by `(...)`. In addition, we can
specify the following arguments:

1. `findall (default=True)` Switches simultaneous matching of all quantities using `re.findall`.
In this case, overlap between matches is not tolerated, i.e. two quantities cannot share the same
block in the file. If this cannot be avoided, set `findall=False` switching to`re.finditer`.
This will perform matching one quantity at a time which is slower but with the benefit that
matching is done independently of other quantities.
2. `repeats (default=False)` Switches finding multiple matches for a quantity. By default,
only the first match is returned.
3. `str_operation (default=None)` An external function to be applied on the matched value
to perform more specific string operations.
4. `sub_parser (default=None)` A nested parser to be applied on the matched block. This can
also be a `TextParser` object with a list of quantities to be parsed or [other `FileParser` objects](#other-fileparser-classes).
5. `dtype (default=None)` The data type of the parsed value.
6. `shape (default=None)` The shape of the parsed data.
7. `unit (default=None)` The pint unit of the parsed data.
8. `flatten (default=True)` Switches splitting the parsed string into a flat list.
9. `convert (default=True)` Switches automatic conversion of parsed value.
10. `comment (default=None)` String preceding a line to ignore.

A `metainfo.Quantity` object can also be passed as first argument in place of name in order
to define the data type, shape and unit for the quantity. `TextParser returns a dictionary
of key-value pairs, where the key is defined by the name of the quantities and the value is
based on the matched re pattern.

To parse a file, simply do:

```shell
mainfile_parser.mainfile = mainfile
mainfile_parser.parse()
```

This will populate the `mainfile_parser` object with parsed data and it can be accessed
like a Python dict with quantity names as keys or directly as attributes:

```python
run = Run()
run.program = Program(
  name='super_code', version=mainfile_parser.get('program_version'))
date = datetime.datetime.strptime(mainfile_parser.date, '%Y/%m/%d')
run.program_compilation_datetime = date.timestamp()

for calculation in mainfile_parser.get('calculation', []):
    system = System(atoms=Atoms())

    system.atoms.lattice_vectors = calculation.get('lattice_vectors')
    sites = calculation.get('sites')
    system.atoms.labels = [site[0] for site in sites]
    system.atoms.positions = [site[1] for site in sites]
    run.system.append(system)

    calc = Calculation(energy=Energy())
    calc.system_ref = system
    calc.energy.total = EnergyEntry(value=calculation.get('energy') * units.eV)
    magic_source = calculation.get('magic_source')
    if magic_source is not None:
        calc.x_example_magic_value = magic_source
    run.calculation.append(calc)
archive.run.append(run)
```

When the parser is run on the given example file:

```shell
python -m exampleparser tests/data/example.out
```

you should get a more comprehensive archive with all the provided information..

## Extending the Metainfo

The NOMAD Metainfo defines the schema of each archive. There are predefined schemas for both
simulation `nomad.datamodel.metainfo.simulation` and experimental data
`nomad.datamodel.metainfo.eln`. The sections `Run`, `System`, and `Calculation` in the
example are taken from the simulation metainfo definitions. While this covers most
of the data usually provided in code input/output files, some data are typically
community-specific and applies only to a certain type of codes or methodologies.For these
cases, we allow for the extension of the definitions like this (`exampleparser/metainfo.py`):

```python
# We extend the existing common definition of a section "Calculation"
class ExampleCalculation(Calculation):
    # We alter the default base class behavior to add all definitions to the existing
    # base class instead of inheriting from the base class
    m_def = Section(extends_base_section=True)

    # We define an additional example quantity. Use the prefix x_<parsername>_ to denote
    # non common quantities.
    x_example_magic_value = Quantity(type=int, description='The magic value from a magic source.')
```

## Testing a parser

We developed an initial parser on some example data, and learned how to print out the
output in an archive format. As a good software development practice, we have to add
testing the parser for future maintenance and to ease the future development.

We use the Python unit test framework `pytest`.
A typical test would take one example file, parse it, and check assertions about the
output:

```python
def test_example():
    parser = ExampleParser()
    archive = EntryArchive()
    parser.parse('tests/data/example.out', archive, logging)

    run = archive.run[0]
    assert len(run.system) == 2
    assert len(run.calculation) == 2
    assert run.calculation[0].x_example_magic_value == 42
```

You can run all tests in the `tests` directory like this:

```shell
python -m pytest -svx tests
```

You should define individual test cases with example files that demonstrate certain
features of the underlying code/format.

## Other FileParser classes
Aside from `TextParser`, other `FileParser` classes are also defined. These include:

1. `DataTextParser` uses the `numpy.loadtxt` function to load a structured data file.
The loaded data can be accessed from property `data`.

2. `XMLParser` uses the ElementTree module to parse an XML file. The `parse` method of
the parser takes in an XPath-style key to access individual quantities. By default,
automatic data type conversion is performed, which can be switched off by setting
`convert=False`.

## Add the parser to NOMAD
NOMAD has to manage multiple parsers and must decide during processing which parsers to run
on which files. To accomplish this, specific parser attributes are matched to a
file. These are specified by interfacing the parser with `MatchingParser`. There are a
couple of ways to do this, first as a plug-in (`nomad.config.__init__.py::plugins`) and
second, directly adding it to the list of parsers (`nomad.parsing.parsers.py::parsers`),
the former being the preferred route. See [how to write parser plug-ins](plugins_dev.md#develop-a-parser-plugin)
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

Not all of these attributes have to be used. Those that are given must all match in order
to use the parser on a file.

The NOMAD infrastructure keeps a [list of parser](../../reference/parsers.md#supported-parsers) objects (in
`nomad/parsing/parsers.py::parsers`). These parsers are considered in the order they
appear in the list. The first matching parser is used to parse a given file.

While each parser project should provide its own tests, a test should be
added to the infrastructure parser tests (`tests/parsing/test_parsing.py`) to guarantee that
the processing runs through.

Once the parser is successfully installed and added, it becomes also available through the
command line interface and normalizers are applied as well:

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
