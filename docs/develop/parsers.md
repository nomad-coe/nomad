# How to write a parser

NOMAD uses parsers to convert raw code input and output files into NOMAD's common archive
format. This is the documentation on how to develop such a parser.

## Getting started

Let's assume we need to write a new parser from scratch.

First we need to install the `nomad-lab` Python package to get the necessary libraries:

```shell
pip install nomad-lab
```

We have prepared an example parser project that you can work with:

```shell
git clone https://github.com/nomad-coe/nomad-parser-example.git --branch hello-world
```

Alternatively, you can fork the example project on GitHub to create your own parser. Clone
your fork accordingly.

The project structure should be:

```text
example/exampleparser/__init__.py
example/exampleparser/__main__.py
example/exampleparser/metainfo.py
example/exampleparser/parser.py
example/LICENSE.txt
example/README.md
example/setup.py
```

Next, you should install your new parser with pip. The `-e` parameter installs the parser
in *development*. This means you can change the sources without having to reinstall:

```shell
cd example
pip install -e .
```

The main code file `exampleparser/parser.py` should look like this:

```python
class ExampleParser(MatchingParser):
    def __init__(self):
        super().__init__(name='parsers/example', code_name='EXAMPLE')

    def run(self, mainfile: str, archive: EntryArchive, logger):
        # Log a hello world, just to get us started. TODO remove from an actual parser.
        logger.info('Hello World')

        run = archive.m_create(Run)
        run.program_name = 'EXAMPLE'
```

A parser is a simple program with a single class. The base class `MatchingParser`
provides the necessary interface to NOMAD. We provide some basic information
about our parser in the constructor. The *main* function `run` simply takes a filepath
and an empty archive as input. Now it's up to you, to open the given file and populate the
given archive accordingly. In the plain *hello world*, we simple create a log entry,
populate the archive with a *root section* `Run`, and set the program name to `EXAMPLE`.

You can run the parser with the included `__main__.py`. It takes a file as argument and
you can run it like this:

```shell
python -m exampleparser tests/data/example.out
```

The output should show the log entry and the minimal archive with one `section_run` and
the respective `program_name`:

```json
{
  "section_run": [
    {
      "program_name": "EXAMPLE"
    }
  ]
}
```

## Parsing test files

Let's do some actual parsing. Here we demonstrate how to parse ASCII files with some
structure information in it, as it is typically used by materials science codes.

On the `master` branch of the example project, we have a more 'realistic' example:

```shell
git checkout master
```

This example imagines a potential code output that looks like this
(`tests/data/example.out`):

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

At the top there is some general information. Below that is a list of simulated systems
with "sites" and "lattice" describing the crystal structure, as well as a computed energy
value as an example for a code specific quantity from a 'magic source'.

In order to convert the information from this file into the archive, we first have to
parse the necessary quantities: the date, system, energy, etc. The `nomad-lab` Python
package provides a `text_parser` module for declarative parsing of text files. You can
define text file parsers like this:

```python
def str_to_sites(string):
    sym, pos = string.split('(')
    pos = np.array(pos.split(')')[0].split(',')[:3], dtype=float)
    return sym, pos


calculation_parser = UnstructuredTextFileParser(quantities=[
    Quantity('sites', r'([A-Z]\([\d\.\, \-]+\))', str_operation=str_to_sites),
    Quantity(
        System.lattice_vectors,
        r'(?:latice|cell): \((\d)\, (\d), (\d)\)\,?\s*\((\d)\, (\d), (\d)\)\,?\s*\((\d)\, (\d), (\d)\)\,?\s*',
        repeats=False),
    Quantity('energy', r'energy: (\d\.\d+)'),
    Quantity('magic_source', r'done with magic source\s*\*{3}\s*\*{3}\s*[^\d]*(\d+)', repeats=False)])

mainfile_parser = UnstructuredTextFileParser(quantities=[
    Quantity('date', r'(\d\d\d\d\/\d\d\/\d\d)', repeats=False),
    Quantity('program_version', r'super\_code\s*v(\d+)\s*', repeats=False),
    Quantity(
        'calculation', r'\s*system \d+([\s\S]+?energy: [\d\.]+)([\s\S]+\*\*\*)*',
        sub_parser=calculation_parser,
        repeats=True)
])
```

The quantities to be parsed can be specified as a list of `Quantity` objects with a name
and a *regular expression (re)* pattern. The matched value should be enclosed in a
group(s) denoted by `(...)`.
By default, the parser uses the `findall` method of `re`, hence overlap between matches is
not tolerated. If overlap cannot be avoided, you should switch to the `finditer` method by
passing `findall=False` to the parser. Multiple matches for the quantity are returned if
`repeats=True` (default). The name, data type, shape and unit for the quantity can also be
initialized by passing a `metainfo.Quantity`.
An external function `str_operation` can also be passed to perform more specific string
operations on the matched value. A local parsing on a matched block can be carried out by
nesting a `sub_parser`. This is also an instance of the `UnstructuredTextFileParser` with
a list of quantities to parse. To access a parsed quantity, you can use the `get` method.

We can apply these parser definitions like this:

```shell
mainfile_parser.mainfile = mainfile
mainfile_parser.parse()
```

This will populate the `mainfile_parser` object with parsed data and it can be accessed
like a Python dict with quantity names as keys:

```python
run = archive.m_create(Run)
run.program_name = 'super_code'
run.program_version = str(mainfile_parser.get('program_version'))
date = datetime.datetime.strptime(
    mainfile_parser.get('date'),
    '%Y/%m/%d') - datetime.datetime(1970, 1, 1)
run.program_compilation_datetime = date.total_seconds()

for calculation in mainfile_parser.get('calculation'):
    system = run.m_create(System)

    system.lattice_vectors = calculation.get('lattice_vectors')
    sites = calculation.get('sites')
    system.atom_labels = [site[0] for site in sites]
    system.atom_positions = [site[1] for site in sites]

    scc = run.m_create(SCC)
    scc.single_configuration_calculation_to_system_ref = system
    scc.energy_total = calculation.get('energy') * units.eV
    scc.single_configuration_calculation_to_system_ref = system
    magic_source = calculation.get('magic_source')
    if magic_source is not None:
        scc.x_example_magic_value = magic_source
```

You can still run the parser on the given example file:

```shell
python -m exampleparser tests/data/example.out
```

Now you should get a more comprehensive archive with all the provided information from
the `example.out` file.

## Extending the Metainfo

The NOMAD Metainfo defines the schema of each archive. There are predefined schemas for
all domains (e.g. `common_dft.py` for electronic structure codes; `common_ems.py` for
experimental data, etc.). The sections `Run`, `System`, and the single configuration
calculations (`SCC`) in the example are taken from `common_dft.py`. While this covers most
of the data usually provided in code input/output files, some data is typically
format-specific and applies only to a certain code or method. For these cases, we allow to
extend the Metainfo like this (`exampleparser/metainfo.py`):

```python
# We extend the existing common definition of a section "single configuration calculation"
class ExampleSCC(SCC):
    # We alter the default base class behavior to add all definitions to the existing
    # base class instead of inheriting from the base class
    m_def = Section(extends_base_section=True)

    # We define an additional example quantity. Use the prefix x_<parsername>_ to denote
    # non common quantities.
    x_example_magic_value = Quantity(type=int, description='The magic value from a magic source.')
```

## Testing a parser

Until now, we simply run our parser on some example data and manually observed the output.
To improve the parser quality and ease the further development, you should get into the
habit of testing the parser.

We use the Python unit test framework `pytest`:

```shell
pip install pytest
```

A typical test would take one example file, parse it, and check assertions about the
output:

```python
def test_example():
    parser = rExampleParser()
    archive = EntryArchive()
    parser.run('tests/data/example.out', archive, logging)

    run = archive.section_run[0]
    assert len(run.section_system) == 2
    assert len(run.section_single_configuration_calculation) == 2
    assert run.section_single_configuration_calculation[0].x_example_magic_value == 42
```

You can run all tests in the `tests` directory like this:

```shell
pytest -svx tests
```

You should define individual test cases with example files that demonstrate certain
features of the underlying code/format.

## Structured data files with NumPsy

The `DataTextParser` uses the `numpy.loadtxt` function to load an structured data file.
The loaded data can be accessed from property `data`.

## XML Parser

The `XMLParser` uses the ElementTree module to parse an XML file. The `parse` method of
the parser takes in an XPath-style key to access individual quantities. By default,
automatic data type conversion is performed, which can be switched off by setting
`convert=False`.

## Add the parser to NOMAD

NOMAD has to manage multiple parsers and must decide during processing which parsers
to run on which files. To decide what parser is used, NOMAD processing relies on specific
parser attributes.

Consider the example where we use the `MatchingParser` constructor to add additional
attributes that determine for which files the parser is intended for:

```python
class ExampleParser(MatchingParser):
    def __init__(self):
        super().__init__(
            name='parsers/example', code_name='EXAMPLE', code_homepage='https://www.example.eu/',
            mainfile_mime_re=r'(application/.*)|(text/.*)',
            mainfile_contents_re=(r'^\s*#\s*This is example output'),
            supported_compressions=["gz", "bz2", "xz"],
            mainfile_alternative=False,
            mainfile_contents_dict={'program': {'version': '1', 'name': 'EXAMPLE'}})
```

- `mainfile_mime_re`: A regular expression on the MIME type of files. The parser is run
  only on files with matching MIME type. The MIME type is *guessed* with libmagic.

- `mainfile_contents_re`: A regular expression that is applied to the first 4k of a file.
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

The NOMAD infrastructure keeps a list of parser objects (in
`nomad/parsing/parsers.py::parsers`). These parsers are considered in the order they
appear in the list. The first matching parser is used to parse a given file.

While each parser project should provide its own tests, a single example file should be
added to the infrastructure parser tests (`tests/parsing/test_parsing.py`).

Once the parser is added, it becomes also available through the command line interface and
normalizers are applied as well:

```shell
nomad parse tests/data/example.out
```

## Developing an existing parser

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

Alternatively, you can also do a full developer setup of the NOMAD infrastructure and
enhance the parser there.
