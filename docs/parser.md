# How to write a parser

NOMAD uses parsers to convert raw code input and output files into NOMAD's common
Archive format. This is documentation on how to develop such a parser.

## Getting started

Let's assume we need to write a new parser from scratch.

First we need the install *nomad-lab* Python package to get the necessary libraries:
```sh
pip install nomad-lab
```

We prepared an example parser project that you can work with.
```sh
git clone https://github.com/nomad-coe/nomad-parser-example.git --branch hello-world
```

Alternatively, you can fork the example project on GitHub to create your own parser. Clone
your fork accordingly.

The project structure should be
```none
example/exampleparser/__init__.py
example/exampleparser/__main__.py
example/exampleparser/metainfo.py
example/exampleparser/parser.py
example/LICENSE.txt
example/README.md
example/setup.py
```

Next you should install your new parser with pip. The `-e` parameter installs the parser
in *development*. This means you can change the sources without the need to re-install.
```sh
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

A parser is a simple program with a single class in it. The base class `MatchingParser`
provides the necessary interface to NOMAD. We provide some basic information
about our parser in the constructor. The *main* function `run` simply takes a filepath
and empty archive as input. Now its up to you, to open the given file and populate the
given archive accordingly. In the plain *hello world*, we simple create a log entry and
populate the archive with a *root section* `Run` and set the program name to `EXAMPLE`.

You can run the parser with the included `__main__.py`. It takes a file as argument and
you can run it like this:
```sh
python -m exampleparser tests/data/example.out
```

The output should show the log entry and the minimal archive with one `section_run` and
the respective `program_name`.
```json
INFO     root                 2020-12-02T11:00:52 Hello World
  - nomad.release: devel
  - nomad.service: unknown nomad service
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
structure information in it. As it is typically used by materials science codes.

The on the `master` branch of the example project, we have a more 'realistic' example:
```sh
git checkout master
```

This example imagines a potential code output that looks like this (`tests/data/example.out`):
```none
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

There is some general information at the top and then a list of simulated systems with
sites and lattice describing crystal structures, a computed energy value, an example for
a code specific quantity from a 'magic source'.

In order to convert the information  from this file into the archive, we first have to
parse the necessary quantities: the date, system, energy, etc. The *nomad-lab* Python
package provides a `text_parser` module for declarative text file parsing. You can
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
and a *regular expression* (re) pattern. The matched value should be enclosed in a group(s)
denoted by `(...)`.
By default, the parser uses the findall method of `re`, hence overlap
between matches is not tolerated. If overlap cannot be avoided, one should switch to the
finditer method by passing *findall=False* to the parser. Multiple
matches for the quantity are returned if *repeats=True* (default). The name, data type,
shape and unit for the quantity can also intialized by passing a metainfo Quantity.
An external function *str_operation* can be also be passed to perform more specific
string operations on the matched value. A local parsing on a matched block can be carried
out by nesting a *sub_parser*. This is also an instance of the `UnstructuredTextFileParser`
with a list of quantities to parse. To access a parsed quantity, one can use the *get*
method.

We can apply these parser definitions like this:
```sh
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

You can still run the parse on the given example file:
```sh
python -m exampleparser tests/data/example.out
```

Now you should get a more comprehensive archive with all the provided information from
the `example.out` file.

** TODO more examples an explanations for: unit conversion, logging, types, scalar, vectors,
multi-line matrices **

## Extending the Metainfo
The NOMAD Metainfo defines the schema of each archive. There are pre-defined schemas for
all domains (e.g. `common_dft.py` for electron-structure codes; `common_ems.py` for
experiment data, etc.). The sections `Run`, `System`, an single configuration calculations (`SCC`)
in the example are taken fom `common_dft.py`. While this covers most data that is usually
provide in code input/output files, some data is typically format specific and only applies
to a certain code or method. For these cases, we allow to extend the Metainfo like
this (`exampleparser/metainfo.py`):
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

Until now, we simply run our parse on some example data and manually observed the output.
To improve the parser quality and ease the further development, you should get into the
habit of testing the parser.

We use the Python unit test framework *pytest*:
```sh
pip install pytest
```

A typical test, would take one example file, parse it, and make assertions about
the output.
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
```sh
pytest -svx tests
```

You should define individual test cases with example files that demonstrate certain
features of the underlying code/format.

## Structured data files with numpy
**TODO: examples**

The `DataTextParser` uses the numpy.loadtxt function to load an structured data file.
The loaded data can be accessed from property *data*.

## XML Parser

**TODO: examples**

The `XMLParser` uses the ElementTree module to parse an xml file. The parse method of the
parser takes in an xpath style key to access individual quantities. By default, automatic
data type conversion is performed, which can be switched off by setting *convert=False*.

## Add the parser to NOMAD

NOMAD has to manage multiple parsers and during processing needs to decide what parsers
to run on what files. To decide what parser is use, NOMAD processing relies on specific
parser attributes.

Consider the example, where we use the `MatchingParser` constructor to add additional
attributes that determine for what files the parser is indented:
```python
class ExampleParser(MatchingParser):
    def __init__(self):
        super().__init__(
            name='parsers/example', code_name='EXAMPLE', code_homepage='https://www.example.eu/',
            mainfile_mime_re=r'(application/.*)|(text/.*)',
            mainfile_contents_re=(r'^\s*#\s*This is example output'))
```

- `mainfile_mime_re`: A regular expression on the mime type of files. The parser is only
run on files with matching mime type. The mime-type is *guessed* with libmagic.
- `mainfile_contents_re`: A regular expression that is applied to the first 4k of a file.
The parser is only run on files where this matches.
- `mainfile_name_re`: A regular expression that can be used to match against the name and path of the file.

Not all of these attributes have to be used. Those that are given must all match in order
to use the parser on a file.

The nomad infrastructure keep a list of parser objects (in `nomad/parsing/parsers.py::parsers`).
These parser are considered in the order they appear in the list. The first matching parser
is used to parse a given file.

While each parser project should provide its own tests, a single example file should be
added to the infrastructure parser tests (`tests/parsing/test_parsing.py`).

Once the parser is added, it become also available through the command line interface and
normalizers are applied as well:
```sh
nomad parser tests/data/example.out
```

## Developing an existing parser
To develop an existing parser, you should install all parsers:
```sh
pip install nomad-lab[parsing]
```

Close the parser project on top:
```sh
git clone <parser-project-url>
cd <parser-dir>
```

Either remove the installed parser and pip install the cloned version:
```sh
rm -rf <path-to-your-python-env>/lib/python3.7/site-packages/<parser-module-name>
pip install -e .
```

Or use `PYTHONPATH` so that the cloned code takes precedence over the installed code:
```sh
PYTHONPATH=. nomad parser <path-to-example-file>
```

Alternatively, you can also do a full developer setup of the NOMAD infrastructure and
develop the parser there.
