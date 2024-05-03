# How to write a parser

NOMAD uses parsers to automatically extract information from raw files and output that information into structured [archives](../../reference/glossary.md#archive). Parsers can decide which files act upon based on the filename, mime type or file contents and can also decide into which schema the information should be populated into.

This documentation shows you how to write a plugin entry point for a parser. You should read the [documentation on getting started with plugins](./plugins.md) to have a basic understanding of how plugins and plugin entry points work in the NOMAD ecosystem.

## Getting started

You can use our [template repository](https://github.com/FAIRmat-NFDI/nomad-plugin-template) to create an initial structure for a plugin containing a parser. The relevant part of the repository layout will look something like this:

```txt
nomad-example
   ├── src
   │   ├── nomad_example
   │   │   ├── parsers
   │   │   │   ├── __init__.py
   │   │   │   ├── myparser.py
   ├── LICENSE.txt
   ├── README.md
   └── pyproject.toml
```

See the documentation on [plugin development guidelines](./plugins.md#plugin-development-guidelines) for more details on the best development practices for plugin, including linting, testing and documenting.

## Parser entry point

The entry point defines basic information about your parser and is used to automatically load the parser code into a NOMAD distribution. It is an instance of a `ParserEntryPoint` or its subclass and it contains a `load` method which returns a `nomad.parsing.Parser` instance that will perform the actual parsing. You will learn more about the `Parser` class in the next sections. The entry point should be defined in `*/parsers/__init__.py` like this:

```python
from pydantic import Field
from nomad.config.models.plugins import ParserEntryPoint


class MyParserEntryPoint(ParserEntryPoint):

    def load(self):
        from nomad_example.parsers.myparser import MyParser

        return MyParser(**self.dict())


myparser = MyParserEntryPoint(
    name = 'MyParser',
    description = 'My custom parser.',
    mainfile_name_re = '.*\.myparser',
)
```

Here you can see that a new subclass of `ParserEntryPoint` was defined. In this new class you can override the `load` method to determine how the `Parser` class is instantiated, but you can also extend the `ParserEntryPoint` model to add new configurable parameters for this parser as explained [here](./plugins.md#extending-and-using-the-entry-point).

We also instantiate an object `myparser` from the new subclass. This is the final entry point instance in which you specify the default parameterization and other details about the parser. In the reference you can see all of the available [configuration options for a `ParserEntryPoint`](../../reference/plugins.md#parserentrypoint).

The entry point instance should then be added to the `[project.entry-points.'nomad.plugin']` table in `pyproject.toml` in order for the parser to be automatically detected:

```toml
[project.entry-points.'nomad.plugin']
myparser = "nomad_example.parsers:myparser"
```

## `Parser` class

The resource returned by a parser entry point must be an instance of a `nomad.parsing.Parser` class. In many cases you will, however, want to use the already existing `nomad.parsing.MatchingParser` subclass that takes care of the file matching process for you. This parser definition should be contained in a separate file (e.g. `*/parsers/myparser.py`) and could look like this:

```python
from typing import Dict

from nomad.datamodel import EntryArchive
from nomad.parsing import MatchingParser


class MyParser(MatchingParser):
    def parse(
        self,
        mainfile: str,
        archive: EntryArchive,
        logger=None,
        child_archives: Dict[str, EntryArchive] = None,
    ) -> None:
        logger.info('MyParser called')
```

If you are using the `MatchingParser` interface, the minimal requirement is
that your class has a `parse` function, which will take as input:

 - `mainfile`: Filepath to a raw file that the parser should open and run on
 - `archive`: The [`EntryArchive` object](../../reference/glossary.md#archive) in which the parsing results will be stored
 - `logger`: Logger that you can use to log parsing events into

Note here that if using `MatchingParser`, the process of identifying which files the `parse` method is run against is take care of by passing in the required parameters to the instance in the `load` mehod. In the previous section, the `load` method looked something like this:

```python
    def load(self):
        from nomad_example.parsers.myparser import MyParser

        return MyParser(**self.dict())
```

There we are passing all of the entry configuration options to the parser instance, including things like `mainfile_name_re` and `mainfile_contents_re`. The `MatchingParser` constructor uses these parameters to set up the file matching appropriately. If you wish to take full control of the file matching process, you can use the `nomad.parsing.Parser` class and override the `is_mainfile` function.

## Match your raw file

If you are using the `MatchingParser` interface you can configure which files
are matched directly in the `ParserEntryPoint`. For example to match only certain file extensions and file contents, you can use the `mainfile_name_re` and `mainfile_contents_re` fields:

```python
myparser = MyParserEntryPoint(
    name = 'MyParser',
    description = 'My custom parser.',
    mainfile_name_re = '.*\.myparser',
    mainfile_contents_re = '\s*\n\s*HELLO WORLD',
)
```

You can find all of the available matching criteria in the [`ParserEntryPoint` reference](../../reference/plugins.md#parserentrypoint)

## Running the parser

If you have the plugin package and `nomad-lab` installed in your Python environment, you can run the parser against a file using the NOMAD CLI:

```shell
nomad parse <filepath> --show-archive
```

The output will return the final archive in JSON format.

Parsing can also be run within a python script (or Jupyter notebook), e.g., to facilate debugging, with the following code:

```python
from nomad.datamodel import EntryArchive
from nomad_example.parsers.myparser import MyParser
import logging

p = ExampleParser()
a = EntryArchive()
p.parse('tests/data/example.out', a, logger=logging.getLogger())

print(a.m_to_dict())
```

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
[predefined plugins containing schema packages](schema_packages.md#schema-packages-developed-by-fairmat) in NOMAD.
However, to better illustrate the connection between a parser and a schema we will define our own schema in this example (See [How to write a schema in python](./schema_packages.md#writing-schemas-in-python-compared-to-yaml-schemas) for additional information on this topic). We define a root section called `Simulation` containing two subsections, `Model` and `Output`. The definitions are found in `exampleparser/metainfo/example.py`:

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
Each of the classes inherit from the base class `ArchiveSection`. This is the abstract class used in NOMAD to define sections and subsections in a schema. The `Model` section is used to store the `sites` and `lattice/cell` information, while the `Output` section is used to store the `energy` quantity.
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
This is the approach for domain-specific schemas such as for [simulation workflows](https://github.com/nomad-coe/nomad-schema-plugin-simulation-workflow.git). Refer to [how to extend schemas](schema_packages.md#extending-existing-sections).

## Other FileParser classes
Aside from `TextParser`, other `FileParser` classes are also defined. These include:

- `DataTextParser`: in addition to matching strings as in `TextParser`, this parser uses the `numpy.loadtxt` function to load structured data files. The loaded `numpy.array` data can then be accessed from the property data.

- `XMLParser`: uses the ElementTree module to parse an XML file. The `parse` method of
the parser takes in an XPath-style key to access individual quantities. By default,
automatic data type conversion is performed, which can be switched off by setting
`convert=False`.

## Parsers developed by FAIRmat

The following is a list of plugins containin parsers developed by FAIRmat:

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
