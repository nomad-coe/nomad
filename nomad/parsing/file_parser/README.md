# NOMAD file parsing module

The parsing module consists of the `UnstructuredTextFileParser`, `DataTextFileParser`
and `XMLParser` classes to enable the parsing of unstructured text, structured data text,
and xml files, respectively. These classes are based on the FileParser class which
provides the common methods for file handling, and querying the parsed results.

## UnstructuredTextFileParser

The most common type of file that are parsed in NOMAD are unstructured text files which
can be handled using the UnstructuredTextFileParser. The parser uses the `re` module to
match a given pattern for a quantity in the text file. To illustrate the use of this parser,
let us consider a file `super_code.out` with the following contents:

```
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

In order to create a nomad archive from this file, we first have to parse the necessary
quantities which includes the date, system, energy, etc. The following python code
illustrates how can this be achieved. Note that we will be using *parser* to refer to the
file parser and to the code parser that writes the archive.

```python
import datetime
import numpy as np

from nomad.parsing.file_parser import UnstructuredTextFileParser, Quantity
from nomad.datamodel import EntryArchive
from nomad.datamodel.metainfo.public import section_run, section_system, section_single_configuration_calculation

p = UnstructuredTextFileParser()

def str_to_sites(string):
    sym, pos = string.split('(')
    pos = np.array(pos.split(')')[0].split(',')[:3], dtype=float)
    return sym, pos

q_system = Quantity(
    'system', r'\s*system \d+([\s\S]+?energy: [\d\.]+)([\s\S]+\*\*\*)*',
    sub_parser=UnstructuredTextFileParser(quantities=[
        Quantity(
            'sites', r'([A-Z]\([\d\.\, \-]+\))',
            str_operation=str_to_sites),
        Quantity(
            section_system.lattice_vectors,
            r'(?:latice|cell): \((\d)\, (\d), (\d)\)\,?\s*\((\d)\, (\d), (\d)\)\,?\s*\((\d)\, (\d), (\d)\)\,?\s*',
            repeats=False),
        Quantity(
            'energy', r'energy: (\d\.\d+)'),
        Quantity(
            'magic_source', r'done with magic source\s*\*{3}\s*\*{3}\s*([\S]+)',
            repeats=False)]),
    repeats=True)

quantities = [
        Quantity('date', r'(\d\d\d\d\/\d\d\/\d\d)', repeats=False),
        Quantity('program_version', r'super\_code\s*v(\d+)\s*', repeats=False),
        q_system]

p.quantities = quantities
# this returns the energy for system 2
p.system[1].get('energy', unit='hartree')
```

The quantities to be parsed can be specified as a list of `Quantity` objects with a name
and a re pattern. The matched value should be enclosed in a group(s). By default,
the parser uses the findall method of `re`, hence overlap
between matches is not tolerated. If overlap cannot be avoided, one should switch to the
finditer method by passing *findall=False* to the parser. Multiple
matches for the quantity are returned if *repeats=True* (default). The name, data type,
shape and unit for the quantity can also intialized by passing a metainfo Quantity.
An external function *str_operation* can be also be passed to perform more specific
string operations on the matched value. A local parsing on a matched block can be carried
out by nesting a *sub_parser*. This is also an instance of the `UnstructuredTextFileParser`
with a list of quantities to parse. To access a parsed quantity, one can use the *get*
method.

The creation of the archive is implemented in the parse method of the code parser which takes
the mainfile, archive and logger as arguments. The file parser, *out_parser* is
created only in the constructor and subsequent parsing on a different *mainfile* can be
performed by assigning it to the file parser.

```python
class SupercodeParser:
    def __init__(self):
        self.out_parser = UnstructuredTextFileParser()
        self.out_parser.quantities = quantities

    def parse(self, mainfile, archive, logger):
        self.out_parser.mainfile = mainfile
        sec_run = archive.m_create(section_run)
        sec_run.program_name = 'super_code'
        sec_run.program_version = str(self.out_parser.get('program_version'))
        date = datetime.datetime.strptime(
            self.out_parser.get('date'), '%Y/%m/%d') - datetime.datetime(1970, 1, 1)
        sec_run.program_compilation_datetime = date.total_seconds()
        for system in self.out_parser.get('system'):
            sec_system = sec_run.m_create(section_system)
            sec_system.lattice_vectors = system.get('lattice_vectors')
            sites = system.get('sites')
            sec_system.atom_labels = [site[0]  for site in sites]
            sec_system.atom_positions = [site[1] for site in sites]

            sec_scc = sec_run.m_create(section_single_configuration_calculation)
            sec_scc.energy_total = system.get('energy')
            sec_scc.single_configuration_calculation_to_system_ref = sec_system
            magic_source = system.get('magic_source')
            if magic_source is not None:
                sec_scc.message_info_evaluation = magic_source

archive = EntryArchive()

parser = SupercodeParser()
parser.parse('temp.dat', archive, None)

print(archive.m_to_json())
```

## DataTextFileParser
The `DataTextFileParser` uses the numpy.loadtxt function to load an structured data file.
The loaded data can be accessed from property *data*.

## XMLParser
The `XMLParser` uses the ElementTree module to parse an xml file. The parse method of the
parser takes in an xpath style key to access individual quantities. By default, automatic
data type conversion is performed, which can be switched off by setting *convert=False*.


