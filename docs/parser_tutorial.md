# How to write a parser for nomad@FAIRDI

## The parser project

First copy an existing parser project as a template. E.g. the vasp parser. Change
the parser metadata in ``setup.py``. You can already install it with ``-e``: from
the main dir of the new parser:

```
pip install -e .
```

The project structure should be
```
myparser/myparser/__init__.py
myparser/test/example_file.out
myparser/LICENSE.txt
myparser/README.txt
myparser/setup.py
```

## Basic skeleton

The following is an example for simple almost empty skeleton for a parser. The
nomad@FAIRDI infrastructure will use the ``ParserInterface`` implementation to use
the parser. The same class can also be used in the ``__main__`` part to run the
parse as a standalone application. You can place it into the ``__init__.py`` of
the main module.

```python
import sys

from nomadcore.simple_parser import SimpleMatcher
from nomadcore.baseclasses import ParserInterface, MainHierarchicalParser

from nomad.parsing import LocalBackend

class VaspOutcarParser(ParserInterface):

    def get_metainfo_filename(self):
        """ The parser specific metainfo. This file must be part of the nomad-meta-info. """
        return 'vasp.nomadmetainfo.json'

    def get_parser_info(self):
        """ Basic info about parser used in archive data and logs. """
        return {
            'name': 'vaspoutcar_parser',
            'version': '1.0.0'
        }

    def setup_version(self):
        """ Can be used to call :func:`setup_main_parser` differently for different code versions. """
        self.setup_main_parser(None)

    def setup_main_parser(self, _):
        """ Setup the actual parser (behind this interface) """
        self.main_parser = MainParser(self.parser_context)


class MainParser(MainHierarchicalParser):
    """
    A MainHierarchicalParser uses a hierarchy of regular expressions to parse
    structured text files.
    """
    def __init__(self, parser_context, *args, **kwargs):
        super().__init__(parser_context, *args, **kwargs)

        # The root_matcher is used to parse the file with a hierarchy of SimpleMatcher
        self.root_matcher = SimpleMatcher(
            name='root',  # a name, just for debugging
            startReStr=r'vasp.',  # the regexp that triggers this matcher
            weak=True,  # True will not consume the matched input, False will
            sections=['section_run'],  # Sections to open and close, 'section_run' is the main section
            subMatchers=[  # Matcher that should be applied after startReStr has matched
                # The group (?P<program_version> ...) will store the matched string
                # in backend for the meta-info quantity 'program_version', which is
                # pars of the 'section_run'
                SimpleMatcher(r'\svasp.(?P<program_version>\d+.\d+.\d+)\s.*')
            ]
        )

if __name__ == "__main__":
    # instantiate the parser via its interface with a LocalBackend
    parser = VaspOutcarParser(backend=LocalBackend)
    # call the actual parsing with the given mainfile
    parser.parse(sys.argv[1])
    # print the results stored in the LocalBackend
    print(parser.parser_context.super_backend)
```

Try to run the file against a test mainfile in `test`:
```
python myparser/__init__.py test/example_file.out
```

or if you installed your new parser in your environment

```
python -m myparser test/example_file.out
```

## Simple matcher
The central thing in a parser using this infrastructure is the SimpleMatcher class.
It is used to defines objects that matches some lines of a file and extracts some values out of them.

The simplest kind of matcher looks like this

```python
    SimpleMatcher(
        r"\s*\|\s*Kinetic energy\s*:\s*(?P<electronic_kinetic_energy_scf__hartree>[-+0-9.eEdD]+) *Ha\s*(?P<fhi_aims_electronic_kinetic_energy_scf_eV>[-+0-9.eEdD]+) *eV")
```


This matcher uses a single regular expression ([regular expressions documentation](https://docs.python.org/2/library/re.html)) to match a line. An online tool to quickly verify regular expressions and to see what they match can be found [here](https://regex101.com/#python).

Note the following things:

* we use a raw string (starting with r"), in it \ does not need to be escaped, otherwise we would need to double every backslash
* we use \s as general space (it matches tab, space,...)
* we extract two expressions using named groups `(?P<nameOfTheGroup>expressionToExtract)`
* group names can have the units of the value they match given after two underscores (`__hartree` means that the values is in hartree).
* it is worthwhile to import SimpleMatcher as SM, to have a more concise code, from now on this is assumed to be the case

The two expressions extracted are automatically assigned to the corresponding meta infos, namely ``electronic_kinetic_energy_scf`` and ``fhi_aims_electronic_kinetic_energy_scf_eV``.
If the value is a scalar then looking at the definition of the meta information the correct type is extracted from the string and the value passed on to the backend.
``electronic_kinetic_energy_scf`` like all code independent values uses SI units, as the matched value is declared to be in hartree (with ``__hartree``) it is automatically converted before storing it.

A matcher can also begin a group of related expressions defined again through other matchers, for example in

```python
SM(name = 'ProgramHeader',
    startReStr = r"\s*Invoking FHI-aims \.\.\.\s*",
    subMatchers = [
        SM(r" *Version *(?P<program_version>[0-9a-zA-Z_.]*)"),
        SM(r" *Compiled on *(?P<fhi_aims_program_compilation_date>[0-9/]+) at (?P<fhi_aims_program_compilation_time>[0-9:]+) *on host *(?P<program_compilation_host>[-a-zA-Z0-9._]+)"),
        SM(name = "nParallelTasks",
        startReStr = r"\s*Using\s*(?P<fhi_aims_number_of_tasks>[0-9]+)\s*parallel tasks\.",
        sections = ["fhi_aims_section_parallel_tasks"],
        subMatchers = [
            SM(name = 'parallelTasksAssignement',
                startReStr = r"\s*Task\s*(?P<fhi_aims_parallel_task_nr>[0-9]+)\s*on host\s*(?P<fhi_aims_parallel_task_host>[-a-zA-Z0-9._]+)\s*reporting\.",
                sections = ["fhi_aims_section_parallel_task_assignement"])
        ]),
        SM(name = 'controlInParser',
        startReStr = r"\s*Parsing control\.in *(?:\.\.\.|\(first pass over file, find array dimensions only\)\.)")
    ])
```

you can see several things:

* you can give a name to the matcher for debugging purposes
* you can nest SimpleMatchers by passing a list of them in the subMatchers argument. If the startReStr regular expression matches then parsing continues in the subMatchers.
* subMatchers by default are executed sequentially, and are optional, so after "Invoking FHI-aims ..." the parser might skip ahead directly to "Compiled on ...". From it the parser cannot "go back" and match "Version ..." but it can match *nParallelTasks* as the subMatchers are sequential.
* After having matched the  "Compiled on ..." part the parser might skip ahead to the *nParallelTask* matcher or *controlInParser*, but not to *parallelTasksAssignement*: the startReStr of a matcher has to match before its subMatchers are considered.
* the *nParallelTasks* matcher opens the section ``fhi_aims_section_parallel_tasks``.
  A section can be seen like a context or a dictionary that starts before the match and is closed when the matcher and all its subMatchers are exited.
  Thus information extracted by *nParallelTasks* (``fhi_aims_number_of_tasks``) or its subMatchers (for example ``fhi_aims_parallel_task_nr``) whose meta information says that the are in the ``fhi_aims_section_parallel_tasks`` (both in this case) get added to this newly created section (again think dictionary).

To recapitulate a parser knows its current position, i.e. which (sub) matcher just matched, and from it knows which matchers can come next.
By default this means going inside the current matcher to one of its direct subMatchers, then going forward at the same level and finally going forward in the superMatchers until the root one is reached.

This matching strategy copes well with print level dependent output (as matchers are optional if absent the parser automatically skips ahead), and group of lines that are not explicitly closed (matching something after the group implicitly closes it) but does not cover all use cases.
For this reason a SimpleMatcher has various arguments to tweak the matching strategy.

* *subFlags* can be set to SubFlags.Unordered instead of the default Sequenced. This makes the parser consider the direct submatchers in any order. This is for example helpful for input files
* if *repeats* is true this matcher is expected to repeat (for example a matcher starting an scf iteration)
* if *endReStr* is provided, the matcher will finish when this regular expression is encountered. Notice that this will not  finish a SimpleMatcher with *repeats=True*, but will only tell when a new repetition should start. If you want to completely stop a repeating SimpleMatcher when a certain regexp is encountered, you should put it inside another SimpleMatcher with a certain *endReStr* and *forward=True*.
* matching a matcher with *weak* = True is attempted unless it is the expected next possible match (does not "steal" the position of the parser). Useful to match patterns that are not distinctive, and could steal the place incorrectly. E.g. the weak flag can be given to a parent matcher if it shares a regular expression with some of it's submatchers. This way the submatchers will have elevated priority over the parent.
* if *required* is true that matcher is expected o match if the enclosing matcher does. This allows to detect errors in parsing early. Currently not enforced.
* if *floating* is true then this section does not steal the position of the parser, it is matched, but then the position of the parser reverts to its old place. This is useful for example for low level debugging/error messages that can happen at any time, but are not supposed to "disturb" the flow of the other output. They are valid from the point they are read until the exit from the enclosing matcher
* *forwardMatch* does not eat the input and forwards it to the adHoc and the subMatchers. This is useful for adHoc parsers, or to have a group that can start with x|y, and then have submatchers for x and y in it. You have to take care that something matches or eats the value if you use *forwardMatch*, otherwise if used in conjunction with repeat it might loop forever.
* *adHoc* lets you define a function that takes over parsing. This function is called with a single argument, the parser. parser.fIn is the PushbackLineFile object from which you can readline() or pushbackLine(lineNotMatched). parser.lastMatch is a dictionary with all groups from startReStr after type and unit conversion, which can be used in simple adHoc functions (e.g. for definition of custom units, see below.
parser.backend is the backend where you can emit the things you parse. With it you have full control on the parsing.
If you return None, then the parser continues normally, but you can also return 2*targetMatcher.index (+ 1 if you are at the end of the matcher) to continue from any other matcher, or -1 to end he parsing.
* *fixedStartValues* can be supplied to set default/fallback values in case optional (by ? quantifier) groups in startReStr are not matched. Is a dictionary
of (metaInfoName: value), unit and type conversion are not applied. Can be also used to emit additional, fixed values not present in startReStr.
* *fixedEndValues* same as fixedStartValues, but applied to endReStr
* The *onClose* argument is a dictionary linking sections into callback functions that are called when a section opened by this SimpleMatcher closes. These callbacks are very similar to the global *onClose* functions that can be defined on the *superContext* object. They have the same call syntax: onClose(backend, gIndex, section), but are specific to this SimpleMatcher and are called before the global counterparts. You can use this e.g. to save the gIndex of *section_single_configuration_calculation* for "frame_sequence_local_frames_ref" when a "section_single_configuration_calculation" is closed within a SimpleMatcher that parses MD or optimization.
* The *onOpen* argument is very similar to *onClose*, just called when section in this SimpleMatcher is opened instead of closed.
* With *startReAction* you can provide a callback function that will get called when the *startReStr* is matched. The callback signature is startReAction(backend, groups). This callback can directly access any matched groups from the regex with argument *groups*. You can use this to e.g. format a line with datetime information into a standard form and push it directly to a corresponding metainfo, or do other more advanced tasks that are not possible e.g. with fixedStartValues.
* Often the parsing is complex enough that you want to create an object to keep track of it and cache various flags, open files,..., to manage the various parsing that is performed. You can store this object as *superContext* of the parser.

## Nomad meta-info

The meta-info has a central role in the parsing.
It represents the conceptual model of the extracted data, every piece of information
extracted by the parsers is described with a meta information.

Nomad Meta Info is described in detail in [here](https://gitlab.rzg.mpg.de/nomad-lab/public-wiki/wikis/nomad-meta-info).

A piece of data is described with a meta information with kindStr = "type_document_content"
(the default, so it can be left out). It must have:

* a *dtypeStr* string that describes the type: floats, integer, boolean, strings, arrays
  of bytes and json dictionaries are supported.
* a *shape* list that gives the dimensions of multidimensional arrays that can contains
  strings to represent symbolic sizes, and is the empty list for scalars.

i.e. at lowest level all data can only be a (multidimensional) array of simple types.

These values can be grouped together through sections (metadata with kindStr = "type_section").
A section can be seen as a dictionary that groups together some values.

When parsing you can open a section (create a new dictionary) and add values to it.
For example ``section_scf_iteration`` groups together all values of an scf iteration.
The next scf iteration will open again the section (create a new dictionary) and add
all the values of the new iteration to it.

Metadata of the same type can be grouped together using meta info
with ``kindStr = "type_abstract_document_content``. This kind of metadata has no
values associated directly with it (and thus is not directly relevant for parsing).
For example all energies inherit from the abstract type energy.

The meta info is described in the [repository nomad-meta-info](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info),
but you can have an [interactive visualization](https://nomad-dev.rz-berlin.mpg.de/ui/index.html)
or [browse](https://nomad-dev.rz-berlin.mpg.de/nmi/info.html) it.

The simple parser infrastructure is closely connected to the meta info.
Normally you will write a code specific meta info in the
[nomad-meta-info repository](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-meta-info)
that imports the common meta info as dependency.

As said all named groups should correspond to a meta info (with kindStr = "type_document_content"),
and the type declared is used to convert the string value extracted by the regular expression.
This works for all scalar values. Sections can be used to organize the values, you can give
the section to open with the sections argument. If an unknown section or group name is used,
the parser will complain and also give a sample meta info that you should add to the meta info file.

The code specific values do not need to be specified in detail, but you should try to connect
them to the common meta info they are related to. For example assigning "settings_XC_functional"
to a code specific setting, allows the analysis part to answer queries like "are all parameters
that influence the XC_functional equal in these two calculations", or "show me the differences
in setting classified by type". This kind of queries is required when building data-sets
that one wants to use for interpolation. Currently there are no methods to answer these queries,
and one simply has to build the datasets himself (often with scripts) to be sure that the
calculations are really consistent.

## What to parse

You should try to parse all input parameters and all lines of the output.

On one side as explained in the previous section this is required to fulfill the Idea of
NOMAD of reusing calculations (and build a tool useful in general for many kind of analysis).
But the other point, that is just as important is to be able of detecting when a parser
fails and needs improvement. This is the only way to keep the parser up to date with
respect to the data that is in the repository.

## Unit conversion

The code independent meta info should use SI units, codes typically do not.
As shown in the examples matchers can convert automatically the value matched by the
their group, if the units of the value read is added to the meta info name with two
underscores. A group name ``g__hartree`` means that the the value read is in hartree and
it should be stored in the meta info named g. Group names cannot contain any special
character, so put more complex units you have to rewrite the expression using just * and
 ^, for example m/s becomes m*s^-1 then you remove the ^ and replace * and - with _,
 so finally m/s becomes ``m_s_1``.

You can find the units that are defined in
``python-common/python/nomadcore/unit_conversions/units.txt``
and some constants in
``python-common/python/nomadcore/unit_conversions/constants.txt``

You might extend these lists if you find important units are missing.
Units names should never use _, use camel case.

Sometime a code uses a unit that is not fixed, but might be defined in the input or some
other part. In this case you can use

```python
    nomadcore.unit_conversion.unit_conversion.register_userdefined_quantity(quantity, units, value=1)
```

from python-common. With it you can (for example in the adHoc callback of a SimpleMatcher)
define the usint usrMyCodeLength (user defined units should always start with "usr" and use
just letters) to be angstrom with

```python
    register_userdefined_quantity("usrMyCodeLength", "angstrom")
```

this call *needs* to be done before any use of that unit. The unit can then be used just
like all others: add __usrMyCodeLength to the group name.

## Backend
The backend is an object can stores parsed data according to its meta-info. The
class :py:class:`nomad.parsing.AbstractParserBackend` provides the basic backend interface.
It allows to open and close sections, add values, arrays, and values to arrays.
In nomad@FAIRDI, we practically only use the :py:class:`nomad.parsing.LocalBackend`. In
NOMAD-coe multiple backend implementations existed to facilitate the communication of
python parsers with the scala infrastructure, including caching and streaming.

## Triggers

When a section is closed a function can be called with the backend, gIndex and the section
object (that might have cached values). This is useful to perform transformations on
the data parsed before emitting it.

The simplest way to achieve this is to define methods called onClose and then the section name in the object that you pass as superContext.
For example

```python
    def onClose_section_scf_iteration(self, backend, gIndex, section):
        logging.getLogger("nomadcore.parsing").info("YYYY bla gIndex %d %s", gIndex, section.simpleValues)
```

defines a trigger called every time an scf iteration section is closed.

## Logging
You can use the standard python logging module in parsers. Be aware that all logging
is stored in a database for analysis. Do not abuse the logging.

## Testing and Debugging
You a writing a python program. You know what to do.

## Added the parser to nomad@FAIRDI

First, you add your parser to the dependencies. Put it into the dependencies folder, then:
```
git submodule add dependencies/parsers/vasp
```

Second, you add your parser to the list of parsers :py:mod:`nomad.parsing`:
```python
parsers = [
    LegacyParser(
        name='parsers/vasp',
        parser_class_name='vaspparser.VaspOutcarParser',
        main_file_re=r'^OUTCAR(\.[^\.]*)?$',
        main_contents_re=(r'^\svasp\..*$')
    ),
]
```

The two regular expressions `main_file_re` and `main_contents_re` are used to match
parsers with potential mainfiles. For each downloaded file, we check if the
`main_file_re` matches the file name and if the `main_contents_re` matches the beginning
of the file's contents. If it matches, the parser will be used to parse it.

Third, add a simple example output file to the nomad@FAIRDI tests in `tests.test_parsing`:
```python
parser_examples = [
    ('parsers/vaspoutcar', 'tests/data/parsers/vasp_outcar/OUTCAR'),
]
```
