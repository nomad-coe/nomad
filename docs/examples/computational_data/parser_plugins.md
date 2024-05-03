# Guide to computational parser plugins

NOMAD uses parsers to convert raw data (for example, output from computational software, instruments,
or electronic lab notebooks) into NOMAD's common Archive format. This page provides a guide to the existing standard computational parsers in NOMAD.

## Parser organization

!!! note
    The majority of NOMAD's computational parsers do not yet exist as plugins. Instead, they are linked to the nomad-lab software as submodules. We will be migrating these parsers to proper plugins in the near future, and may reorganize the projects described below. We will update this page accordingly.

The NOMAD computational parsers can be found within your local NOMAD distribution under
`<path_to_nomad-lab>/dependencies/parsers/<parserproject>`, where `<parser_project>` corresponds to the following organization:

* [atomistic](https://github.com/nomad-coe/atomistic-parsers) - Parsers for output from classical molecular simulations, e.g., from Gromacs, Lammps, etc.
* [database](https://github.com/nomad-coe/database-parsers) - Parsers for various databases, e.g., OpenKim.
* [eelsdb](https://github.com/nomad-coe/nomad-parser-eelsdb) - Parser for the EELS database (https://eelsdb.eu/; to be integrated in the database project).
* [electronic](https://github.com/nomad-coe/electronic-parsers) - Parsers for output from electronic structure calculations, e.g., from Vasp, Fhiaims, etc. <!-- TODO ab Initio instead of electronic structure?  -->
* [nexus](https://github.com/nomad-coe/nomad-parser-nexus) - Parsers for combining various instrument output formats and electronic lab notebooks.
* [workflow](https://github.com/nomad-coe/workflow-parsers) - Parsers for output from task managers and workflow schedulers.

You can also examine the source code of the parsers by following the above links to the corresponding GitHub repository for each project. Within each project folder you will find a `test/` directory, containing the [parser tests](../../howto/plugins/plugins.md#testing), and also a directory containing the parsers' source code,
`<parserproject>parser` or `<parserproject>parsers`, depending on if one or more
parsers are contained within the project, respectively. In the case of multiple parsers, the files
for individual parsers are contained within a corresponding subdirectory: `<parserproject>parsers/<parsername>`
For example, the Quantum Espresso parser files are found in `dependencies/parsers/electronic/electronicparsers/quantumespresso/`.


## Developing your own parser plugin

### Prerequisites

The general docs contain information about the nuts and bolts of developing a plugin:

- [How to write a plugin](../../howto/plugins/plugins.md): Some basic information about different types of plugins, plugin anatomy, and creating a plugin project.

- [How to write a parser](../../howto/plugins/parsers.md): The basics of how NOMAD parsers work, how to configure the files that your parser will match, and how to utilize existing parser classes.

!!! Attention
    This page is under construction as we convert NOMAD's standard computational parsers to parser plugins. Along the way, we will add content below to guide you in the development of your own computational parser plugins.

<!-- TODO Add best practices + tips for parser implementation/design -->
<!-- ### Best practices for computational parser design -->

<!-- ### Tips for implementation of computaional parsers -->





