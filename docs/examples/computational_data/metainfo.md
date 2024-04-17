# Guide to Computational MetaInfo

<!-- TODO replace everything below with the new schema description for data, and move to the simulation schema for data repo, then simply place a link here -->

### Overview of metadata organization for computation

NOMAD stores all processed data in a well defined, structured, and machine readable format, known as the `archive`.
The schema that defines the organization of (meta)data within the archive is known as the [MetaInfo](../../reference/glossary.md#metainfo). See [Explanation > Data structure](../../explanation/data.md) for general information about data structures and schemas in NOMAD.

The following diagram is an overarching visualization of the most important archive sections for computational data:

```
archive
├── run
│    ├── method
│    │      ├── atom_parameters
│    │      ├── dft
│    │      ├── forcefield
│    │      └── ...
│    ├── system
│    │      ├── atoms
│    │      │     ├── positions
│    │      │     ├── lattice_vectors
│    │      │     └── ...
│    │      └── ...
│    └── calculation
│           ├── energy
│           ├── forces
│           └── ...
└── workflow2
     ├── method
     ├── inputs
     ├── tasks
     ├── outputs
     └── results
```

The most important section of the archive for computational data is the `run` section, which is
divided into three main subsections: `method`, `system`, and `calculation`. `method` stores
information about the computational model used to perform the calculation.
`system` stores attributes of the atoms involved in the calculation, e.g., atom types, positions, lattice vectors, etc. `calculation` stores the output of the calculation, e.g., energy, forces, etc.
<!-- TODO Comment from ND - I would highlight the semantics of each section, since this is the main information to be communicated quickly. e.g. -->

The `workflow` section of the archive then stores information about the series of tasks performed
to accumulate the (meta)data in the run section. The relevant input parameters for the workflow are
stored in `method`, while the `results` section stores output from the workflow beyond observables
of single configurations.
For example, any ensemble-averaged quantity from a molecular dynamics
simulation would be stored under `workflow/results`. Then, the `inputs`, `outputs`, and `tasks` sections define the specifics of the workflow.
For some standard workflows, e.g., geometry optimization and molecular dynamics, the NOMAD [normalizers](../../explanation/processing.md#normalizing)
For non-standard workflows, the parser (or more appropriately the corresponding normalizer) must
populate these sections accordingly.
See [Standard and Custom Computational Workflows in NOMAD](./workflows.md) for more information about the structure of the workflow section, as well as instructions on how to upload custom workflows to link individual Entries in NOMAD.
<!-- TODO Comment from ND - Wouldn't it be easier to say that its subsection reference other sections in run? I think this better summarizes the general rule. -->
<!-- TODO add graph showing how inputs, outputs, and tasks are connected  -->
<!-- TODO add reference page of standard computational workflows and link to the above sentence. -->
<!-- TODO Link to normalizer docs will automatically populate these specifics. The parser must only create the appropriate workflow section.  -->
<!-- TODO Should give an example somewhere -->
<!-- TODO specify which workflow sections have to be set by the parser: workflow2 or these standard workflows. -->

!!! warning "Attention"
    We are currently performing a complete refactoring of the computational MetaInfo schema. Details and updates about this task, and how it may benefit your future usage of NOMAD, will be added below.

<!-- TODO Start adding the description of the Data Schema here -->