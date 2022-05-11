This is a simple example for a basic ELN. It demonstrates the use of a NOMAD schema
to define different types of entries. Based on this schema the ELN allows you to create
Samples, Chemicals, and Instruments. The Sample entry type also allows to define
processes.

The schema is meant as a starting point. You can download the schema file and
extend the schema. For your own ELNs, you have to consider the data that you acquire
in your lab and create your own schema based on the concepts shown in this example.

Consult our [documentation on the NOMAD Archive and Metainfo](https://nomad-lab.eu/prod/v1/docs/archive.html) to learn more about schemas.

This example uploads contains the following entries
- A schema in NOMAD's *archive.yaml* format: *schema.archive.yaml*
that defines Three types of ELN entries: sample, instrument, and chemical
- Three chemicals (as defined in the schema): *Copper_II_Selenide.archive.json*,
*Tin_II_Selenide.archive.json*,
*Zinc_Selenide.archive.json*
- An instrument *PVD-P*.archive.json
- A sample (*sample.archive.json*) with two processes (PVD evaporation, hotplate annealing) as sub-sections,
and references to instrument and chemicals.
- A *.csv* file. This is not directly parser by NOMAD, but the sample ELN uses it to
parse data for the PVD evaporation process.