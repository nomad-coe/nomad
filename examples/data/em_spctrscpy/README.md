# An Electronic Lab Notebook (ELN) Example for Electron Microscopy (EM)

## Introduction

This example shows how the NOMAD ELN functionalities can be used to collect
required metadata for groups, fields, and attributes for creating a dataset
which is compliant with the NeXus NXem application definition.

Specifically, this example shows how I/O functionalities of hyperspy can be used
to load data of spectroscopy experiments from three exemplar file formats, namely
Bruker BCF, Velox EMD, and DigitalMicrograph DM3 into NOMAD OASIS. The implementation
shows how instances of NeXus base classes like NXspectrum_set_em_xray,
NXspectrum_set_em_eels, and NXimage_set_em_adf can be created and registered inside
instances of NXevent_data_em. The example shows how all these base classes can be
composed and stored inside a NeXus/HDF5 file.

This makes the example relevant for researchers who work within the fields of
scanning electron microscopy (SEM) EDS/EDX, transmission electron microscopy
(EDS/STEM), and electron energy loss spectroscopy (EELS).

The NeXus NXem data model is documented here [NXem](https://fairmat-experimental.github.io/nexus-fairmat-proposal)

This example also shows how the NOMAD OASIS ELN functionalities can be used
to supplement a NeXus file with metadata which are currently or usually not stored
in vendor or community files including details about the specimen, the instrument,
users, and post-processing steps.

This example upload contains the following entries:
- A schema in NOMAD's *archive.yaml* format: *nxem.schema.archive.yaml*
- A schema instance file used by NOMAD *nxem.archive.json* which is filled for educational purpose with values for the example.
- The primary consumer of this json file is NOMAD and its internal data management system.
- Another schema instance file used by the nomad-parser-nexus *eln_data.yaml*. This file contains all entered
quantities from the ELN GUI (after the save button was stored). The example is also filled for educational purposes
with values matching those in nxem.archive.json.
Files are updated each time the save button in the ELN GUI is clicked.
The eln_data.yaml file is used by the [NOMAD-PARSER-NEXUS](https://github.com/nomad-coe/nomad-parser-nexus).

This example is configured to take an example dataset and call the nomad-parser-nexus dataconverter
which creates a NeXus file compliant with NXem. Once completed, this file is available in the
upload section/staging area. This makes also these files explorable with H5Web visualization
through the files menu.

This example comes with a measured datasets which are meant for testing and exploring.
The datasets include two examples kindly shared by Adrien Teurtrie and Cécile Hebert at EPFL
(Bruker BCF, Velox EMD), and a Digital Micrograph DM3 file with EELS data measured and kindly
shared by Hannah Nerl and Christoph Koch.

## Creating NeXus files

When you modify the ELN and click the save button, the data from the ELN will be
parsed and combined with the content from the vendor file to create the NXS file.
You can replace these files with your own and accordingly use the ELN to enter your
own metadata. Upon saving, a NeXus/HDF5 file in NXapm format will be created for your specific dataset.

With this functionality, you can use this example as a template to translate your own
datasets into NeXus. The drag-and-drop functionality of the upload section can be
used to pass your vendor file onto the respective file upload fields of the ELN.
After clicking the save button, the newly entered metadata and files will be processed
on-the-fly and a new NeXus file, compliant with NXem will be generated.

## Where to go from now

With an example in your upload **you can explore** the content in H5Web.
Furthermore, you can work with the data by starting for instance
a generic jupyterlab container via the **Analytics tab in the** NOMAD OASIS
menu bar. This container is a part of the Nomad Remote Tools Hub (NORTH) service.

**Once running, the container offers a jupyter-lab server and notebook.**

## Summary

The example is meant as a starting point. You can download the schema file and extend the
schema to collect further metadata (e. g. for adding optional quantities defined in NXem) based
on what is relevant for your laboratory and use case. Also you can explore the implementation
of the example to customize it for your own needs. We would be happy if you could support us
with improving this example and the associated schemes by leaving us comments via the 
nexus-fairmat-proposal pages or by contacting us via the various channels.
We are also very happy to guide you on how to customize these functionalities
for your own laboratory and needs.

Consult our [documentation on the NOMAD Archive and Metainfo](https://nomad-lab.eu/prod/v1/docs/archive.html)
to learn more about schemes.

## Questions, comments, suggestions?

For general questions regarding the EM tools and if you're interested in building one for your
own research workflow or your colleagues and group you are very welcome to contact
[Markus Kühbach](https://www.fair-di.eu/fairmat/fairmat_/fairmatteam) from the FAIRmat consortium.


## Known bugs

