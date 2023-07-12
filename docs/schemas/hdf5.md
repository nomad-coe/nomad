---
title: Large data
---
The NOMAD schemas and processed data system is designed to describe and manage
intricate hierarchies of connected data. This is ideal for metadata and lots of small
data quantities, but does not work for large quantities. Quantities are atomic and
are always manages as a whole; there is currently no functionality to stream or
splice large quantities. Consequently, tools that produce or work with such data
cannot scale.

## HDF5Reference

Large data quantities should be managed within HDF5 raw files. These can be files
that are either produced during processing (e.g. created by a parser) or that are
uploaded and contain the large quantities already. To describe these data, schemas
and processed data can include references into HDF5 raw files. This allows to
describe large data quantities with metadata that is described with schemas and
contained in the processed data.

The data type `HDF5Reference` is a regular type that can be used for quantities
in schemas. The values are similar to [reference values](basics.md#different-forms-of-references) and contain
the HDF5 file and a path to a group or field in the HDF5, e.g. `../upload/raw/large_data.hdf5#group/large_field`.

For HDF5 files that are generated during processing, it is good practice to structure
the HDF5 inline with the schema. E.g. a `HDF5Reference` quantity with a sub-section
path `data.process.log` should be stored in a field names `log` in the sub-group `process`
that is part of the root group `data`.

In future version of NOMAD, processed data might be maintained partially or in full in
HDF5 files. Structuring HDF5 files and processed data alike, might simplify later migration.

NOMAD clients (e.g. NOMAD UI) can pick up on these `HDF5Reference` quantities and
provide respective functionality (e.g. showing a H5Web view).

!!! attention

    This part of the documentation is still work in progress.

## Metadata for large quantities

!!! attention

    This will be implemented and documented soon.