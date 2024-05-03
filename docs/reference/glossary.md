# Glossary

<!--
TODO consider the following items:
- Coauthor
- Owner (if this really differs from author)
- Reviewer
- Oasis
- Package (metainfo)
- Workflow
- App (maybe too soon to add since we don't know yet what this really entails)
-->

This is a list of terms that have a specific meaning for NOMAD and are used through
out the application and this documentation.

### Annotation

*Annotations* are part of data [schemas](#schema) and they describe aspects that are not
directly defining the type or shape of data. They often allow to alter how certain data is
managed, represented, or edited. See [annotations in the schema documentation](../howto/customization/elns.md#annotations).

### App

Apps allow you to build customized user interfaces for specific research domains, making it easier to navigate and understand the data. This typically means that certain domain-specific properties are highlighted, different units may be used for physical properties, and specialized dashboards may be presented. This becomes crucial for NOMAD installations to be able to scale with data that contains a mixture of experiments and simulations, different techniques, and physical properties spanning different time and length scales.

### Archive

NOMAD processes (parses and normalizes) all data.
The entirety of all [processed data](#processed-data) is referred to as the
*Archive*.
Sometimes the term *archive* is used to refer to the processed data of a
particular [entry](#entry), e.g. "the archive of that entry".

The term *archive* is an old synonym for [processed data](#processed-data). Since the term
has different meaning for different people and is very abstract, we are slowly deprecating
its use in favor of *processed data* (e.g. NOMAD Repository and NOMAD Archive).

### Author

An *author* is typically a natural person that has uploaded a piece of data into NOMAD and
has authorship over it. Often *authors* are [users](#user), but not always.
Therefore, we have to distinguish between authors and users.

### ELN

Electronic Lab Notebooks (*ELNs*) are a specific kind of [entry](#entry) in NOMAD. These
entries can be edited in NOMAD, in contrast to entries that are created by uploading
and processing data. ELNs offer form fields and other widgets to modify the contents of
an entry. As all entries, *ELNs* are based on a [schema](#schema); how [quantities](#quantity)
are edited (e.g. which type of widget) can be controlled through [annotations](#annotation).

### Entry

Data in NOMAD is organized in *entries* (as in "database *entry*"). Entries have an
*entry id*. Entries can be searched for and entries have individual pages on the NOMAD GUI. Entries are always
associated with [raw files](#raw-file), where one of these files is the [mainfile](#mainfile).
Raw files are processed to create the [processed data](#processed-data) (or the [archive](#archive))
for an entry.

### Dataset

Users can organize [entries](#entry) into *datasets*. Datasets are not created automatically,
don't confuse them with [uploads](#upload). Datasets can be compared to albums, labels, or tags
on other platforms. Datasets are used to reference a collection of data and users can get a DOI for their
datasets.

### Mainfile

Each [entry](#entry) has one [raw file](#raw-file) that defines it. This is called the
*mainfile* of that entry. Typically most, if not all, [processed data](#processed-data)
of an entry is retrieved from that mainfile.

### Metadata

In NOMAD *metadata* refers to a specific technical sub-set of [processed data](#processed-data).
The metadata of an [entry](#entry) comprises ids, timestamps, hashes, authors, datasets,
references, used schema, and other information.

### Metainfo

The term *metainfo* refers to the sum of all [schemas](#schema). In particular it is
associated with all pre-defined schemas that are used to represent all
[processed data](#processed-data) in a standardized way. Similar to an *ontology*,
the metainfo provides additional meaning by associated in each piece of data with
name, description, categories, type, shape, units, and more.

### Normalizer

A *normalizer* is a small tool that can refine the [processed data](#processed-data) of an
[entry](#entry). Normalizers can read and modify processed data and thereby either normalize
(change) some of the data or add normalized derived data. Normalizers are run after
[parsers](#parser) and are often used to do processing steps that need to be applied to
the outcome of many parsers and are therefore not part of the parsers themselves.

There are normalizer *classes* and normalize *functions*. The normalizer classes are
run after parsing in a particular order and if certain conditions are fulfilled.
Normalize functions are part of [schemas](#schema) (i.e. [section definitions](#section-and-subsection)).
They are run at the end of processing on all the sections that instantiate the respective
section definition.

### Parser

A *parser* is a small program that takes a [mainfile](#mainfile) as input and produces
[processed data](#processed-data). Parsers transform information from a particular source
format into NOMAD's structured schema-based format. Parsers start with a mainfile, but
can open and read data from other files (e.g. those referenced in the mainfile). Typically,
a parser is associated with a certain file-format and is only applied to files of that
format.

### Plugin

NOMAD installations can be customized through plugins, which are Git repositories containing an installable python package that will add new features upon being installed. Plugins can contain one or many plugin entry points, which represent individual customizations.

### Plugin entry point

Plugin entry points are used to configure and load different types of NOMAD customizations. There are several entry point types, including entry points for parsers, schema packages and apps. A single plugin may contain multiple entry points.

### Processed data

NOMAD processes (parses and normalizes) all data. The *processed data* is the outcome of this process.
Therefore, each NOMAD [entry](#entry) is associated with *processed data* that contains all the parsed and normalized
information. Processed data always follows a [schema](#schema). Processed data can be retrieved (via API) or downloaded as `.json` data.

### Processing

NOMAD processes (parses and normalizes) all data. During processing, all provided files
are considered. First, files are matched to [parsers](#parser). Second, files that
match with a parser, become [mainfiles](#mainfile) and an [entry](#entry) is created.
Third, we run the parser to create [processed data](#processed-data). Fourth, the
processed data is further refined by running [normalizers](#normalizer). Last, the
processed data is saved and indexed. The exact processing time depends on the size of the
uploaded data and users can track the processing state of each entry in the GUI.

### Quantity

All [processed data](#processed-data) is structured into sections and quantities. Sections
provide hierarchy and organization, *quantities* refer to the actual pieces of data.
In NOMAD, a quantity is the smallest referable unit of [processed data](#processed-data).
Quantities can have many types and shapes; examples are strings, numbers, lists, or matrices.

In a [schema](#schema), quantities are defined by their name, description, type, shape, and unit.
Quantities in processed data are associated with a respective quantity definition from
the respective schema.

### Raw file

A *raw file* is any file that was provided by a NOMAD user. A raw-file might produce an
[entry](#entry), if it is of a supported file-format, but does not have to. Raw files
always belong to an [upload](#upload) and might be associated with an [entry](#entry)
(in this case, raw-files are also [mainfiles](#mainfile)).

The sum of all raw files is also referred to as the *Repository*. This is an old term
from when NOMAD was implemented as several services (e.g. NOMAD Repository and NOMAD Archive).

### Results (section `results`)

The *results* are a particular section of [processed data](#processed-data). They comprise
a summary of the most relevant data for an [entry](#entry).

While all [processed data](#processed-data) can be downloaded and is accessible via API,
an entry's *results* (combined with its [metadata](#metadata)) is also searchable and can
be read quicker and in larger amounts.

### Schema

*Schemas* define possible data structures for [processed data](#processed-data). Like
a book they organize data hierarchically in *sections* and *subsections*. Schemas
are similar to *ontologies* as they define possible relationships between data organized
within them.

A schema is a collection of [section](#section-and-subsection) and [quantity](#quantity)
definitions. Schemas are organized in [schema packages](#schema-package), i.e. collections of definitions. All schemas combined form the [metainfo](#metainfo).

### Schema package

*Schema packages* contain a collection of [schema](#schema) definitions. Schema packages may be defined as [YAML files](../howto/customization/basics.md) or in Python as [plugin entry points](../howto/plugins/schema_packages.md).

### Section and Subsection

All [processed data](#processed-data) is structured into sections and quantities. *Sections*
provide hierarchy and organization, *quantities* refer to the actual pieces of data.

In a [schema](#schema), a section are defined by their name, description, all possible
*subsections*, and [quantities](#quantity). Section definitions can also inherit
all properties (subsections, quantities) from other section definitions using them as
*base sections*.

### Upload

NOMAD organizes [raw-files](#raw-file) (and all [entries](#entry) created from them)
in *uploads*. Uploads consist of a directory structure of raw-files and a list of
respective entries.

Uploads are created by a single [user](#user), the *owner*. Uploads have two states.
Initially, they are mutable and have limited visibility. Owners can invite other to
collaborate and those users can add/remove/change data.
The owner can publish an upload at some point, where the upload becomes immutable and visible
to everyone. Uploads are the smallest unit of data that can be individually shared and published.

### User

A *user* is anyone with a NOMAD account. It is different from an [author](#author) as
all users can be authors, but not all authors have to be users. All data in NOMAD is always
owned by a single user (others can be collaborators and co-authors).