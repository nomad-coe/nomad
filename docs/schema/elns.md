# Schemas for ELNs

A schema defines all possible data structures. With small editions to our schemas,
we can instruct NOMAD to provide respective editors for data. This allows us
to build Electronic Lab Notebooks (ELNs) as tools to acquire data in a formal and
structured way. For schemas with ELN annotations, users can create new entries
in NOMAD GUI and edit the archive (structured data) of these entries directly in the GUI.

## Annotations

Definitions in a schema can have annotations. These annotations provide additional information that NOMAD can use to alter its behavior around these definitions. Annotations
are named blocks of key-value pairs:

```yaml
definitions:
  sections:
    MyAnnotatedSection:
      m_annotations:
        annotation_name:
          key1: value
          key2: value
```

Many annotations control the representation of data in the GUI. This can be for plots
or data entry/editing capabilities. There are three main categories of annotations
relevant to ELNs.

- `eln` annotations shape the data editor, i.e. allow you to control which type of forms to use to edit quantities.
- With `tabular` annotation data from linked `.csv` or `Excel` files can be parsed and added
to data.
- `plot` annotation allows you to plot numerical data (e.g. added via tables) directly
in the ELN (or)

## Example ELN
The is the commented ELN schema from our ELN example upload that can be created from
NOMAD's upload page:
```yaml
--8<-- "examples/data/eln/schema.archive.yaml"
```

{{ pydantic_model('nomad.datamodel.metainfo.annotations.ELNAnnotation', heading='## ELN Annotation') }}

{{ pydantic_model('nomad.datamodel.metainfo.annotations.BrowserAnnotation', heading='## Browser Annotation') }}


## Tabular Annotations
In order to import your data from a `.csv` or `Excel` file, NOMAD provides three distinct (and separate) ways, that
with each comes unique options for importing and interacting with your data. To better understand how to use
NOMAD parsers to import your data, three commented sample schemas are presented below. Also, each section follows
and extends a general example explained thereafter.

Two main components of any tabular parser schema are:
1) implementing the correct base-section(s), and
2) providing a `data_file` `Quantity` with the correct `m_annotations` (only exception for the entry mode). 

Please bear in mind that the schema files should 1) follow the NOMAD naming convention
(i.e. `My_Name.archive.yaml`), and 2) be accompanied by your data file in order for NOMAD to parse them.
In the examples provided below, an `Excel` file is assumed to contain all the data, as both NOMAD and
`Excel` support multiple-sheets data manipulations and imports. Note that the `Excel` file name in each schema
should match the name of the `Excel` data file, which in case of using a `.csv` data file, it can be replaced by the
`.csv` file name.

`TableData` (and any other section(s) that is inheriting from `TableData`) has a customizable checkbox Quantity 
(i.e. `fill_archive_from_datafile`) to turn the tabular parser `on` or `off`.
If you do not want to have the parser running everytime you make a change to your archive data, it is achievable then via
unchecking the checkbox. It is customizable in the sense that if you do not wish to see this checkbox at all,
you can configure the `hide` parameter of the section's `m_annotations` to hide the checkbox. This in turn sets 
the parser to run everytime you save your archive. 

Be cautious though! Turning on the tabular parser (or checking the box) on saving your data will cause
losing/overwriting your manually-entered data by the parser! 

#### Column-mode Sample:
The following sample schema creates one quantity off the entire column of an excel file (`column mode`).
For example, suppose in an excel sheet, several rows contain information of a chemical product (e.g. `purity` in one
column). In order to list all the purities under the column `purity` and import them into NOMAD, you can use the
following schema by substituting `My_Quantity` with any name of your choice (e.g. `Purity`),
`tabular-parser.data.xlsx` with the name of the `csv/excel` file where the data lies, and `My_Sheet/My_Column` with
sheet_name/column_name of your targeted data. The `Tabular_Parser` can also be changed to any arbitrary name of your 
choice.

Important notes:

- `shape: ['*']` under `My_Quantity` is essential to parse the entire column of the data file.
- The `data_file` Quantity can have any arbitrary name (e.g. `xlsx_file`) and can be referenced within the `tabular_parser`
annotation of other sections which are of type `TableData` via `path_to_data_file` in  (please see [Tabular Parser](#tabular-parser) section)
- `My_Quantity` can also be defined within another subsection (see next sample schema)

```yaml
--8<-- "examples/data/docs/tabular-parser-col-mode.archive.yaml"
```
#### Row-mode Sample:
The sample schema provided below, creates separate instances of a repeated section from each row of an excel file
(`row mode`). For example, suppose in an excel sheet, you have the information for a chemical product
(e.g. `name` in one column), and each row contains one entry of the aforementioned chemical product.
Since each row is separate from others, in order to create instances of the same product out of all rows
and import them into NOMAD, you can use the following schema by substituting `My_Subsection`,
`My_Section` and `My_Quantity` with any appropriate name (e.g. `Substance`, `Chemical_product`
and `Name` respectively).

Important notes:

- This schema demonstrates how to import data within a subsection of another subsection, meaning the
targeted quantity should not necessarily go into the main `quantites`.
- Setting `mode` to `row` signals that for each row in the sheet_name (provided in `My_Quantity`),
one instance of the corresponding (sub-)section (in this example, `My_Subsection` sub-section as it has the `repeats`
option set to true), will be appended. Please bear in mind that if this mode is selected, then all other quantities
should exist in the same sheet_name.

```yaml
--8<-- "examples/data/docs/tabular-parser-row-mode.archive.yaml"
```
#### Entry-mode Sample:
The following sample schema creates one entry for each row of an excel file (`entry mode`).
For example, suppose in an excel sheet, you have the information for a chemical product (e.g. `name` in one column),
and each row contains one entry of the aforementioned chemical product. Since each row is separate from others, in
order to create multiple archives of the same product out of all rows and import them into NOMAD, you can use the
following schema by substituting `My_Quantity` with any appropriate name (e.g. `Name`).

Important note:

- For entry mode, the convention for reading data from csv/excel file is to provide only the column name and the
data are assumed to exist in the first sheet

```yaml
--8<-- "examples/data/docs/tabular-parser-entry-mode.archive.yaml"
```

Here are all parameters for the two annotations `Tabular Parser` and `Tabular`.

{{ pydantic_model('nomad.datamodel.metainfo.annotations.TabularParserAnnotation', heading='### Tabular Parser') }}
{{ pydantic_model('nomad.datamodel.metainfo.annotations.TabularAnnotation', heading='### Tabular') }}

{{ pydantic_model('nomad.datamodel.metainfo.annotations.PlotAnnotation', heading='## Plot Annotation') }}

## Built-in base sections for ELNs

Coming soon ...

## Custom normalizers

For custom schemas, you might want to add custom normalizers. All files are parsed
and normalized when they are uploaded or changed. The NOMAD metainfo Python interface
allows you to add functions that are called when your data is normalized.

Here is an example:

```python
--8<-- "examples/archive/custom_schema.py"
```

To add a `normalize` function, your section has to inherit from `ArchiveSection` which
provides the base for this functionality. Now you can overwrite the `normalize` function
and add you own behavior. Make sure to call the `super` implementation properly to
support schemas with multiple inheritance.

If we parse an archive like this:

```yaml
--8<-- "examples/archive/custom_data.archive.yaml"
```

we will get a final normalized archive that contains our data like this:

```json
{
  "data": {
    "m_def": "examples.archive.custom_schema.SampleDatabase",
    "samples": [
      {
        "added_date": "2022-06-18T00:00:00+00:00",
        "formula": "NaCl",
        "sample_id": "2022-06-18 00:00:00+00:00--NaCl"
      }
    ]
  }
}
```

## Third-party integration

NOMAD offers integration with third-party ELN providers, simplifying the process of connecting
and interacting with external platforms. Three main external ELN solutions that are integrated into NOMAD
are: [elabFTW](https://www.elabftw.net/), [Labfolder](https://labfolder.com/) and [chemotion](https://chemotion.net/). 
The process of data retrieval and data mapping onto NOMAD's schema
varies for each of these third-party ELN provider as they inherently allow for certain ways of communicating with their
database. Below you can find a <b>How-to</b> guide on importing your data from each of these external
repositories.


### elabFTW integration

elabFTW is part of [the ELN Consortium](https://github.com/TheELNConsortium)
and supports exporting experimental data in ELN file format. ELNFileFormat is a zipped file
that contains <b>metadata</b> of your elabFTW project along with all other associated data of
your experiments.

<b>How to import elabFTW data into NOMAD:</b>

Go to your elabFTW experiment and export your project as `ELN Archive`. Save the file to your filesystem under
your preferred name and location (keep the `.eln` extension intact).
To parse your ebalFTW data into NOMAD,
go to the upload page of NOMAD and create a new upload. In the `overview` page, upload your exported file (either by 
drag-dropping it into the <i>click or drop files</i> box or by navigating to the path where you stored the file).
This causes triggering NOMAD's parser to create as many new entries in this upload as there are experiments in your
elabFTW project.

You can inspect the parsed data of each of your entries (experiments) by going to the <b>DATA</b> 
tab of each entry page. Under <i>Entry</i> column, click on <i>data</i> section. Now a new lane titled 
`ElabFTW Project Import` should be visible. Under this section, (some of) the metadata of your project is listed.
There two sub-sections: 1) <b>experiment_data</b>, and 2) <b>experiment_files</b>.

<b>experiment_data</b> section contains detailed information of the given elabFTW experiment, such as 
links to external resources and extra fields. <b>experiment_files</b> section is a list of sub-sections
containing metadata and additional info of the files associated with the experiment.


### Labfolder integration

Labfolder provides API endpoints to interact with your ELN data. NOMAD makes API calls to
retrieve, parse and map the data from your Labfolder instacne/database to a NOMAD's schema.
To do so, the necessary information are listed in the table below:

<i>project_url</i>:
        The URL address to the Labfolder project. it should follow this pattern:
        'https://your-labfolder-server/eln/notebook#?projectIds=your-project-id'. This is used to setup
        the server and initialize the NOMAD schema.

<i>labfolder_email</i>:
        The email (user credential) to authenticate and login the user. <b>Important Note</b>: this
        information <b>is discarded</b> once the authentication process is finished.

<i>password</i>:
        The password (user credential) to authenticate and login the user. <b>Important Note</b>: this
        information <b>is discarded</b> once the authentication process is finished.

<b>How to import Labfolder data into NOMAD:</b>

To get your data transferred to NOMAD, first go to NOMAD's upload page and create a new upload.
Then click on `CREATE ENTRY` button. Select a name for your entry and pick `Labfolder Project Import` from
the `Built-in schema` dropdown menu. Then click on `CREATE`. This creates an entry where you can
insert your user information. Fill the `Project url`, `Labfolder email` and `password` fields. Once completed,
click on the `save icon` in the
top-right corner of the screen. This triggers NOMAD's parser to populate the schema of current ELN.
Now the metadata and all files of your Labfolder project should be populated in this entry.

The `elements` section lists all the data and files in your projects. There are 6 main data types
returned by Labfolder's API: `DATA`, `FILE`, `IMAGE`, `TABLE`, `TEXT` and `WELLPLATE`. `DATA` element is
a special Labfolder element where the data is structured in JSON format. Every data element in NOMAD has a special
`Quantity` called `labfolder_data` which is a flattened and aggregated version of the data content.
`IMAGE` element contains information of any image stored in your Labfolder project. `TEXT` element
contains data of any text field in your Labfodler project.

### Chemotion integration

Coming soon