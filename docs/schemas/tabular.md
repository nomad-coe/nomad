In order to import your data from a `.csv` or `Excel` file, NOMAD provides three distinct (and separate) ways, that
with each comes unique options for importing and interacting with your data. In order to better understand how to use
NOMAD tabular parser to import your data, follow three sections below. In each section you
can find a commented sample schema with a step-by-step guide on how to import your tabular data.

Tabular parser, implicitly, parse the data into the same NOMAD entry where the datafile is loaded. Also, explicitly,
this can be defined by putting the corresponding annotations under `current_entry` (check the examples below).
In addition, tabular parser can be set to parse the data into new entry (or entries). For this, the proper annotations
should be appended to `new_entry` annotation in your schema file.

Two main components of any tabular parser schema are:
1) implementing the correct base-section(s), and
2) providing a `data_file` `Quantity` with the correct `m_annotations`.

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

## Column-mode
The following sample schema creates one quantity off the entire column of an excel file (`column mode`).
For example, suppose in an excel sheet, several rows contain information of a chemical product (e.g. `purity` in one
column). In order to list all the purities under the column `purity` and import them into NOMAD, you can use the
following schema by substituting `My_Quantity` with any name of your choice (e.g. `Purity`),
`tabular-parser.data.xlsx` with the name of the `csv/excel` file where the data lies, and `My_Sheet/My_Column` with
sheet_name/column_name of your targeted data. The `Tabular_Parser` can also be changed to any arbitrary name of your
choice.

Important notes:

- `shape: ['*']` under `My_Quantity` is essential to parse the entire column of the data file.
- The `data_file` `Quantity` can have any arbitrary name (e.g. `xlsx_file`)
- `My_Quantity` can also be defined within another subsection (see next sample schema)
- Use `current_entry` and append `column_to_sections` to specify which sub_section(s) is to be filled in
this mode. `Leaving this field empty` causes the parser to parse the entire schema under column mode.

```yaml
--8<-- "examples/data/docs/tabular-parser-col-mode.archive.yaml"
```

<b>Step-by-step guide to import your data using column-mode:</b>

After writing your schema file, you can create a new upload in NOMAD (or use an existing upload),
and upload both your `schema file` and the `excel/csv` file together (or zipped) to your NOMAD project. In the
`Overview` page of your NOMAD upload, you should be able to see a new entry created and appended to the `Process data`
section. Go to the entry page, click on `DATA` tab (on top of the screen) and in the `Entry` lane, your data
is populated under the `data` sub_section.

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
- Setting `row_to_sections` under `current_entry` signals that for each row in the sheet_name (provided in `My_Quantity`),
one instance of the corresponding (sub-)section (in this example, `My_Subsection` sub-section as it has the `repeats`
option set to true), will be appended. Please bear in mind that if this mode is selected, then all other quantities
in this sub_section, should exist in the same sheet_name.


```yaml
--8<-- "examples/data/docs/tabular-parser-row-mode.archive.yaml"
```

<b>Step-by-step guide to import your data using row-mode:</b>

After writing your schema file, you can create a new upload in NOMAD (or use an existing upload),
and upload both your `schema file` and the `excel/csv` file together (or zipped) to your NOMAD project. In the
`Overview` page of your NOMAD upload, you should be able to see as many new sub-sections created and appended
to the repeating section as there are rows in your `excel/csv` file.
Go to the entry page of the new entries, click on `DATA` tab (on top of the screen) and in the `Entry` lane,
your data is populated under the `data` sub_section.

#### Entry-mode Sample:
The following sample schema creates one entry for each row of an excel file (`entry mode`).
For example, suppose in an excel sheet, you have the information for a chemical product (e.g. `name` in one column),
and each row contains one entry of the aforementioned chemical product. Since each row is separate from others, in
order to create multiple archives of the same product out of all rows and import them into NOMAD, you can use the
following schema by substituting `My_Quantity` with any appropriate name (e.g. `Name`).

Important note:

- To create new entries based on your entire schema, set `row_to_entries` to `- root`. Otherwise, you can
provide the relative path of specific sub_section(s) in your schema to create new entries.
- Leaving `row_to_entries` empty causes the parser to parse the entire schema using <b>column mode</b>!


```yaml
--8<-- "examples/data/docs/tabular-parser-entry-mode.archive.yaml"
```

<b>Step-by-step guide to import your data using entry-mode:</b>

After writing your schema file, you can create a new upload in NOMAD (or use an existing upload),
and upload both your `schema file` and the `excel/csv` file together (or zipped) to your NOMAD project. In the
`Overview` page of your NOMAD upload, you should be able to see as many new entries created and appended
to the `Process data` section as there are rows in your `excel/csv` file.
Go to the entry page of the new entries, click on `DATA` tab (on top of the screen) and in the `Entry` lane,
your data is populated under the `data` sub_section.

<b>Advanced options to use/set in tabular parser:</b>

- If you want to populate your schema from multiple `excel/csv` files, you can
define multiple data_file `Quantity`s annotated with `tabular_parser` in the root level of your schema
(root level of your schema is where you inherit from `TableData` class under `base_sections`).
Each individual data_file quantity can now contain a list of sub_sections which are expected to be filled
using one- or all of the modes mentioned above. Check the `MyOverallSchema` section in
`Complex Schema` example below. It contains 2 data_file quantities that each one, contains separate instructions
to populate different parts of the schema. `data_file_1` is responsible to fill `MyColSubsection` while `data_file_2`
fills all sub_sections listed in `row_to_sections` and `entry_to_sections` under `new_entry`.

- When using the entry mode, you can create a custom `Quantity` to hold a reference to each new entries
generated by the parser. Check the `MyEntrySubsection` section in the `Complex Schema` example below.
The `refs_quantity` is a `ReferenceEditQuantiy` with type `#/MyEntry` which tells the parser to
populate this quantity with a reference to the fresh entry of type `MyEntry`. Also, you may use
`tabular_pattern` annotation to explicitly set the name of the fresh entries.

- If you have multiple columns with exact same name in your `excel/csv` file, you can parse them using row mode.
For this, define a repeating sub_section that handles your data in different rows and inside each row, define another
repeating sub_section that contains your repeating columns. Check `MySpecialRowSubsection` section in the
`Complex Schema` example below. `data_file_2` contains a repeating column called `row_quantity_2` and
we want to create a section out of each row and each column. This is done by
creating one row of type `MySpecialRowSubsection` and populate
`MyRowQuantity3` quantity from `row_quantity_3` column in the `csv` file, and appending each column of
`row_quantity_2` to `MyRowQuantity2`.

```yaml
--8<-- "examples/data/docs/tabular-parser-complex.archive.yaml"
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
retrieve, parse and map the data from your Labfolder instance/database to a NOMAD's schema.
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

NOMAD supports importing your data from Chemotion repository via `chemotion` parser. The parser maps
your data that is structured under chemotion schema, into a predefined NOMAD schema. From your Chemotion
repo, you can export your entire data as a zip file which then is used to populate NOMAD schema.

<b>How to import Chemotion data into NOMAD:</b>

Go to your Chemotion repository and export your project. Save the file to your filesystem under
your preferred name and location (`your_file_name.zip`).
To get your data parsed into NOMAD,
go to the upload page of NOMAD and create a new upload. In the `overview` page, upload your exported file (either by
drag-dropping it into the <i>click or drop files</i> box or by navigating to the path where you stored the file).
This causes triggering NOMAD's parser to create one new entry in this upload.

You can inspect the parsed data of each of this new entry by navigating to the <b>DATA</b>
tab of the current entry page. Under <i>Entry</i> column, click on <i>data</i> section. Now a new lane titled
`Chemotion Project Import` should be visible. Under this section, (some of) the metadata of your project is listed.
Also, there are various (sub)sections which are either filled depending on whether your datafile
contains information on them.

If a section contains an image (or attachment) it is appended to the same section under `file` Quantity.
