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