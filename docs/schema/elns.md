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


## ELN Annotations
The `eln` annotations can contain the following keys:

{{ get_schema_doc('eln') }}

The `eln` `component` can be one of the following components:

{{ get_schema_doc('component') }}

As part of the GUI, you'll find an overview about all
ELN edit annotations and components [here]({{ nomad_url() }}/../gui/dev/editquantity).


## Tabular Annotations
Tabular annotation accepts the following keys:

{{ get_schema_doc('tabular') }}

## Plot Annotations
Plot annotation is a wrapper for [plotly](https://plotly.com) library. One can use the following keys for plot annotation:

{{ get_schema_doc('plot') }}

which can be customized by using plotly commands. See [plot examples]({{ nomad_url() }}/../gui/dev/plot).

## Build-in base sections for ELNs

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