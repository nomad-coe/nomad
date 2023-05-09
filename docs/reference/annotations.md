# Schema annotations

Definitions in a schema can have annotations. These annotations provide additional information that NOMAD can use to alter its behavior around these definitions. Annotations are named blocks of key-value pairs:

```yaml
definitions:
  sections:
    MyAnnotatedSection:
      m_annotations:
        annotation_name:
          key1: value
          key2: value
```

Many annotations control the representation of data in the GUI. This can be for plots or data entry/editing capabilities.

{{ pydantic_model('nomad.datamodel.metainfo.annotations.ELNAnnotation', heading='## eln') }}

## Tabular data
{{ pydantic_model('nomad.datamodel.metainfo.annotations.TabularParserAnnotation', heading='### tabular_parser') }}
{{ pydantic_model('nomad.datamodel.metainfo.annotations.TabularAnnotation', heading='### tabular') }}

{{ pydantic_model('nomad.datamodel.metainfo.annotations.PlotAnnotation', heading='## plot') }}

{{ pydantic_model('nomad.datamodel.metainfo.annotations.BrowserAnnotation', heading='## browser') }}
