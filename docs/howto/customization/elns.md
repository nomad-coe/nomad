# How to define and use ELNs in NOMAD

## Schemas for ELNs

A schema defines all possible data structures. With small additions to our schemas, we can instruct NOMAD to provide respective editors for data. This allows us to build Electronic Lab Notebooks (ELNs) as tools to acquire data in a formal and structured way. For schemas with ELN annotations, users can create new entries in the NOMAD repository and edit the archive (structured data) of these entries directly in the GUI.

### Annotations

Definitions in a schema can have annotations. With these annotations you can provide additional information that NOMAD can use to alter its behavior around these definitions. Annotations are named blocks of key-value pairs:

```yaml
definitions:
  sections:
    MyAnnotatedSection:
      m_annotations:
        annotation_name:
          key1: value
          key2: value
```

Many annotations control the representation of data in the GUI. This can be for plots or data entry/editing capabilities. There are three main categories of annotations relevant to ELNs. You find a reference of all annotations [here](../../reference/annotations.md).

### Example ELN schema
The is the commented ELN schema from our ELN example upload that can be created from NOMAD's upload page:
```yaml
--8<-- "examples/data/eln/schema.archive.yaml"
```

---
**NOTE: Defining Labels for Quantities**

When defining labels for quantities, utilize the
[display annotations](../../reference/annotations.md#display-annotations)
and ensure that you follow the conventions as described 
[here](./basics.md#conventions-for-labels).
---
