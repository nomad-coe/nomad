# How to extend the search

## The search indices

NOMAD uses Elasticsearch as the underlying search engine. The respective indices
are automatically populated during processing and other NOMAD operations. The indices
are built from some of the archive information of each entry. These are mostly the
sections `metadata` (ids, user metadata, other "administrative" and "internal" metadata)
and `results` (a summary of all extracted (meta)data). However, these sections are not
indexed verbatim. What exactly and how it is indexed is determined by the Metainfo
and the `elasticsearch` Metainfo extension.

### The `elasticsearch` Metainfo extension

Here is the definition of `results.material.elements` as an example:

```python
class Material(MSection):
    ...
    elements = Quantity(
        type=MEnum(chemical_symbols),
        shape=["0..*"],
        default=[],
        description='Names of the different elements present in the structure.',
        a_elasticsearch=[
            Elasticsearch(material_type, many_all=True),
            Elasticsearch(suggestion="simple")
        ]
    )
```

Extensions are denoted with the `a_` prefix as in `a_elasticsearch`.
Since extensions can have all kinds of values, the `elasticsearch` extension is rather
complex and uses the `Elasticsearch` class.

There can be multiple values. Each `Elasticsearch` instance configures a different part
of the index. This means that the same quantity can be indexed multiple time. For example,
if you need a text- and a keyword-based search for the same data. Here
is a version of the `metadata.mainfile` definition as another example:

```python
mainfile = metainfo.Quantity(
    type=str, categories=[MongoEntryMetadata, MongoSystemMetadata],
    description='The path to the mainfile from the root directory of the uploaded files',
    a_elasticsearch=[
        Elasticsearch(_es_field='keyword'),
        Elasticsearch(
            mapping=dict(type='text', analyzer=path_analyzer.to_dict()),
            field='path', _es_field='')
    ]
)
```

### The different indices

The first (optional) argument for `Elasticsearch` determines where the data is indexed.
There are three principle places:

- the entry index (`entry_type`, default)
- the materials index (`material_type`)
- the entries within the materials index (`material_entry_type`)

#### Entry index

This is the default and is used even if another (additional) value is given. All data
is put into the entry index.

#### Materials index

This is a separate index from the entry index and contains aggregated material
information. Each document in this index represents a material. We use a hash over some
material properties (elements, system type, symmetry) to define what a material is and
which entries belong to which material.

Some parts of the material documents contain the material information that is always
the same for all entries of this material. Examples are elements, formulas, symmetry.

#### Material entries

The materials index also contains entry-specific information that allows to filter
materials for the existence of entries with certain criteria. Examples are publish status,
user metadata, used method, or property data.

### Adding quantities

In principle, all quantities could be added to the index, but for convention and
simplicity, only quantities defined in the sections `metadata` and `results` should be
added. This means that if you want to add custom quantities from your parser, for example,
you will also need to customize the results normalizer to copy or reference parsed data.

## The search API

The search API does not have to change. It automatically supports all quantities with the
`elasticsearch` extension. The keys that you can use in the API are the Metainfo paths of
the respective quantities, e.g. `results.material.elements` or `mainfile` (note that the
`metadata.` prefix is always omitted). If there are multiple `elasticsearch` annotations
for the same quantity, all but one define a `field` parameter, which is added to the
quantity path, e.g. `mainfile.path`.

## The search web interface

!!! warning "Attention"
        Coming soon ...