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

{{ pydantic_model('nomad.datamodel.metainfo.annotations.ELNAnnotation', heading='## ELN annotations') }}

{{ pydantic_model('nomad.datamodel.metainfo.annotations.BrowserAnnotation', heading='## Browser') }}

{{ pydantic_model('nomad.datamodel.metainfo.annotations.BrowserAnnotation', heading='## browser') }}

### `label_quantity`

This annotation goes in the section that we want to be filled with tabular data, not in the single quantities.
It is used to give a name to the instances that might be created by the parser. If it is not provided, the name of the section itself will be used as name.
Many times it is useful because, i. e., one might want to create a bundle of instances of, say, a "Substrate" class, each instance filename not being "Substrate_1", "Substrate_2", etc., but being named after a quantity contained in the class that is, for example, the specific ID of that sample.


```yaml
MySection:
  more:
    label_quantity: my_quantity
  quantities:
    my_quantity:
      type: np.float64
      shape: ['*']
      description: "my quantity to be filled from the tabular data file"
      unit: K
      m_annotations:
        tabular:
          name: "Sheet1/my header"
        plot:
          x: timestamp
          y: ./my_quantity
```

!!! important
    The quantity designated as `label_quantity` should not be an array but a integer, float or string, to be set as the name of a file. If an array quantity is chosen, the parser would fall back to the use of the section as name.
## Tabular data

{{ pydantic_model('nomad.datamodel.metainfo.annotations.TabularAnnotation', heading='### `tabular`') }}

Each and every quantity to be filled with data from tabular data files should be annotated as the following example.
A practical example is provided in [How To](../schemas/tabular.md#preparing-the-tabular-data-file) section.

```yaml
my_quantity:
  type: np.float64
  shape: ['*']
  description: "my quantity to be filled from the tabular data file"
  unit: K
  m_annotations:
    tabular:
      name: "Sheet1/my header"
    plot:
      x: timestamp
      y: ./my_quantity
```

### `tabular_parser`

One special quantity will be dedicated to host the tabular data file. In the following examples it is called `data_file`, it contains the `tabular_parser` annotation, as shown below.

{{ pydantic_model('nomad.datamodel.metainfo.annotations.TabularParserAnnotation', heading = '') }}

### Available Combinations

|Tutorial ref.|`file_mode`|`mapping_mode`|`sections`|How to ref.|
|---|---|---|---|---|
|1|`current_entry`|`column`|`root`|[HowTo](../schemas/tabular.md#1-column-mode-current-entry-parse-to-root)|
|2|`current_entry`|`column`|my path|[HowTo](../schemas/tabular.md#2-column-mode-current-entry-parse-to-my-path)|
|<span style="color:red">np1</span>|`current_entry`|`row`|`root`|<span style="color:red">Not possible</span>|
|3|`current_entry`|`row`|my path|[HowTo](../schemas/tabular.md#3-row-mode-current-entry-parse-to-my-path)|
|<span style="color:red">np2</span>|`single_new_entry`|`column`|`root`|<span style="color:red">Not possible</span>|
|4|`single_new_entry`|`column`|my path|[HowTo](../schemas/tabular.md#4-column-mode-single-new-entry-parse-to-my-path)|
|<span style="color:red">np3</span>|`single_new_entry`|`row`|`root`|<span style="color:red">Not possible</span>|
|5|`single_new_entry`|`row`|my path|[HowTo](../schemas/tabular.md#5-row-mode-single-new-entry-parse-to-my-path)|
|<span style="color:red">np4</span>|`multiple_new_entries`|`column`|`root`|<span style="color:red">Not possible</span>|
|<span style="color:red">np5</span>|`multiple_new_entries`|`column`|my path|<span style="color:red">Not possible</span>|
|6|`multiple_new_entries`|`row`|`root`|[HowTo](../schemas/tabular.md#6-row-mode-multiple-new-entries-parse-to-root)|
|7|`multiple_new_entries`|`row`|my path|[HowTo](../schemas/tabular.md#7-row-mode-multiple-new-entries-parse-to-my-path)|

```yaml
data_file:
  type: str
  description: "the tabular data file containing data"
  m_annotations:
    tabular_parser:
      parsing_options:
        comment: '#'
      mapping_options:
      - mapping_mode: column
        file_mode: single_new_entry
        sections:
        - my_section/my_quantity
```

<!-- The available options are:

|**name**|**type**|**description**|
|---|---|---|
|`parsing_options`|group of options|some pandas `Dataframe` options.|
|`mapping_options`|list of groups of options|they allow to choose among all the possible modes of parsing data from the spreadsheet file to the NOMAD archive file. Each group of options can be repeated in a list. | -->


{{ pydantic_model('nomad.datamodel.metainfo.annotations.PlotAnnotation', heading='## Plot') }}

