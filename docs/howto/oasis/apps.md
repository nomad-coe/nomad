# How to configure custom apps

Apps provide customized views of data for specific domains, making it easier for
the users to navigate and understand the data. This typically means that certain
domain-specific properties are highlighted, different units may be used for
physical properties, and specialized dashboards may be presented. This becomes
crucial for NOMAD installations to be able to scale with data that contains a
mixture of experiments and simulations, different techniques, and physical
properties spanning different time and length scales.

Apps only affect the way data is *displayed* for the user: if you wish to affect
the underlying data structure, you will need to define a custom [Python schema](../customization/plugins_dev.md#develop-a-schema-plugin)
or [YAML schema](../customization/basics.md). It is common that a custom schema has
an app associated with it, but apps can also provide different views of the same
underlying data.

Apps are defined with a static YAML configuration file, which means that no
special programming skills are needed and app definitions can be shared easily.

## App example

Here is an example of a simple app definition in YAML. A full breakdown of the
configuration options are given in the [reference below](#app-configuration-reference).

```yaml
# Label of the App
label: 'My App'
# Path used in the URL, must be unique
path: 'myapp'
# Used to categorize apps in the explore menu
category: 'Simulations'
# Brief description used in the app menu
description: 'An app customized for me.'
# Longer description that can also use markdown
readme: 'Here is a much longer description of this app'
# Controls which columns are shown in the results table
columns:
  selected:
    - 'entry_type'
  options:
    entry_type:
      label: 'Entry type'
      align: 'left'
    upload_create_time:
      label: 'Upload time'
      align: 'left'
# Dictionary of search filters that are always enabled for queries made within
# this app. This is especially important to narrow down the results to the
# wanted subset. Any available search filter can be targeted here.
filters_locked:
  upload_create_time:
    gte: 0
# Controls the filter menus shown on the left
filter_menus:
  options:
    material:
      label: 'Material'
      level: 0
    elements:
      label: 'Elements / Formula'
      level: 1
      size: 'xl'
# Controls the default dashboard shown in the search interface
dashboard:
  widgets:
  - type: histogram
    showinput: false
    autorange: true
    nbins: 30
    scale: linear
    quantity: results.material.n_elements
    layout:
      lg:
        minH: 3
        minW: 3
        h: 4
        w: 12
        y: 0
        x: 0
```

## Customizing default apps in a NOMAD installation

Each NOMAD installation has a set of built-in apps, which are controlled through
the [ui.apps](../../reference/config.md#ui) field in the `nomad.yaml` configuration file. These are
the apps that are defined by default in a NOMAD installation:

{{ default_apps_list()}}

In `nomad.yaml`, it is easy to to select which apps to include or exclude like
this:

```yaml
ui:
  apps:
    include: ['entries', 'materials']
```

It is also possible to customize specific parts of an existing app definition:

```yaml
ui:
  apps:
    options:
      entries:
        columns:
          exclude: ['upload_create_time']
```

Completely new apps can also be defined by adding new entries to the
`ui.apps.options` dictionary:

```yaml
ui:
  apps:
    options:
      myapp:
        label: 'My App'
        ...
```

If no explicit rules are added in `ui.apps.include` or `ui.apps.exclude`, these
new options will be included by default.

## Using schema quantities in an app

Each app may define additional **filters** that should be enabled in it. Filters
have a special meaning in the app context: filters are pieces of (meta)info than
can be queried in the search interface of the app, but also targeted in the rest
of the app configuration as explained below in. By default, none of the
quantities from custom schemas are loaded as filters, and they need to be
explicitly added.

!!! note

    Note that not all of the quantities from a custom schema can be exposed as
    filters. At the moment we only support targeting **scalar** quantities from
    custom schemas.

Each schema has a unique name within the NOMAD ecosystem, which is needed to
target them in the configuration. The name depends on the resource in which the
schema is defined in:

- Python schemas are identified by the python path for the class that inherits
from `EntryData`. For example, if you have a python package called `myschema`,
which has a module called `schema.py`, which contains the class `MySchema`, then
the schema name will be `myschema.schema.MySchema`.
- YAML schemas are identified by the entry id of the schema file together with
the name of the section defined in the YAML schema. For example
if you have uploaded a schema YAML file containing a section definition called
`MySchema`, and it has been assigned an `entry_id`, the schema name will be
`entry_id:<entry_id>.MySchema`.

The quantities from schemas may be included or excluded as filter by using the
[`filters`](#filters) field in the app config. This option supports a
wildcard/glob syntax for including/excluding certain filters. For example, to
include all filters from the Python schema defined in the class
`myschema.schema.MySchema`, you could use:

```yaml
myapp:
  filters:
    include:
      - '*#myschema.schema.MySchema'
```

The same thing for a YAML schema could be achieved with:

```yaml
myapp:
  filters:
    include:
      - '*#entry_id:<entry_id>.MySchema'
```

Once quantities from a schema are included in an app as filters, they can be
targeted in the rest of the app. The app configuration often refers to specific
filters to configure parts of the user interface. For example,
one could configure the results table to show a new column using one of the
schema quantities with:

```yaml
myapp:
  columns:
    include:
      - 'data.mysection.myquantity#myschema.schema.MySchema'
      - 'entry_id'
    options:
      data.mysection.myquantity#myschema.schema.MySchema:
      ...
```

The syntax for targeting quantities depends on the resource:

- For python schemas, you need to provide the path and the python schema name separated
by a hashtag (#), for example `data.mysection.myquantity#myschema.schema.MySchema`.
- For YAML schemas, you need to provide the path and the YAML schema name separated
by a hashtag (#), for example `data.mysection.myquantity#entry_id:<entry_id>.MySchema`.
- Quantities that are common for all NOMAD entries can be targeted by using only
the path without the need for specifying a schema, e.g. `results.material.symmetry.space_group`.

## App configuration reference

{{ pydantic_model('nomad.config.models.config.App')}}
