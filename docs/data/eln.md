This guide describes how to manually create entries and enter information
via ELNs (electronic lab notebooks). NOMAD ELNs allow you to acquire
consistently structured data from users to augment uploaded files.

!!! attention

    This part of the documentation is still work in progress.

## Create a basic ELN entry

Go to `PUBLISH` / `Uploads`. Here you can create an upload with the `CREATE A NEW UPLOAD`
button. This will bring you to the upload page.

Click the `CREATE ENTRY` button. This will bring-up a dialog to choose an ELN schema.
All ELNs (as any entry in NOMAD) needs to follow a schema. You can choose from uploaded
custom schemas or NOMAD build-in schemas. You can choose the `Basic ELN` to create a
simple ELN entry.

The name of your ELN entry, will be the filename for your ELN without the `.archive.json`
ending that will be added automatically. You can always find and download your ELNs
on the `FILES` tab.

The `Basic ELN` offers you simple fields for a *name*, *tags*, a *date/time*, and a rich text
editor to enter your notes.

## Add your own ELN schema

To make NOMAD ELNs more useful, you can define your own schema to create you own data
fields, create more sub-sections, reference other entries, and much more.

You should have a look at our ELN example upload. Go to `PUBLISH` / `Uploads` and
click the `ADD EXAMPLE UPLOADS` button. The `Electronic Lab Notebook` example, will
contain a schema and entries that instantiate different parts of the schema.
The *ELN example sample (`sample.archive.json`) demonstrates what you can do.

Follow the [How to write a schema](../schemas/basics.md) and [How to define ELN](../schemas/elns.md)
guides to create you own customized of ELNs.
