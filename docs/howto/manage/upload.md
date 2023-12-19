# How to upload and publish data for supported formats

This guide describes how to upload data in NOMAD [supported file formats](../../reference/parsers.md). You find a
list of supported formats on top of each upload page, see below.

## Preparing files

You can upload files one by one, but you can also provider larger `.zip` or `.tar.gz`
archive files, if this is easier to you. Also the file upload via frp or command line with
curl or with wget generates an archive files. The specific layout of these files is up to you.
NOMAD will simply extract them and consider the whole directory structure within.

## Create an upload and add files

Open [NOMAD](https://nomad-lab.eu/prod/v1){:target="_blank"} and log in; if you don't have a NOMAD account, please create one.

Go to `PUBLISH` / `Uploads`. Here you can create an upload with the `CREATE A NEW UPLOAD`
button. This will bring you to the upload page.

Before you start, make sure that the size of your data does not exceed the [upload limits](#upload-limits). If it does, please contact us.


You can drop your files on (or click) the `CLICK OR DROP FILES` button. On top you will
see a list of supported file formats and details on the files to upload.
You can also go to the `FILES` tab. Here you can create directories and drop files into directories.

## Processing files

NOMAD interprets your files. It checks each file and recognizes the main output file of the
supported codes. NOMAD creates an entry for this **mainfile** that represents the respective
data of this code run, experiment, etc.

While you can browse all files of an upload from its **upload page**, NOMAD only
allows to search for such recognized **mainfiles**. As long as your upload does not contain any
files that are recognized by NOMAD, you cannot publish the data.

However, all files that are associated to a recognized *mainfile* by being in the
same directory are displayed as **auxiliary** files next to the entry represented
by the **mainfile**.



!!! note
    **A note for VASP users**.
    On the handling of **POTCAR** files: NOMAD takes care of it; you don't
    need to worry about it. We understand that POTCAR files are not supposed to be visible to
    the public according to your VASP license. Thus, in agreement with Georg Kresse, NOMAD extracts
    the most important information of POTCAR files and stores it in the files named
    `POTCAR.stripped`. These files can be accessed and downloaded by anyone, while the original
    POTCAR files are automatically removed.


## Add user metadata

NOMAD automatically extracts as much information as possible from your files but you
can still specify additional metadata. This is what we call *user metadata*. This includes
you and your co-authors (use the *edit members* function on the *upload page*) as well
as comments, additional web-references, and datasets (use the *edit metadata* function on
the *upload page*).

User metadata can also be provided in an uploaded file. This can be a `.json` or
`.yaml` file. It has to be named `nomad.json` or `nomad.yaml`. Here is a JSON example:

```json
{
    "comment": "Data from a cool research project",
    "references": ["http://archivex.org/mypaper"],
    "coauthors": [
        "<email-or-username>",
        "<email-or-username>"
    ],
    "datasets": [
        "<dataset-name>"
    ],
    "entries": {
        "path/to/entry_dir/vasp.xml": {
            "comment": "An entry specific comment."
        }
    }
}
```

This file is only applied during the initial processing of an entry. So make sure you either
upload it first or with everything else as part of an archive file.

## Publish and get a DOI

After clicking the `PUBLISH` button, the uploaded files will become immutable, but you can still
edit the metadata.

As part of the *edit metadata* functionality, you can create and assign *datasets*.
Go to `PUBLISH` / `Datasets` in the menu to see all your datasets. Here you can assign
a DOI to created *datasets*. For a *dataset* with DOI, you can only add more entries, but
not remove entries.


## Upload limits

- One upload cannot exceed **32 GB** in size.
- Only **10 non published uploads** are allowed per user.
- Only uploads with at least one recognized entry can be published. See also [supported codes/formats](../../reference/parsers.md) below.


## Strategies for large amounts of data

Before attempting to upload large amounts of data, run some experiments with a representative
and small subset of your data. Use this to simulate a larger upload that you can review and edit
in the normal way. You do not have to publish this test upload; simply delete it before publish,
once you are satisfied with the results.

Ask for assistance and [Contact us](https://nomad-lab.eu/about/support){:target="_blank"} in advance. This will
allow us to react to your specific situation and eventually prepare additional measures.
Allow enough time before you need your data to be published. Adding multiple hundreds of
GBs to NOMAD isn't a trivial feat and will take some time and effort on all sides.

The [upload limits](#limits) above are necessary to keep NOMAD data manageable and we cannot easily
grant exceptions to these rules. This means you have to split your data into 32 GB uploads.
Uploading these files, observing the processing, and publishing the data can be automatized through NOMAD APIs.

When splitting your data, it is important to not split subdirectories that contain files of the same single entry. NOMAD can only bundle those related files to an entry if
they are part of the same upload (and directory). Therefore, there is no single recipe to
follow, and a script to split your data depends heavily on how your data is organized.

If you provide data for a potentially large amount of entries, it might be advisable
to provide *user metadata* via file. See [user metadata](#user-metadata) above for details.

To further automate, you can also upload and directly publish data. After performing some
smaller test uploads, you should consider skipping our staging and publish the upload
right away. This can save you some time and additional API calls. The upload endpoint
has a parameter `publish_directly`. You can modify the upload command you get on the upload page as follows:

```
curl "http://nomad-lab.eu/prod/v1/uploads/?token=<your-token>&publish_directly=true" -T <local_file>
```

HTTP makes it easy for you to upload files via browser and curl, but it is not an
ideal protocol for the stable transfer of large and many files. Alternatively, we can organize
a separate manual file transfer to our servers. We will put your prepared upload
files (.zip or .tag.gz) on a predefined path on the NOMAD servers. NOMAD allows to *"upload"*
files directly from its servers via an additional `local_path` parameter:

```
curl -X PUT "http://nomad-lab.eu/prod/v1/api/uploads/?token=<your-token>&local_path=<path-to-upload-file>"
```
