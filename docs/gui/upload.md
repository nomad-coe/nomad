# Uploading and publishing data

To contribute your data to the repository, please, login to our [upload page](https://nomad-lab.eu/prod/rae/gui/uploads)
(you need to register first, if you do not have a NOMAD account yet).

*A note for returning NOMAD users!* We revised the upload process with browser based upload
alongside new shell commands. The new Upload page allows you to monitor upload processing
and verify processing results before publishing your data to the Repository.

The [upload page](https://nomad-lab.eu/prod/rae/gui/uploads) acts as a staging area for your data. It allows you to
upload data, to supervise the processing of your data, and to examine all metadata that
NOMAD extracts from your uploads. The data on the upload page will be private and can be
deleted again. If you are satisfied with our processing, you can publish the data.
Only then, data will become publicly available and cannot be deleted anymore.
You will always be able to access, share, and download your data. You may curate your data
and create datasets to give them a hierarchical structure. These functions are available
from the Your data page by selecting and editing data.

You should upload many files at the same time by creating .zip or .tar files of your folder structures.
Ideally, input and output files are accompanied by relevant auxiliary files. NOMAD will
consider everything within a single directory as related.

Once published, data cannot be erased. Linking a corrected version to a corresponding older
one ("erratum") will be possible soon. Files from an improved calculation, even for the
same material, will be handled as a new entry.

You can publish data as being open access or restricted for up to three years (with embargo).
For the latter you may choose with whom you want to share your data. We strongly support the
idea of open access and thus suggest to impose as few restrictions as possible from the very
beginning. In case of open access data, all uploaded files are downloadable by any user.
Additional information, e.g. pointing to publications or how your data should be cited,
can be provided after the upload. Also DOIs can be requested. The restriction on data
can be lifted at any time. You cannot restrict data that was published as open access.

Unless published without an embargo, all your information will be private and only visible
to you (or NOMAD users you explicitly shared your data with). Viewing private data will
always require a login.

By uploading you confirm authorship of the uploaded calculations. Co-authors must be specified
after the upload process. This procedure is very much analogous to the submission of a
publication to a scientific journal.

Upload of data is free of charge.

## Limits

The following limitations apply to uploading:

- One upload cannot exceed 32 GB in size
- Only 10 non published uploads are allowed per user

## On the supported codes

NOMAD is interpreting your files. It will check each file and recognize if it is the
main output file of one of the supported codes. NOMAD will create a entry for this *mainfile*
that represents the respective data of this code run, experiment, etc. NOMAD only
shows that for such recognized entries. If you uploads do not contain any files that
NOMAD recognizes, you upload will be shown as empty and no data can be published.

However, all files that are associated to a recognized *mainfile* by residing in the
same directory, will be presented as *auxiliary* files along side the entry represented
by the *mainfile*.

### A note for VASP users

On the handling of **POTCAR** files: NOMAD takes care of it; you don't
need to worry about it. We understand that according to your VASP license, POTCAR files are
not supposed to be visible to the public. Thus, in agreement with Georg Kresse, NOMAD will
extract the most important information of POTCAR files and store it in the files named
`POTCAR.stripped`. These files can be assessed and downloaded by anyone, while the original
POTCAR files are removed when the upload is published. This is done automatically; you don't
need to do anything.

## Preparing an upload file

You can upload .zip and .tar.gz files to NOMAD. The directory structure within can
be arbitrary. Keep in mind that files in a single directory are all associated (see above).
Ideally you only keep the files of a single (or closely related) code runs, experiments, etc.
in one directory.

You should not place files in additional archives within the upload file. NOMAD will not
extract any zips in zips and similar entrapments.

## Uploading large amounts of data

This problem is many fold. In the remainder the following topics are discussed.

- NOMAD restrictions about upload size and number of unpublished simultaneous uploads
- Managing metadata (comments, references, co-authors, datasets) for a large number of entries
- Safely transferring the data to NOMAD

### General strategy

Before you attempt to upload large amounts of data, do some experiments with a representative
and small subset of your data. Use this to simulate a larger upload,
checking and editing it the normal way. You do not have to publish this test upload;
simply delete it before publish, once you are satisfied with the results.

Ask for assistance. [Contact us](https://nomad-lab.eu/about/contact) in advance. This will
allow us to react to your specific situation and eventually prepare additional measures.

Keep enough time before you need your data to be published. Adding multiple hundreds of
GBs to NOMAD isn't a trivial feat and will take some time and effort from all sides.

### Upload restrictions

The upload restrictions are necessary to keep NOMAD data in manageable chunks and we cannot
simply grant exceptions to these rules.

This means you have to split your data into 32 GB uploads. Uploading these files, observing
the processing, and publishing the data can be automatized through NOMAD APIs.

When splitting your data, it is important to not split sub-directories if they represent
all files of a single entry. NOMAD can only bundle those related files to an entry if
they are part of the same upload (and directory). Therefore, there is no single recipe to
follow and a script to split your data will depend on how your data is organized.

### Avoid additional operations on your data

Changing the metadata of a large amounts of entries can be expensive and will also mean
more work with our APIs. A simpler solution is to add the metadata directly to your uploads.
This way NOMAD can pick it up automatically, no further actions required.

Each NOMAD upload can contain a `nomad.json` file at the root. This file can contain
metadata that you want to apply to all your entries. Here is an example:

```
{
    "comment": "Data from a cool research project",
    "references": ['http://archivex.org/mypaper'],
    "coauthors": [
        '<co-author-ids>',
        '<co-author-ids>'
    ]
    "datasets": [
        '<dataset-id>'
    ],
    "entries": {
        "path/to/calcs/vasp.xml": {
            "commit": "An entry specific comment."
        }
    }
}
```

Another measure is to directly publish your data upon upload. After performing some
smaller test upload, you should consider to skip our staging and publish the upload
right away. This can save you some time and additional API calls. The upload endpoint
has a parameter `publish_directly`. You can modify the upload command
that you get from the upload page like this:

```
curl "http://nomad-lab.eu/prod/rae/api/uploads/?token=<your-token>&publish_directly=true" -T <local_file>
```

### Save transfer of files

HTTP makes it easy for you to upload files via browser and curl, but it is not an
ideal protocol for the stable transfer of large and many files. Alternatively, we can organize
a separate manual file transfer to our servers. We will put your prepared upload
files (.zip or .tag.gz) on a predefined path on the NOMAD servers. NOMAD allows to *"upload"*
files directly from its servers via an additional `local_path` parameter:

```
curl -X PUT "http://nomad-lab.eu/prod/rae/api/uploads/?token=<your-token>&local_path=<path-to-upload-file>"
```
