# Using the web interface

For our upcoming tutorial on metadata management with NOMAD, we will produce a serious
of videos that demonstrate the web interface of NOMAD v1. The videos will show how to

- upload and publish data
- search and download data
- use the archive

This will come in February 2022.


## Uploading and publishing data


### Preparing files

You can upload files one by one, but you can also provider larger `.zip` or `.tar.gz`
archive files, if this is easier to you. Also the file upload frp, the command line with
curl or wget, is based an archive files. The specific layout of these files is up to you.
NOMAD will simply extract them and consider the whole directory structure within.


### Supported codes

NOMAD is interpreting your files. It will check each file and recognize if it is the
main output file of one of the supported codes. NOMAD will create a entry for this *mainfile*
that represents the respective data of this code run, experiment, etc.

While you can browse all files of an upload from its *upload page*, NOMAD only
allows to search for such recognized entries. As long as your upload does not contain any
files that NOMAD recognizes, you cannot be publish the upload.

However, all files that are associated to a recognized *mainfile* by residing in the
same directory, will be presented as *auxiliary* files along side the entry represented
by the *mainfile*.

**A note for VASP users**

On the handling of **POTCAR** files: NOMAD takes care of it; you don't
need to worry about it. We understand that according to your VASP license, POTCAR files are
not supposed to be visible to the public. Thus, in agreement with Georg Kresse, NOMAD will
extract the most important information of POTCAR files and store it in the files named
`POTCAR.stripped`. These files can be assessed and downloaded by anyone, while the original
POTCAR files are automatically removed.


### User metadata

NOMAD will automatically extract as much information as possible from your files. But you
can still provide additional metadata. This is what we call *user metadata*. This includes
you and your co-authors (use the *edit members* function on the *upload page*) as well
as comments, additional web-references, and datasets (use the *edit metadata* function on
the *upload page*).

User metadata can also be provided in a file that is uploaded. This can be a `.json` or
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
        "path/to/calcs/vasp.xml": {
            "commit": "An entry specific comment."
        }
    }
}
```

This file is only applied on initial processing of an entry. So make sure you either
upload it first, or upload everything as part of one archive file.


### Limits

- One upload cannot exceed **32 GB** in size.
- Only **10 non published uploads** are allowed per user.
- Only uploads at least one recognized entry can be published. See also [supported codes](#supported-codes) below.


### Strategies for large amounts of data

Before you attempt to upload large amounts of data, do some experiments with a representative
and small subset of your data. Use this to simulate a larger upload,
checking and editing it the normal way. You do not have to publish this test upload;
simply delete it before publish, once you are satisfied with the results.

Ask for assistance. [Contact us](https://nomad-lab.eu/about/contact) in advance. This will
allow us to react to your specific situation and eventually prepare additional measures.
Keep enough time before you need your data to be published. Adding multiple hundreds of
GBs to NOMAD isn't a trivial feat and will take some time and effort from all sides.

The [upload limits](#limits) above are necessary to keep NOMAD data manageable and we cannot easily
grant exceptions to these rules. This means you have to split your data into 32 GB uploads. Uploading these files, observing the processing, and publishing the data can be automatized through NOMAD APIs.

When splitting your data, it is important to not split sub-directories if they represent
all files of a single entry. NOMAD can only bundle those related files to an entry if
they are part of the same upload (and directory). Therefore, there is no single recipe to
follow and a script to split your data will depend on how your data is organized.

If you provide data for a potentially large amount of entries, it might be advisable
to provide *user metadata* via file. See [user metadata](#user-metadata) above for details.

To further automate, you can also upload and directly publish data. After performing some
smaller test upload, you should consider to skip our staging and publish the upload
right away. This can save you some time and additional API calls. The upload endpoint
has a parameter `publish_directly`. You can modify the upload command
that you get from the upload page like this:

```
curl "http://nomad-lab.eu/prod/rae/api/uploads/?token=<your-token>&publish_directly=true" -T <local_file>
```

HTTP makes it easy for you to upload files via browser and curl, but it is not an
ideal protocol for the stable transfer of large and many files. Alternatively, we can organize
a separate manual file transfer to our servers. We will put your prepared upload
files (.zip or .tag.gz) on a predefined path on the NOMAD servers. NOMAD allows to *"upload"*
files directly from its servers via an additional `local_path` parameter:

```
curl -X PUT "http://nomad-lab.eu/prod/rae/api/uploads/?token=<your-token>&local_path=<path-to-upload-file>"
```