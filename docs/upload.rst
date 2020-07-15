======================================
Uploading Data to the NOMAD Repository
======================================

To contribute your data to the repository, please, login to our `upload page <../uploads>`_ (you need to register first, if you do not have a NOMAD account yet).

*A note for returning NOMAD users!* We revised the upload process with browser based upload
alongside new shell commands. The new Upload page allows you to monitor upload processing
and verify processing results before publishing your data to the Repository.

The `upload page <../uploads>`_ acts as a staging area for your data. It allows you to
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

**A note for VASP users** on the handling of **POTCAR** files: NOMAD takes care of it; you don't
need to worry about it. We understand that according to your VASP license, POTCAR files are
not supposed to be visible to the public. Thus, in agreement with Georg Kresse, NOMAD will
extract the most important information of POTCAR files and store it in the files named
``POTCAR.stripped``. These files can be assessed and downloaded by anyone, while the original
POTCAR files are only available to the uploader and assigned co-authors.
This is done automatically; you don't need to do anything.

Once published, data cannot be erased. Linking a corrected version to a corresponding older one ("erratum") will be possible soon.
Files from an improved calculation, even for the same material, will be handled as a new entry.

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