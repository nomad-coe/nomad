## Why a federated data infrastructure?

There are several benefits for using multiple NOMAD installations:

- **Sovereignty**: Data is stored and managed locally. This is important for data privacy and security.
- **Resource**: Local resources can be used to manage and analyse data. This is important for large data sets.
- **Customization**: Local installations can be customized (plugins). This is important to support local workflows and special file formats.

## What is a NOMAD Oasis?

The software that runs NOMAD is Open-Source and can be used independently of the NOMAD
_central installation_ at [http://nomad-lab.eu](http://nomad-lab.eu){:target="\_blank"}.
We call any NOMAD installation that is not the _central_ one a NOMAD Oasis.

## Use cases

There are several use-cases how the NOMAD software could be used. Of course other
uses and hybrids are imaginable:

- Academia: Use the Oasis for local management of unpublished research data
- Industry: Use of Oasis to manage private data and full internal use of published data in compliance with strict privacy policies
- Mirror: Use the Oasis as a mirror that hosts a copy of all published NOMAD data
- FAIRmat: Use Oasis to form a network of repositories to build a federated data infrastructure
  for materials science.
  This is what we do in the [FAIRmat project](https://www.fair-di.eu/fairmat/consortium){:target="\_blank"}.

<figure markdown>
  ![oasis use-cases](images/oasis-use-cases.png){ width=700 }
  <figcaption>NOMAD Oasis use-cases</figcaption>
</figure>

## Data transfer

!!! warning "Attention"

    The data transfer implementation is still under development and not
    yet relyable for all use-cases. See "Current Limitations and Plans" below.

NOMAD data is inherently file based. You upload raw-files (or in case of ELNs edit raw files
in the web interface) and the NOMAD processes these files to extract data which
again is stored in files. Therefore, data transfer is basically "just" copying
files from one installation to another.

We implemented a _bundle_ format for published data. A bundle is a zip file that contains
all raw files, archive files of an _upload_, and a bit of metadata, e.g. with
information on upload name, dates, authors, datasets, etc. This bundle format
is our tranfer file format. Bundles are exported on one end and imported on
the other.

Upon import, the reciving NOMAD replicates the state of the upload on the
sending NOMAD. This means that the upload is published and the data is indexed.
Users can treat the upload as if they had uploaded it on the receiving
NOMAD themselves.

There are two principle ways to transfer data:

- **Automated**: The Oasis user interface offers an option to transfer data to the configured
  central NOMAD installation. This uses an API endpoint on the Oasis that creates a bundle and pushes it
  to the central NOMAD installation via the API of the central NOMAD.

- **Manual**: You can use the NOMAD CLI to export uploads as bundle files and
  upload them manually to another NOMAD installation.

- **Re-upload**: Of course you could also skip the bundle format entrily and just
  re-upload the raw files to another NOMAD installation.

## Current Limitations and Plans

The basic bundle format, export, and import, and a transfer button in the UI are implemented.
However, there are still limitations and plans for future development:

- **Trust issues**: Some of the bundle metadata is normally created by
  NOMAD. Timestamps like upload times are the biggest concerns here.
  This means that the receiving NOMAD needs to trust that the timestamps in a bundle
  are not forged. Therefore, it needs to verify that the bundle
  was actually created by a trusted Oasis. Currently, Oasis have to be
  configured to use a trusted user (usually the Oasis admin) to submit bundles
  to the central NOMAD.
  A version were timestamps are ignored and the upload is dated to the time
  it was received, would also be possible implementation, but is not implemented
  yet.

- **Plugins and versions**: Two installations (e.g. Oasis and central NOMAD)
  might not have the same plugins installed and not have the same version
  of all schemas available. Received data might not be fully processable.
  Therefore, it is planned, but not yet implemented, to make all required schemas
  part of the bundle.

- **Metadata transfers**: Currently, only full uploads with all raw files
  and all archive files can be transferred.
  It is planned to allow transfers for only the "metadata". Here metadata
  refers to the archive files (or even just certain result sections from
  said archive files). This would still allow to index and find entries on the
  central NOMAD without transfering potentially huge raw files. This would allow
  the central NOMAD to act as a portal for all data in the federation. This is not
  yet implemented.

## References to the implementation

There are no detailed how-tos yet. Here are some pointers to reference documentation that might be relevant:

- [The API endpoint to download an upload in bundle format.](https://nomad-lab.eu/prod/v1/api/v1/extensions/docs#/uploads%2Fbundle/get_upload_bundle_uploads__upload_id__bundle_get)
- [The API endpoint to upload a bundle.](https://nomad-lab.eu/prod/v1/api/v1/extensions/docs#/uploads%2Fbundle/post_upload_bundle_uploads_bundle_post)
- [The API endpoint to publish an upload, including publishing to the configured central NOMAD.](https://nomad-lab.eu/prod/v1/api/v1/extensions/docs#/uploads%2Faction/post_upload_action_publish_uploads__upload_id__action_publish_post)
- [The CLI commands to export and import an upload as a bundle.](http://127.0.0.1:8000/reference/cli.html#nomad-admin-uploads-export)
