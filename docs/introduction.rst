Introduction
============

This documentation is for the nomad software. This software can be used to
store, processing, and manage (computational)
material science data. It comprises a set of services and libraries that is
refered to as the *nomad infrastructure*.
The original nomad software was was developed as part of the
NOMAD-coe project (refered to as NOMAD-coe). This software nomad@FAIR
or just nomad is part of the FAIRDE.eu innitiative.

There are different use-modes for the nomad software, but the most common use is
to run the nomad infrastructure on a cloud and provide clients access to
web-based GUIs and REST APIs. This nomad infrastructure logically comprises the
*nomad repository* for uploading, searching, and downloading raw calculation output
from most relevant computionational material science codes. A second part of nomad
is the archive. It allows to access all uploaded data in a common data format
via a data schema called *meta-info* that includes common and code specific
specifications of structured data. Further services are available from
(nomad-coe.eu)[http://nomad-coe.eu], e.g. the  *nomad encyclopedia*, *analytics toolkit*,
and *advanced graphics*.

Architecture
------------

The following depicts the *nomad@FAIR* architecture with respect to software compenents
in terms of python modules, gui components, and 3rd party services (e.g. databases,
search enginines, etc.). It comprises a revised version of the repository and archive.

.. figure:: components.png
   :alt: nomad components

   The main modules of nomad

Nomad uses a series of 3rd party technologies that already solve most of nomads
processing, storage, availability, and scaling goals:

celery
^^^^^^
http://celeryproject.org (incl. rabbitmq) is a popular combination for realizing
long running tasks in internet applications. We use it to drive the processing of uploaded files.
It allows us to transparently distribute processing load.

elastic search
^^^^^^^^^^^^^^
Elastic search is used to store repository data (not the raw files).
Elastic search allows for flexible scalable search and analytics.

mongodb
^^^^^^^
Mongo is used to store and track the state of the processing of uploaded files and therein c
ontained calculations.

elastic stack
^^^^^^^^^^^^^
The *elastic stack* (previously *ELK* stack) is a central logging, metrics, and monitoring
solution that collects data within the cluster and provides a flexible analytics frontend
for said data.

Data model
----------

.. figure:: data.png
   :alt: nomad's data model

   The main data classes in nomad

See :py:mod:`nomad.processing`, :py:mod:`nomad.users`, and :py:mod:`nomad.repo`
for further information.

Processing
----------

.. figure:: proc.png
   :alt: nomad's processing workflow

   The workflow of nomad's processing tasks

See :py:mod:`nomad.processing` for further information.

Design principles
-----------------

- simple first, complicated only when necessary
- adopting generic established 3rd party solutions before implementing specific solutions
- only uni directional dependencies between components/modules, no circles
- only one language: Python (except, GUI of course)

General concepts
----------------

terms
^^^^^

There are is some terminology consistently used in this documentastion and the source
code:
- upload: A logical unit that comprises one (.zip) file uploaded by a user.
- calculation: A computation in the sense that is was created by an individual run of a CMS code.
- raw file: User uploaded files (e.g. part of the uploaded .zip), usually code input or output.
- upload file/uploaded file: The actual (.zip) file a user uploaded
- mainfile: The mainfile output file of a CMS code run.
- aux file: Additional files the user uploaded within an upload.
- repo entry: Some quantities of a calculation that are used to represent that calculation in the repository.
- archive data: The normalized data of one calculation in nomad's meta-info-based format.

ids and hashes
^^^^^^^^^^^^^^

Throughout nomad, we use different ids and hashes to refer to entities. If something
is called *id*, it is usually a random uuid and has no semantic conection to the entity
it identifies. If something is calles a *hash* than it is a hash build based on the
entitiy it identifies. This means either the whole thing or just some properties of
said entities.

The most common hashes are the *upload_hash* and *calc_hash*. The upload hash is
a hash over an uploaded file, as each upload usually refers to an indiviudal user upload
(usually a .zip file). The calc_hash is a hash over the mainfile path within an upload.
The combination of upload_hash and calc_hash is used to identify calculations. They
allow us to id calculations independently of any random ids that are created during
processing. To create hashes we use :func:`nomad.utils.hash`.