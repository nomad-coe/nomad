Introduction
============

**NOvel Materials Discorvery (NOMAD)** comprises storage, processing, management, discovery, and
analytics of computational material science data from over 40 community *codes*.
The original NOMAD software, developed by the
[NOMAD-coe](http://nomad-coe.eu) project, is used to host over 50 million total energy
calculations in a single central infrastructure instances that offers a variety
of services (repository, archive, encyclopedia, analytics, visualization).

This is the documentation of **nomad@FAIR**, the Open-Source continuation of the
original NOMAD-coe software that reconciles the original code base,
integrate it's services,
allows 3rd parties to run individual and federated instance of the nomad infrastructure,
provides nomad to other material science domains, and applies the FAIR principles
as proliferated by the (FAIR Data Infrastructure e.V.)[http://fairdi.eu].

There are different use-modes for the nomad software, but the most common use is
to run the nomad infrastructure on a cloud and provide clients access to
web-based GUIs and REST APIs. This nomad infrastructure logically comprises the
*nomad repository* for uploading, searching, and downloading raw calculation input and output
from all relevant computational material science codes. A second part of nomad
is the *archive*. It provides all uploaded data in a common data format
called *meta-info* and includes common and code specific
schemas for structured data. Further services are available from
(nomad-coe.eu)[http://nomad-coe.eu], e.g. the  *nomad encyclopedia*, *analytics toolkit*,
and *advanced graphics*.

Architecture
------------

The following depicts the *nomad@FAIR* architecture with respect to software components
in terms of python modules, gui components, and 3rd party services (e.g. databases,
search engines, etc.). It comprises a revised version of the repository and archive.

.. figure:: components.png
   :alt: nomad components

   The main modules of nomad

Nomad uses a series of 3rd party technologies that already solve most of nomads
processing, storage, availability, and scaling goals:

celery
^^^^^^
http://celeryproject.org (incl. rabbitmq) is a popular combination for realizing
long running tasks in internet applications. We use it to drive the processing of uploaded files.
It allows us to transparently distribute processing load while keeping processing state
available to inform the user.

elastic search
^^^^^^^^^^^^^^
Elastic search is used to store repository data (not the raw files).
Elastic search allows for flexible scalable search and analytics.

mongodb
^^^^^^^
Mongo is used to store and track the state of the processing of uploaded files and therein
contained calculations.

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
