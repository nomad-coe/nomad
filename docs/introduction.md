# Introduction

**NOvel Materials Discorvery (NOMAD)** comprises storage, processing, management, discovery, and
analytics of computational material science data from over 40 community *codes*.
The original NOMAD software, developed by the
[NOMAD-coe](http://nomad-coe.eu) project, is used to host over 50 million total energy
calculations in a single central infrastructure instances that offers a variety
of services (*repository*, *archive*, *encyclopedia*, *analytics*, *visualization*).

.. figure:: nomad.png
   :alt: nomad's overall structure

This is the documentation of **nomad@FAIRDI**, the Open-Source continuation of the
original NOMAD-coe software that reconciles the original code base,
integrate it's services, allows 3rd parties to run individual and federated instance of
the nomad infrastructure, provides nomad to other material science domains, and applies
the FAIRDI principles as proliferated by the [FAIRDI Data Infrastructure e.V.](http://fairdi.eu).
A central and publically available instance of the nomad software is run at the
[MPCDF](https://www.mpcdf.mpg.de/) in Garching, Germany.

The nomad software runs SAAS on a server and is used via web-based GUI and ReSTful
API. Originally developed and hosted as individual services, **nomad@FAIRDI**
provides all services behind one GUI and API into a single coherent, integrated, and
modular software project.

This documentation is only about the nomad *software*; it is about architecture,
how to contribute, code reference, engineering and operation of nomad. It is not a
nomad user manual.


## Architecture

The following depicts the *nomad@FAIRDI* architecture with respect to software components
in terms of python modules, gui components, and 3rd party services (e.g. databases,
search engines, etc.). It comprises a revised version of the repository and archive.

.. figure:: components.png
   :alt: nomad components

Besides various scientific computing, machine learning, and computational material
science libraries (e.g. numpy, skikitlearn, tensorflow, ase, spglib, matid, and many more),
Nomad uses a set of freely available or Open Source technologies that already solve most
of its processing, storage, availability, and scaling goals. The following is a non
comprehensive overview of used languages, libraries, frameworks, and services.


### Python 3

The *backend* of nomad is written in Python. This includes all parsers, normalizers,
and other data processing. We only use Python 3 and there is no compatibility with
Python 2. Code is formatted close to [pep8](https://www.python.org/dev/peps/pep-0008/),
critical parts use [pep484](https://www.python.org/dev/peps/pep-0484/) type-hints.
[Pycodestyle](https://pypi.org/project/pycodestyle/),
[pylint](https://www.pylint.org/), and
[mypy](http://mypy-lang.org/) (static type checker) are used to ensure quality.
Tests are written with [pytest](https://docs.pytest.org/en/latest/contents.html).
Logging is done with [structlog](https://www.structlog.org/en/stable/) and *logstash* (see
Elasticstack below). Documentation is driven by [Sphinx](http://www.sphinx-doc.org/en/master/).


### celery

[Celery](http://celeryproject.org) (+ [rabbitmq](https://www.rabbitmq.com/))
is a popular combination for realizing long running tasks in internet applications.
We use it to drive the processing of uploaded files.
It allows us to transparently distribute processing load while keeping processing state
available to inform the user.


### elastic search

[Elasticsearch](https://www.elastic.co/webinars/getting-started-elasticsearch)
is used to store repository data (not the raw files).
Elasticsearch allows for flexible scalable search and analytics.


### mongodb

[Mongodb](https://docs.mongodb.com/) is used to store and track the state of the
processing of uploaded files and therein contained calculations. We use
[mongoengine](http://docs.mongoengine.org/) to program with mongodb.


### PostgreSQL

A relational database is used to store all user provided metadata: users, datasets
(curated sets of uploaded data), references, comments, DOIs, coauthors, etc.
Furthermore, it is still used to store some of the calculation metadata derived
via parsing. *This will most likely move out of Postgres in the future.* We
use [SQLAlchemy](https://docs.sqlalchemy.org/en/latest/) as on ORM framework.


### flask, et al.

The ReSTful API is build with the [flask](http://flask.pocoo.org/docs/1.0/)
framework and its [ReST+](https://flask-restplus.readthedocs.io/en/stable/) extension. This
allows us to automatically derive a [swagger](https://swagger.io/) description of the nomad API,
which in  turn allows us to generate programming language specific client libraries, e.g. we
use [bravado](https://github.com/Yelp/bravado) for Python and
[swagger-js](https://github.com/swagger-api/swagger-js) for Javascript.
Fruthermore, you can browse and use the API via [swagger-ui](https://swagger.io/tools/swagger-ui/).


### Elasticstack

The [elastic stack](https://www.elastic.co/guide/index.html)
(previously *ELK* stack) is a central logging, metrics, and monitoring
solution that collects data within the cluster and provides a flexible analytics frontend
for said data.


### Javascript, React, Material-UI

The frontend (GUI) of **nomad@FAIRDI** build on top of the
[React](https://reactjs.org/docs/getting-started.html) component framework.
This allows us to build the GUI as a set of re-usable components to
achieve a coherent representations for all aspects of nomad, while keeping development
efforts manageable. React uses [JSX](https://reactjs.org/docs/introducing-jsx.html)
(a ES6 variety) that allows to mix HTML with Javascript code.
The component library [Material-UI](https://material-ui.com/)
(based on Google's popular material design framework) provides a consistent look-and-feel.


### docker

To run a **nomad@FAIRDI** instance, many services have to be orchestrated:
the nomad app, nomad worker, mongodb, Elasticsearch, PostgreSQL, RabbitMQ,
Elasticstack (logging), the nomad GUI, and a reverse proxy to keep everything together.
Further services might be needed (e.g. JypiterHUB), when nomad grows.
The container platform [Docker](https://docs.docker.com/) allows us to provide all services
as pre-build images that can be run flexibly on all types of platforms, networks,
and storage solutions. [Docker-compose](https://docs.docker.com/compose/) allows us to
provide configuration to run the whole nomad stack on a single server node.


### kubernetes + helm

To run and scale nomad on a cluster, you can use [kubernetes](https://kubernetes.io/docs/home/)
to orchestrated the  necessary containers. We provide a [helm](https://docs.helm.sh/)
chart with all necessary service and deployment descriptors that allow you to setup and
update nomad with few commands.


### GitLab

Nomad as a software project is managed via [GitLab](https://docs.gitlab.com/).
The **nomad@FAIRDI** project is hosted [here](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR).
GitLab is used to manage versions, different branches of development, tasks and issues,
as a [registry for Docker images](https://docs.gitlab.com/ee/user/project/container_registry.html),
and [CI/CD platform](https://docs.gitlab.com/ee/ci/).


## Data model

.. figure:: data.png
   :alt: nomad's data model

The entities that comprise the nomad data model are *users*, *datasets*, *uploads*,
*calculations* (calc), and *materials*. *Users* upload multiple related *calculations*
in one *upload*. *Users* can curate *calculations* into *datasets*. *Caclulations*
belong to one *material* based on the simulated system.

### Users
- The user `email` is used as a primary key to uniquely identify users
  (even among different nomad installations)

### Uploads
- An upload contains related calculations in the form of raw code input and output files
- Uploader are encouraged to upload all relevant files
- The directory structure of an upload might be used to relate calculations to each other
- Uploads have a unique randomly choosen `upload_id` (UUID)
- The `uploader` is the user that provided the upload. There is always one immutable `uploader`
- Currently, uploads can be provided as `.zip` or `.tar.gz` files.

### Entries (Calculations, Code runs)
- There are confusing names. Internally, in the nomad source code, the term `calc` is used.
An entry represents a single set of input/output used and produces by an individual run of a
DFT code. If nomad is applied to other domains, i.e. experimental material science, entries might represent
experiments or other entities.
- An entry (calculation) has a unique `calc_id` that is based on the upload's id and the `mainfile`
- The `mainfile` is a upload relative path to the main output file.
- Each calculation, when published, gets a unique `pid`. Pids are ascending intergers. For
each `pid` a shorter `handle` is created. Handles can be registered with a handle system,
e.g. the central nomad installation at MPCDF is registered at a MPCDF/GWDW handle system.
- The `calc_hash` is computed from the main and other parsed raw files.
- Entry data comprises *user metadata* (comments, references, datasets, coauthors),
*calculation metadata* (code, version, system and symmetry, used DFT method, etc.),
the *archive data* (a hierarchy of all parsed quantities), and the uploaded *raw files*.

### Datasets
- Datasets are user curated sets of calculations.
- Users can assign names and nomad can register a DOI for a dataset.
- A calculation can be put into multiple datasets.

### Materials
- Materials aggregate calculations based on common system properties
(e.g. system type, atoms, lattice, space group, etc.).

### Data
We distinguish various forms of calculation data:
- raw data: The raw files provided by nomad users
- (repository) metadata: All data necessary to search and inspect nomad entries.
- archive data: The data extracted from raw files by nomad parsers and normalizers.
This data is represented in the *meta-info* format.
- materials data: Aggregated information about calculations that simulated the *same* material.

.. figure:: datamodel_dataflow.png
   :alt: nomad's data flow

### Metadata
Metadata refers to those pieces of data, those quantities/attributes that we use
to represent, identify, and index uploads and calculations in the API, search, GUI, etc.
There are three catergories of metadata:
- entry metadata: attributes that are necessary to uniquely identify entities (see also :ref:`id-reference-label`), that describe the upload, processing, etc. This data is derived by the nomad infrastructure.
- user metadata: attributes provided by the user, e.g. comments, references, coauthors, datasets, etc. This data is provided by the user.
- domain metadata: metadata parsed from raw files that describe calculations on a high level, e.g. code name, basis set, system type, etc. This data is derived from the uploaded data.

Those sets of metadata along with the actual raw and archive data are often transformed,
passed, stored, etc. by the various nomad modules.

.. figure:: datamodel_metadataflow.png
   :alt: nomad's metadata flow

### Implementation
The different entities have often multiple implementations for different storage systems.
For example, aspects of calculations are stored in files (raw files, calc metadata, archive data),
Postgres (user metadata), Elasticsearch (metadata), and mongodb (processing state).
Different transformation between different implementations exist. See
:py:mod:`nomad.datamodel` for further information.


## Processing

.. figure:: proc.png
   :alt: nomad's processing workflow

See :py:mod:`nomad.processing` for further information.
