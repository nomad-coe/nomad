[![pipeline status](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/master/pipeline.svg)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/master)
[![coverage report](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/master/coverage.svg)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/master)

This project implements the new *nomad@FAIRDI* infrastructure. It is currently used
to enable users to upload data, process the data, maintain a version of the NOMAD
archive and meta-info, provide search, inspection, and download to all NOMAD raw and
archive data. As a long term strategy, this project will integrate, refactor, and re-write
more and more of the existing NOMAD CoE components.

The overall goal of *nomad@FAIRDI* is to provide common interfaces to the main services of NOMAD:
*Repository*, *Archive*, and *Encyclopedia*. These interfaces comprise a graphical web-based
UI that allows users to upload data, supervise data processing, inspect and download
metadata, raw-files, and archive data, provide visual tools to explore the data, and
to learn more about advanced use modes, like API and Analytics Toolkit. The second interface
is a unified REST API with various endpoints that represent the core NOMAD services.
This will allow users the automated use of NOMAD for managing their data, and using
data on NOMAD for analytics. A specific way of using the API is through the NOMAD
Analytics Toolkit, which is revamped as a
[separate project](https://gitlab.mpcdf.mpg.de/nomad-lab/analytics-jupyterhub).

Furthermore, this projects aims at establishing NOMAD as a distributed platform for
material science data sharing and management. This includes the on-site deployment of
NOMAD as a standalone service (*oasis*), the federated use of NOMAD through a
serious of full and partial *mirrors*, the integration of 3rd party material science
databases (i.e. [Aflow](http://www.aflow.org/), [OQMD](http://oqmd.org/),
[Materials Project](https://materialsproject.org/)), and support for open APIs and
standards like the [Optimade](http://www.optimade.org/) API.

## Getting started

Read the docs. The documentation is part of the source code. It covers aspects like
introduction, architecture, development setup/deployment, contributing, and API reference.

### Read the docs on the latest deployed version

You can access the running system and its documentation here:

[https://repository.nomad-coe.eu/uploads/api/docs](https://repository.nomad-coe.eu/uploads/api/docs/index.html)

### Generate the docs from the source

First, clone this repo and init its submodules:
```
git clone git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-FAIR.git
cd nomad-FAIR
git submodule init --depth 1
```

Second, create and source your own virtual python environment:
```
pip install virtualenv
virtualenv -p `which python3` .pyenv
source .pyenv/bin/activate
```

Third, install the development dependencies, including the documentation system
[sphinx](http://www.sphinx-doc.org/en/master/index.html):
```
pip install --upgrade pip
pip install --upgrade setuptools
pip install -r requirements.txt
```

Forth, generate the documentation:
```
cd docs
make html
```

Conintue with reading the documentation for further setup and contribution guidelines:
```
cd .build/html
python -m http.server 8888
```
Open [http://localhost:8888/html/setup.html](http://localhost:8888/html/setup.html) in
your browser.

## Change log
Omitted versions are plain bugfix releases with only minor changes and fixes.

### v0.6.0
- Keycloak based user management
- no dependencies with the NOMAD CeE Repository

### v0.5.2
- allows to download large files over longer time period
- streamlined deployment without API+GUI proxy
- minor bugfixes

### v0.5.1
- integrated parsers Dmol3, qbox, molcas, fleur, and onetep
- API endpoint for query based raw file download
- improvements to admin cli: e.g. clean staging files, reprocess uploads based on codes
- improved error handling in the GUI
- lots of parser bugfixes
- lots of minor bugfixes

### v0.5.0
The first production version of nomad@fairdi as the upload API and gui for NOMAD
- Production ready software and deployments (term agreements, better GUI docs)
- Raw file API with support to list directories. This replaces the `files` calculation
  metadata key. It was necessary due to arbitrary large lists of *auxfiles* in some
  calculations.
- Search interface that contains all features of the CoE Repository GUI.
- Refactored search API that allows to search for entries (paginated + scroll),
  metrics based on quantity aggregations (+ paginated entries), quantity aggregations
  with all values via `after` key (+ paginated entries).
- reprocessing of published results (e.g. after parser/normalizer improvements)
- mirror functionality
- refactored command line interface (CLI)
- potential GUI user tracking capabilities
- many minor bugfixes

### v0.4.7
- more migration scripts
- minor bugfixes

### v0.4.6
- admin commands to directly manipulate upload data
- additional migration scripts
- fixed system normalizer to understand indexed atom labels correctly
- many minor bugfixes

### v0.4.5
- improved uploads view with published uploads
- support for publishing to the existing nomad CoE repository
- many minor bugfixes

### v0.4.4
- improved GUI navigation
- support for multiple domains
- info API endpoint
- metainfo browser
- support for latest exciting version
- bugfixes in system normalization
- many minor bugfixes

### v0.4.3
- more flexible celery routing
- config via nomad.yml
- repo_db can be disabled
- publishing of calculations with failed processing
- cli for managing running processing tasks

### v0.4.2
- bugfixes regarding the migration
- better migration configurability and reproducibility
- scales to multi node kubernetes deployment
