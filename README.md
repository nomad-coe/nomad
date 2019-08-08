[![pipeline status](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/master/pipeline.svg)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/master)
[![coverage report](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/master/coverage.svg)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/master)

This project tries and test approaches that might lead to an improved architecture for
nomad@FAIR.

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
git submodules init --depth 1
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
