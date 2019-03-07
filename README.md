[![pipeline status](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/master/pipeline.svg)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/master)
[![coverage report](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/master/coverage.svg)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/master)

This project tries and test approaches that might lead to an improved architecture for
nomad@FAIR.

## Getting started

Read the docs. The documentation is part of the source code. It covers aspects like
introduction, architecture, development setup/deployment, contributing, and API reference.

### Read the docs on the latest deployed version

Currently, there is only a *staging* version running at garching. There is
no real production system yet. You have to expect frequent down times and restarts.
You can access the running system and its documentation here:

[http://enc-staging-nomad.esc.rzg.mpg.de/nomad/docs](http://enc-staging-nomad.esc.rzg.mpg.de/fairdi/nomad/v0.3.0/docs)

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

### v0.4.2
- bugfixes regarding the migration
- better migration configurability and reproducibility
- scales to multi node kubernetes deployment