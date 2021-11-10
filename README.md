[![pipeline status](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/master/pipeline.svg)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/master)
[![coverage report](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/master/coverage.svg)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/master)

This project implements the new *nomad@FAIRDI* infrastructure. Contrary to its NOMAD CoE
predecessor, it implements the NOMAD Repository and NOMAD Archive functionality within
a single cohesive application. This project provides all necessary artifacts to develop,
test, deploy, and operate the NOMAD Respository and Archive, e.g. at
[https://nomad-lab.eu](https://nomad-lab.eu).

In the future, this project's aim is to integrate more NOMAD CoE components, like the NOMAD
Encyclopedia and NOMAD Analytics Toolkit, to fully integrate NOMAD with one GUI and consistent
APIs. Furthermore, this projects aims at establishing NOMAD as a distributed platform for
material science data sharing and management. This includes the on-site deployment of
NOMAD as a standalone service (*oasis*), the federated use of NOMAD through a
serious of full and partial *mirrors*, the integration of 3rd party material science
databases (i.e. [Aflow](http://www.aflow.org/), [OQMD](http://oqmd.org/),
[Materials Project](https://materialsproject.org/)), and support for open APIs and
standards like the [Optimade](http://www.optimade.org/) API.


## Getting started

### Using NOMAD as a Python package

You can install the `nomad` Python package from source distribution with pip. Please
note, that this will only install part of NOMAD's dependencies that will only allow
your to use NOMAD's client library, e.g. to access the NOMAD Archive.
```
pip install nomad-lab
```

To **use the NOMAD parsers for example**, install the `parsing` extra:
```
pip install nomad-lab[parsing]
nomad parse --show-archive <your-file-to-parse>
```

### For NOMAD developer

Read the [docs](https://nomad-lab.eu/prod/rae/docs/index.html). The documentation is also part
of the source code. It covers aspects like introduction, architecture, development setup/deployment,
contributing, and API reference.


## Change log

Omitted versions are plain bugfix releases with only minor changes and fixes.

### v0.10.9
- new AI Toolkit GUI page
- many minor parser fixes and improvements

### v0.10.7
- adding OpenMX parser

### v0.10.6
- support for NOMAD fields in optimade

### v0.10.4
- new "basic" parser to cover codes without proper parser
- removed old nomad-coe parser dependencies
- many minor parser fixes and improvements

### v0.10.3
- fixes in the VASP parser
- new turbemole parser
- property placeholders while loading entry page
- improved UI navigation with breadcrumbs

### v0.10.2
- fixes small parser and normalizer issues
- fixes broken embargo lifting
- fixes default keycloak configuration for authenticated access via ArchiveQuery

### v0.10.0
- The entries page shows visualizations for key properties of the underlying data
- A new more consistent API (/api/v1) alongside the old API (/api)
- OPTIMADE implementation based on optimade-python-tools
- Re-written parsers for VASP, FHI-aims, exciting, ABINIT, and Crystal

### v0.9.9
- A rdf-API that provides dcat datasets and catalog for NOMAD entries.
- Support to directly publish upon upload via API.

### v0.9.8
- A new library for parsing text-based raw files.
- A new main menu in the GUI.
- Upload OASIS uploads to central NOMAD.
- Updated documentation.

### v0.9.3
- Encyclopedia with dedicated materials search index.

### v0.9.0
- The encyclopedia runs on top of the new infrastructure. The GUIs are integrated via
  bi-lateral navigation between entries.
- The parsers have new documentation and development instructions and can be easily developed
  based on NOMAD's PyPi package nomad-lab
- We introduce workflow metadata (starting with geometry optimizations). This includes search
  for workflow parameter, and 3-tiered archive storage for quick access to results for analytics.
- A new GUI to browse  the Archive and Metainfo
- The Artificial Intelligence Toolkit (at least its tutorial page) is part of the GUI
- The OASIS comprises Repository, Archive, Metainfo, Encyclopedia, and Artificial Intelligence Toolkit.

### v0.8.7
- a new variant of the Metainfo browser

### v0.8.1
- switched to support Python 3.7
- client library as pypi package `nomad-lab`

### v0.8.0
- new underlying datamodel that allows to maintain multiple domains
- multiple domains supported the GUI
- new metainfo implementation
- API endpoint to access the metainfo
- new archive based on new metainfo
- client library that serves archive data as objects (with tab completion) not dictionaries
- properties and user tab in the search GUI
- improved performance on most parsers
- NOMAD source distribution

### v0.7.9
- Everything to run a simple NOMAD OASIS based on the central user-management
- minor bugfixes

### v0.7.7
- Shows dataset contents with embargo data, but hides the entry details (raw-files, archive)
- minor bugfixes

### v0.7.5
- AFLOWLIB prototypes (archive)
- primitive label search
- improved search performance based on excluded fields
- improved logs
- minor bugfixes

### v0.7.3
- fixed aborted raw-file downloads
- improved representation of data availability (staging, embargo, public) in GUI
- user data uploads ordered by upload time
- user data shows uploads with name
- configurable embargo period length
- minor bugfixes

### v0.7.2
- API curl, Python, and results for entries and search queries shown in the GUI
- minor bugfixes

### v0.7.1
- Download of archive files based on search queries
- minor bugfixes

### v0.7.0
- User metadata editing and datasets with DOIs
- Revised GUI lists (entries, grouped entries, datasets, uploads)
- Keycloak based user management
- Rawfile preview
- no dependencies with the old NOMAD CoE Repository

### v0.6.2
- GUI performance enhancements
- API /raw/query endpoint takes file pattern to further filter download contents and
  strips potential shared path prefixes for a cleaner download .zip
- Stipped common path prefixes in raw file downloads
- minor bugfixes

### v0.6.0
- GUI URL, and API endpoint that resolves NOMAD CoE legacy PIDs
- Support for datasets in the GUI
- more flexible search python module and repo API
- support for external_id
- support for code-based raw_id
- Optimade API 0.10.0
- GUI supports Optimade filter query and other quantities
- minor bugfixes

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