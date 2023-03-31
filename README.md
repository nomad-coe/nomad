[![pipeline status](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/develop/pipeline.svg)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/develop)
[![backend coverage report](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/develop/coverage.svg?job=python+tests&key_text=backend+coverage&key_width=130)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/develop)
[![frontend coverage report](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/badges/develop/coverage.svg?job=gui+tests&key_text=frontend+coverage&key_width=130)](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/commits/develop)

**NOMAD** is a web-based research data management software for materials science.
You find the official project homepage and documentation here [https://nomad-lab.eu](https://nomad-lab.eu).
NOMAD is used to provide an open service for managing and publish research data of the same name.
On-premise installations of NOMAD (Oasis) allow research groups to locally manage data with
customized NOMAD version and with their own compute and storage resources.

## Contributing

There are two forks of this repository, one on **GitHUB** and one on MPCDF's **GitLab**.

NOMAD's [GitHUB project](https://github.com/nomad-coe/nomad) always contains the current `develop` branch.
It can be used to report issues, fork the project, and to create pull requests. After review, pull requests
will be pushed to the GitLab project and merged there. Use the regular GitHUB flow to contribute
as an external developer here.

NOMAD's [GitLab project](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR) at [MPCDF](https://www.mpcdf.mpg.de/)
is used for the main development activities. It runs all CI/CD pipelines and official deployments.
It is openly readable, but requires an MPCDF account for active contributions. If you
are a member of the FAIRmat or NOMAD CoE project, contribute here.

Most sub-modules, e.g. NOMAD's parsers, are hosted in individual projects on GitHUB
within the [nomad-coe organization](https://github.com/nomad-coe).

## Getting started

For a general project overview visit the official project page [https://nomad-lab.eu](https://nomad-lab.eu). For specific use of the NOMAD software follow these links to our documentation:

- [get started as a developer](https://nomad-lab.eu/prod/v1/docs/develop/setup.html)
- [install and use NOMAD as Python package (to use our APIs or parsers)](https://nomad-lab.eu/prod/v1/docs/pythonlib.html)
- [install NOMAD Oasis](https://nomad-lab.eu/prod/v1/docs/oasis.html)

## Change log

Omitted versions are plain bugfix releases with only minor changes and fixes. The
file [`CHANGELOG.md`](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR/-/blob/develop/CHANGELOG.md)
contains much more detailed information about changes and fixes in the released versions.

### v1.1.9
- ELN improvements for safer file handling
- many smaller fixes and changes

### v1.1.8
- updated app branding to fit the new http://nomad-lab.eu
- use Python 3.9
- Pydantic models for config

### v1.1.7
- refactored NORTH configuration and k8s deployment
- categorized build-in ELN schemas
- categorized explore menu
- categorized example uploads

### v1.1.6
- application dashboards with solar cell as an example
- detailed selection of schemas and reference targets in ELN GUI elements
- reference navigation in archive browser sections and on the entry overview
- new workflow model and updated workflow parsing/normalization for computational codes
- examples from FAIRmats experimental material science area
- search for custom quantities
- initial implementation of bundle import/export from CLI
- formalized and documented configuration and metainfo annotations
- new build system for docker image and Python package

### v1.1.0
- example uploads
- custom schema support
- ELN functionality
- nexus schema
- north hub
- file browser

### v1.0.6
- upgraded to Elasticsearch 7.x

### v1.0.4
- tabular data schema

### v1.0.3
- refactored DCAT to use fast api, added DOIs
- refactored ArchiveQuery client
- documentation and fixes for Oasis with keycloak
- many minor GUI bugfixes

### v1.0.0
- new search interface
- new v1 API (entries, materials, upload, datasets, sync)
- refactored metainfo and parsers
- new upload UI and incremental uploads

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
