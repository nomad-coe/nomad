---
title: 'NOMAD: A distributed web-based platform for managing materials science research data'
tags:
  - Python
  - javascript
  - database
  - materials science
  - data management
  - electronic lab notebook
  - FAIR
authors:
  - name: Markus Scheidgen
    email: markus.scheidgen@physik.hu-berlin.de
    orcid: 0000-0002-8038-2277
    corresponding: true
    equal-contrib: true
    affiliation: 1
  - name: Lauri Himanen
    email: lauri.himanen@gmail.com
    orcid: 0000-0002-3130-8193
    equal-contrib: true
    affiliation: 1
  - name: Alvin Noe Ladines
    email: alvin.noe.ladines@physik.hu-berlin.de
    orcid: 0000-0003-0077-2097
    equal-contrib: true
    affiliation: 1
  - name: David Sikter
    email: david.sikter@physik.hu-berlin.de
    orcid: 0000-0002-2102-7160
    equal-contrib: true
    affiliation: 1
  - name: Mohammad Nakhaee
    email: mohammad.nakhaee@physik.hu-berlin.de
    orcid: 0000-0003-4146-0129
    equal-contrib: true
    affiliation: 1
  - name: Ádám Fekete
    email: adam@fekete.co.uk
    orcid: 0000-0002-6263-897X
    equal-contrib: true
    affiliation: 1
  - name: Theodore Chang
    email: theodore.chang@physik.hu-berlin.de
    orcid: 0000-0002-4911-0230
    equal-contrib: true
    affiliation: 1
  - name: Amir Golparvar
    email: amir.golparvar@physik.hu-berlin.de
    orcid: 0000-0002-2412-6615
    equal-contrib: true
    affiliation: 1
  - name: José A. Márquez
    email: jose.marquez@physik.hu-berlin.de
    orcid: 0000-0002-8173-2566
    affiliation: 1
  - name: Sandor Brockhauser
    email: sandor.brockhauser@physik.hu-berlin.de
    orcid: 0000-0002-9700-4803
    affiliation: 1
  - name: Sebastian Brückner
    email: sebastian.brueckner@ikz-berlin.de
    orcid: 0000-0002-5969-847X
    affiliation: 2
  - name: Luca M. Ghiringhelli
    email: luca.ghiringhelli@physik.hu-berlin.de
    orcid: 0000-0001-5099-3029
    affiliation: 1
  - name: Felix Dietrich
    email: felix.dietrich@tum.de
    orcid: 0000-0002-2906-1769
    affiliation: 3
  - name: Daniel Lehmberg
    email: d.lehmberg@tum.de
    orcid: 0000-0002-4012-5014
    affiliation: 3
  - name: Thea Denell
    email: denell@physik.hu-berlin.de
    orcid: 0009-0000-7185-9363
    affiliation: 1
  - name: Andrea Albino
    email: andrea.albino@physik.hu-berlin.de
    orcid: 0000-0001-9280-7431
    affiliation: 1
  - name: Hampus Näsström
    email: hampus.naesstroem@physik.hu-berlin.de
    orcid: 0000-0002-3264-1692
    affiliation: 1
  - name: Sherjeel Shabih
    email: sherjeel.shabih@physik.hu-berlin.de
    orcid: 0009-0008-6635-4465
    affiliation: 1
  - name: Florian Dobener
    email: florian.dobener@physik.hu-berlin.de
    orcid: 0000-0003-1987-6224
    affiliation: 1
  - name: Markus Kühbach
    email: markus.kuehbach@physik.hu-berlin.de
    orcid: 0000-0002-7117-5196
    affiliation: 1
  - name: Rubel Mozumder
    email: mozumder@physik.hu-berlin.de
    orcid: 0009-0007-5926-6646
    affiliation: 1
  - name: Joseph F. Rudzinski
    email: joseph.rudzinski@physik.hu-berlin.de
    orcid: 0000-0003-3403-640X
    affiliation: 1
  - name: Nathan Daelman
    email: ndaelman@physik.hu-berlin.de
    orcid: 0000-0002-7647-1816
    affiliation: 1
  - name: José M. Pizarro
    email: jose.pizarro@physik.hu-berlin.de
    orcid: 0000-0002-6751-8192
    affiliation: 1
  - name: Martin Kuban
    email: kuban@physik.hu-berlin.de
    orcid: 0000-0002-1619-2460
    affiliation: 1
  - name: Cuauhtemoc Salazar
    email: cuauhtemoc.salazar@physik.hu-berlin.de
    orcid: 0000-0002-9635-5062
    affiliation: 1
  - name: Pavel Ondračka
    email: ondracka@mail.muni.cz
    orcid: 0000-0003-0729-629X
    affiliation: 4
  - name: Hans-Joachim Bungartz
    email: bungartz@cit.tum.de
    orcid: 0000-0002-0171-0712
    affiliation: 3
  - name: Claudia Draxl
    email: claudia.draxl@physik.hu-berlin.de
    orcid: 0000-0003-3523-6657
    affiliation: 1
affiliations:
  - name: Department of Physics, Humboldt-Universität zu Berlin, Berlin, Germany
    index: 1
  - name: Leibniz Institut für Kristallzüchtung, Berlin, Germany
    index: 2
  - name: School of Computation, Information and Technology, Technical University of Munich, Munich, Germany
    index: 3
  - name: Department of Plasma Physics and Technology, Masaryk University, Brno, Czech Republic
    index: 4
date: 20 March 2023
bibliography: paper.bib
---

# Summary

Materials science research is becoming increasingly data-driven, which requires more effort to manage, share, and publish data.
NOMAD is a web-based application that provides data management for materials science research data.
In addition to core data management functions like uploading and sharing files, NOMAD allows for entering structured data using customizable forms providing the software with electronic laboratory notebook (ELN) functionalities.
It automatically extracts rich meta-data from supported file formats, normalizes and converts data from these formats, and provides a faceted search with materials science-specific filters based on extracted meta-data.
NOMAD integrates data analysis and machine learning tools.
Installations of NOMAD can be connected to share data between research institutes and can publish data to an open central NOMAD service.
The NOMAD software is distributed as a Docker image to create data management services and as a Python package to automate the client's use of these services.

# Statement of need

In materials science, researchers use many methods, instruments, tools, and workflows
to produce large volumes of heterogeneous data artifacts. The contained data often
describes related research objects (materials, samples, or properties) and it is believed
that all combined data hold great potential for data re-use and machine learning [@scheffler:2022; @sbailo:2022].
This is clearly being acknowledged not only by the research community but also by funding agencies,
which are increasingly demanding coordinated efforts in availability and longevity of open data by preserving
and documenting all produced research data and meta-data.

While individual researchers struggle with organizing and analyzing more and more data
artifacts, communities face new challenges in making data findable, accessible, inter-operable,
and reproducible (FAIR) [@wilkinson:2016]. A key factor to FAIR data is to combine data with meta-data
and to put all data into machine and human comprehensible representations [@ghiringhelli:2017; @ghiringhelli:2023].

Materials scientists require effective solutions for managing their research data, but they should not have to develop their own individual solutions. Hence, there is great demand in services (and software to run such services) that provide the mentioned features and make data FAIR.
This is evident in the great number of published datasets on services like [NOMAD](https://nomad-lab.eu) [@draxl:2018] (the main deployment of the NOMAD software), and an increasing number of materials science databases that all (re-)implement very similar functionality to publish their data.

NOMAD addresses these needs in two ways. First, NOMAD improves the data-driven workflows of individuals
and small labs by formalizing data acquisition, organizing and sharing data, homogenizing and normalizing data
for analysis, and integrating with analysis tools.
This way, NOMAD provides the incentives and tools for research individuals
to put the necessary efforts into preparing FAIR (meta-)data. Secondly, NOMAD allows for sharing and
publishing prepared data and can be used by communities as a repository for FAIR data.

# Usage of NOMAD and related software

The NOMAD software is used to operate a public and free NOMAD service
that allows everyone to share and publish materials science research data (https://nomad-lab.eu).
This public NOMAD service contains over 12 million individual materials science simulations and around 50 thousand entries describing materials experiments or synthesis. NOMAD is publicly available since 2014 and includes data from
over 500 international authors.

The NOMAD software can also be independently hosted by universities and other institutions when
the use of the central service is not possible. Such self-managed installations are called NOMAD
Oases to distinguish them from the public NOMAD service. A NOMAD Oasis might be required when an
institution needs to significantly customize the software for a specific need, the data volumes
are too large to be conveniently transferred over the public internet, or when there are concerns
about privacy or security. It should be noted that there is the possibility to transfer data
between different installations, and in order to adhere to the FAIR principles, the data (or
at least meta-data) in these Oases would ideally be transferred to the public NOMAD service.
NOMAD Oasis is used by an increasing number of research institutes. NOMAD Oasis can be used freely
as per our OSI license following the instruction in the [NOMAD documentation](https://nomad-lab.eu/docs).

The [NFDI consortium FAIRmat](https://www.fairmat-nfdi.eu/fairmat)
uses NOMAD software as the basis for its federated FAIR data infrastructure [@scheffler:2022].

[OPTIMADE](https://www.optimade.org/) [@andersen:2021] is an API specification (with associated
software implementation) for materials science databases. NOMAD provides an implementation
of the OPTIMADE specification and is an active member of the OPTIMADE consortium.

Other materials science databases (and the respective software) focus on publishing data that were
produced with a specific framework and carefully curated by the group behind the database.
Typical examples are databases of
high-throughput simulations that try to systematically explore theoretical materials. Three
of the larger databases of this kind are the [Materials Project](https://materialsproject.org/)
[@jain:2013], [AFLOW](https://aflow.org/) [@curtarolo:2012], and [OQMD](https://oqmd.org/) [@saal:2013].
The raw data of these databases have also been published on NOMAD. The project [AiiDA](https://aiida.net/) [@huber:2020]
allows scientists to design and run simulation workflows. AiiDA data can be published
to AiiDA's [materialscloud](https://www.materialscloud.org/). There are also examples for experimental
materials science databases, e.g. [HTEM](https://data.nrel.gov/submissions/75) [@zakutayev:2018].

NOMAD relies on many open-source packages; a few more notable ones from the materials
science domains are: *MatID*, a software package to identify material structure system types and symmetries [@himanen:2018],
*ASE*, a software package to manipulate material structures in Python [@larsen:2017],
pymatgen, an open-source python library for materials analysis [@ong:2013], and *NeXus*,
a file-format standard, schemas, and tools for experimental materials science data [@konnecke:2015].

# Acknowledgements

NOMAD software development is funded by the German National Research Data
 Infrastructure (NDFI) consortia FAIRmat (Deutsche Forschungsgemeinschaft DFG, 460197019) and the NOMAD CoE (EU Horizon 2020, 951786), 
previous financial support was provided by the NOMAD CoE (EU Horizon 2020, 676580) and the Max-Planck Network BigMax. 
The Max Planck Computing and Data Facility (MPCDF) is hosting NOMAD's GitHub and operating the public NOMAD service. 

# References
