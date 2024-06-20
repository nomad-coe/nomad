# Domain-specific examples for X-ray photoelectron spectroscopy

!!! warning "Attention"
    We are currently working to update this content.

## Contextualization for the technique and the scientific domain
A variety of file formats are used in the research field of X-ray photoelectron spectroscopy and related techniques. The `pynxtools-xps` plugin of the `pynxtools` parsing library solves the challenge of how these formats can be parsed and normalized into a common representation that increases interoperability and adds semantic expressiveness.

`pynxtools-xps`, which is a plugin for [pynxtools](https://github.com/FAIRmat-NFDI/pynxtools), provides a tool for reading data from various propietary and open data formats from technology partners and the wider XPS community and standardizing it such that it is compliant with the [NeXus](https://www.nexusformat.org/) application definitions [`NXmpes`](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXmpes.html) and [`NXxps`](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXxps.html), which is an extension of `NXmpes`.

- [Documentation](https://fairmat-nfdi.github.io/pynxtools-xps/)
- [GitHub repository](https://github.com/FAIRmat-NFDI/pynxtools-xps/)
- [Issue tracker](https://github.com/FAIRmat-NFDI/pynxtools-xps/issues/)

## How to upload XPS data to NOMAD
Documentation on how to upload XPS data sets from different sources can be found [here](https://fairmat-nfdi.github.io/pynxtools-xps/tutorial/nomad.html)

## Supported file formats
A list of the supported file formats can be found in the `pynxtools-xps` [documentation](https://fairmat-nfdi.github.io/pynxtools-xps/).
