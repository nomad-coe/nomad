# Domain-specific examples for STS / STM (scanning tunneling spectroscopy / microscopy)


Build upon your understanding of NOMAD's features with domain-specific examples and explanations.

### Contextualization for the technique and the scientific domain
A variety of file formats (coming technology instrumens) are used in the research field of scanning tunneling microscopy (STM) and scanning tunneling spectroscopy (STS) to investigate topological propertise of surface of subjected material. The [pynxtools-stm](https://github.com/FAIRmat-NFDI/pynxtools-stm) plugin (note: The single plugin handles both STM as well STS techniques) of the [pynxtools](https://github.com/FAIRmat-NFDI/pynxtools) parsing library solves the challenges of how these formats can be parsed and normalized into a common representation that increases interoperability and adds semantic expressiveness.

The [pynxtools-stm](https://fairmat-nfdi.github.io/pynxtools-stm/) provides a indispensable tool to transfor the STM as well as STS experimental data (sometime refered as raw data or machine data) to common standarized structure defined in [NXsts](https://fairmat-nfdi.github.io/nexus_definitions/classes/contributed_definitions/NXsts.html#nxsts) ([GitHub page](https://github.com/FAIRmat-NFDI/nexus_definitions/blob/fairmat/contributed_definitions/NXsts.nxdl.xml)) application definition build with the help of [NeXus ontology](https://www.nexusformat.org/) ([GitHub page](https://github.com/FAIRmat-NFDI/nexus_definitions/tree/fairmat)). One of the main goals of such effort is to make the data comming from diverse sources comparable, searchable and shearable under the hood of NONAD research data management platform.

For full benefits and usages of the reader please following links:

- [Full Reader Documentation](https://fairmat-nfdi.github.io/pynxtools-stm/)
- [GitHub Repository](https://github.com/FAIRmat-NFDI/pynxtools-stm)
- [Issue Tracker](https://github.com/FAIRmat-NFDI/pynxtools-stm/issues)


## How to upload XPS data to NOMAD
Documentation on how to upload STM / STS data sets from different sources can be found [here](https://fairmat-nfdi.github.io/pynxtools-sts/tutorials/nomad-tutorial.html)

## Supported file formats
A list of the supported file formats can be found in the `pynxtools-stm` [documentation](https://fairmat-nfdi.github.io/pynxtools-stm/).
