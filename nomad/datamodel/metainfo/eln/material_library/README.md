# Material Library

This an example ELN schema inspired by a lab working with thin-file material library samples.
The schema has three top-level entry types (i.e. sections that extend `EntryData`):
sample, chemical, instrument. It shows some of NOMAD's ELN capabilities:

- References to files: Chemicals can have safety sheet PDFs, Instruments refer to calibration documentation.
- References between different entry types: Sample processes can refer to chemicals and
instruments.
- Custom schema specific processing of referenced files: The PVD Evaporation process
can reference a `.csv` file that is parsed and populates values of the respective section.
The XRD measurement can reference a `.xrdml` file which is parsed with the xrdtools library.
- A boolean flag on processes `creates_layer` automatically creates a new layer sub-section
in the sample and links it back the the process that *created* it.


