# Introduction

This example presents the capabilities of the NOMAD platform to store and standardize multi photoemission spectroscopy (MPES) experimental data. It contains three major examples:

- Taking a pre-binned file, here stored in a h5 file, and converting it into the standardized MPES NeXus format. There exists a [NeXus application definition for MPES](https://manual.nexusformat.org/classes/contributed_definitions/NXmpes.html#nxmpes) which details the internal structure of such a file.
- Binning of raw data (see [here](https://www.nature.com/articles/s41597-020-00769-8) for additional resources) into a h5 file and consecutively generating a NeXus file from it.
- An analysis example using data in the NeXus format and employing the [pyARPES](https://github.com/chstan/arpes) analysis tool to reproduce the main findings of [this paper](https://arxiv.org/pdf/2107.07158.pdf).

# Viewing uploaded data

Below, you find an overview of your uploaded data.
Click on the `> /` button to get a list of your data or select **FILES** from the top menu of this upload.
You may add your own files to the upload or experiment with the pre-existing electronic lab book example.
The ELN follows the general structure of NOMAD ELN templates and you may refer to the [documentation](https://nomad-lab.eu/prod/v1/staging/docs/archive.html) or a [YouTube tutorial](https://youtu.be/o5ETHmGmnaI) (~1h)
for further information.
When the ELN is saved a NeXus file will be generated from the provided example data.
You may also view your supplied or generated NeXus files here with the H5Web viewer.
To do so open the **FILES** tab and just select a `.nxs` file.

# Analysing the data

The examples work through the use of NOMAD remote tools hub (NORTH) containers, i.e. besides using and dealing with the uploaded MPES data, an analysis container can be started. If you want to execute the examples locally, please refer to the `INSTALL.md` file in this directory.
Please note that the binning from raw files is rather memory intense and you should have at least 40 GB of available
memory to be able to execute this particular example.
If your local container suddenly stops and shows an error 137 it is due to memory limitations.

In the container you'll find three example notebooks containing examples mentioned above.

To start an analysis, note your upload id (which you find on top of this explanation) and select **ANALYZE** from the top menu, then **NOMAD Remote Tools Hub**.
In the appearing list you'll find the `mpes` Example, click on it and click **LAUNCH**.
After a few moments a new tab will open which displays a jupyter environment providing the required analysis tools.
To find the examples navigate to uploads inside the jupyter hub browser and select the folder with your noted upload id.
There you'll find the example `ipynb` notebooks.
Double-clicking one of the notebooks will open the example in the jupyter main window.
From here you find detailed instructions inside each of the notebooks.

# Where to go from here?

If you're interested in using this pipeline and NOMAD in general you'll find support at [FAIRmat](https://www.fairmat-nfdi.eu/fairmat/consortium).

If you have questions regarding the experiments in the examples please refer to the [Dynamics of Correlated Materials](https://pc.fhi-berlin.mpg.de/docm/) or [Structural & Electronic Surface Dynamics](https://pc.fhi-berlin.mpg.de/sesd/) working group at the Fritz-Haber Institute Berlin.

For general questions regarding the MPES pipeline and if you're interested in building one for your
own research workflow you may contact [Florian Dobener](https://www.fairmat-nfdi.eu/fairmat/fairmat_/fairmatteam) from the FAIRmat consortium.
