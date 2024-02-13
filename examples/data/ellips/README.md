# Introduction

This example presents the capabilities of the NOMAD platform to store and standardize ellipsometry data. It shows the generation of a NeXus file according to the [NXellipsometry](https://manual.nexusformat.org/classes/contributed_definitions/NXellipsometry.html#nxellipsometry) application definition and a successive analysis of an example data set (measured values of the ellipsometric angles Psi and Delta of SiO2 on Si).

# Viewing uploaded data

Below, you find an overview of your uploaded data.
Click on the `> /` button to get a list of your data or select **FILES** from the top menu of this upload.
You may add your own files to the upload or experiment with the pre-existing electronic lab notebook (ELN) example.
The ELN follows the general structure of NOMAD ELN templates and you may refer to the [documentation](https://nomad-lab.eu/prod/v1/staging/docs/archive.html) or a [YouTube tutorial](https://youtu.be/o5ETHmGmnaI) (~1h)
for further information.
When the ELN is saved a NeXus file will be generated from the provided example data.
You may also view your supplied or generated NeXus files here with the H5Web viewer.
To do so open the **FILES** tab and just select a `.nxs` file.

## Filelist

- Measurement data: test-data.dat
- Metadata file: eln_data.yaml
- Notebook with instructions: Ellipsometry workflow example.ipynb
- NeXus file: SiO2onSi.ellips.nxs (will be created when running the notebook)

# Analyzing the data

The examples work through the use of NOMAD remote tools hub (NORTH) containers, i.e. besides using and dealing with the uploaded ellipsometry data, an analysis container can be started. If you want to execute the examples locally you may also use your local python and jupyterlab installation. Please refer to the documentation of [pynxtools](https://github.com/FAIRmat-NFDI/pynxtools.git), analysis tool [pyElli](https://github.com/PyEllips/pyElli) and [h5web](https://github.com/silx-kit/h5web) on how to install it on your machine.

To start an analysis, note your upload id (which you find on top of this explanation) and select **ANALYZE** from the top menu, then **NOMAD Remote Tools Hub**.
In the appearing list you'll find the `ellips` container, click on it and click **LAUNCH**.
After a few moments a new tab will open which displays a jupyter environment providing the required analysis tools.
To find the examples navigate to uploads inside the jupyter hub browser and select the folder with your noted upload id.
There you'll find the example `ipynb` notebook.
Double-clicking the notebook will open the example in the jupyter main window.

# Where to go from here?

If you're interested in using this pipeline and NOMAD in general you'll find support at [FAIRmat](https://www.fairmat-nfdi.eu/fairmat/consortium).

For questions regarding the experiment or this specific example contact [Carola Emminger](https://www.fairmat-nfdi.eu/fairmat/fairmat_/fairmatteam).

If you want to learn more about the analysis tool, please refer to its [github page](https://github.com/PyEllips/pyElli), where you may raise an [issue](https://github.com/PyEllips/pyElli/issues), look at the [documentation](https://pyelli.readthedocs.io/en/latest/) or just get in contact with the developers.
