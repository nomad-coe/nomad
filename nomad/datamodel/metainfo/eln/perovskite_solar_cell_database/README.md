# Perovskite solar cells database

This is nomad like schema of the [perovskite database solar cell](https://www.perovskitedatabase.com/). This project was first described in a paper published in Nature Energy in 2021 with the title:

An open-access database and analysis tool for perovskite solar cells based on the FAIR data principles by Jacobsson et al.
The paper is published open sourced and available at: [paper link](https://www.nature.com/articles/s41560-021-00941-3).

The NOMAD schema contains annotations so it can be used as an ELN to manage your local data.

# The perovskite database NOMAD metainfo

The Nomad metainfo of the perovskite database is at the moment designed replicating the description of the original [database](https://www.perovskitedatabase.com/Resources). The quantities defined in the sections have a description that conincides with the ones of the original database, but this might evolve in the future in NOMAD. For example, the original peorvskite compositions were given with the cations abbreviated. We have parsed most of these formulas containing abbreviations to chemical formulas so they are serachable by elements and by their formula.
Thanks *Jinzhao Li* for help translating the abbreviations!

# ELN functionality and extending the Perovskite DB metainfo

Aditional quantities have been included to be able to use aditional functionalities NOMAD and augemnt in the future the databsse. For example in the section jv, there is a quantity called `data_file` which can be used to deposit a jv file which can be parsed (`jv_parser.py`) to fill automatically many of the other quantities of the `jv` section. We note that this parser is very specific for an particular JV file.

In the section `eqe` one can drop an EQE file and several quantities get populated. The `eqe_parser.py` parser is adapted by the original matlab code 'SCRIPT_Vocrad_EQE_fit_Urbachtail' published by [Krückemeier et al.](https://doi.org/10.1002/aenm.201902573), and later translated to python by [Christian Wolff](https://github.com/chrmwo/VOC-loss-analysis). It looks for the EQE file and does some routines to extract quantities like the `bandgap_eqe` or the `urbach_energy`.

