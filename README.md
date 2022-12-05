# TIDES

Implementation of the contribution of non-coplanarity into the tidal evolution of [MESA](https://docs.mesastar.org/en/r15140/)

This code contains the tidal evolution as described in [Repetto & Nelemans, 2014](https://ui.adsabs.harvard.edu/abs/2014MNRAS.444..542R/abstract)
by using the `use_other_sync_spin_to_orbit` and `use_other_edot_tidal` flags of the **MESAbinary** module

The structure of this research project follows the philosophy of [Cookiecutter Data Science](https://github.com/drivendata/cookiecutter-data-science).

Data inside `data/raw` folder can be shared on request.

## Repository structure

```
.
├── README.md
├── config
├── data
│   └── raw
├── notebooks
├── reports
│   ├── paper
│   └── figures
├── reports
│   └── paper
└── src
    ├── models
    │   └── mesa-tides
```

The `config` folder contains configuration files for jupyter notebooks. Also, the python environment used is located in the `config` folder to
reproduce all figures shown.

Inside the `data` folder the output of our simulations is located. To maintain the storage of the repository in a low level, these raw output is not
contained within this repository, but can be replicated with the source code found in other folder

Analysis of simulations and comparisons with different codes are located in the `notebooks` folder

The modifications to the **MESA**  source code can be found in the `src` folder. All these modifications use the `run_binary_extras` module

The `reports` folder has the source files to compile the PDF of the publication
