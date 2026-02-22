## Repository structure

This repository contains separate sub-directories for different analyses

### `multimer_cd8_bayes_model/`
Code and data for Bayesian modeling of epitope-specific CD8 T-cell responses measured by pHLA multimer assays. This subdirectory includes:

- raw input data files (Excel)
- processed intermediate data objects (`.rds`)
- R scripts for data processing and figure generation
- Stan model code for fitting the mixed-effects logistic regression model
- output directories for model results and figures