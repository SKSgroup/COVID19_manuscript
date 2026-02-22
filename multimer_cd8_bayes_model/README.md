# Bayesian model analysis of T-cell response to SARS-CoV-2 

This subdirectory contains code and data for Bayesian modeling of epitope-specific CD8 T-cell responses measured by pHLA multimer assays.

The workflow includes:

1. processing raw Excel input files into analysis-ready tables and `.rds` objects,
2. constructing a Stan data list (`standata`) for model fitting,
3. fitting a mixed-effects logistic regression model in Stan,
4. using fitted model outputs to generate figures for the manuscript.

## Directory structure

- `code/R/`  
  R scripts for data processing, model preparation, fitting, and figure generation.

- `code/Stan/`  
  Stan model code.

- `data/raw/`  
  Raw input files (Excel).

- `data/processed/`  
  Processed intermediate `.rds` objects and Stan input data (`standata`).

- `results/`  
  Model fit outputs and other generated analysis outputs.

- `figures/`  
  Generated figures for the manuscript.

## Running the code (RStudio project / working directory)

This subdirectory contains an RStudio project file: `multimer_cd8_bayes_model.Rproj`.

Please open this `.Rproj` file in RStudio before running the scripts.  
This ensures that the **working directory is the root of this subdirectory** (`multimer_cd8_bayes_model/`), so that relative paths such as:

- `data/raw/...`
- `data/processed/...`
- `code/Stan/...`

resolve correctly.

If you are not using RStudio, set the working directory manually to the `multimer_cd8_bayes_model/` directory before running scripts.

## Typical workflow (example order)

The scripts in `code/R/` are intended to be run in sequence (names may change as the workflow evolves). For example:

1. `01_process_raw_data.R`  
   Read raw Excel files, clean/format variables, create processed `.rds` objects, and create `standata`.

2. *(future scripts)* `02_...`, `03_...`, etc.  
   Model fitting and figure generation.

## Software requirements

TBD

