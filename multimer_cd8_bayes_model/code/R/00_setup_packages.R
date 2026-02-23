# ================================================================================
# R packages and other software needed to run code
# ================================================================================

# ================================================================================
# cmdstanR (R-interface to cmdstan), and cmdstan, needed to run Stan Bayesian model
# Recommended to restart R before steps below to avoid namespace issues
# ================================================================================

# Install cmdstanR interface if needed
if (!"cmdstanr" %in% installed) {
  install.packages(
    "cmdstanr",
    repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
  )
} else {
  message("cmdstanr already installed.")
}

# Check that the required C++ toolchain for CmdStan is available:
library(cmdstanr)
check_cmdstan_toolchain()

# If C++ toolchain OK ("The C++ toolchain required for CmdStan is setup properly!"):
# Install cmdstan using function from cmdstanR:
install_cmdstan(cores = 2)

# If issue with toolchain: check instructions here: 
# https://mc-stan.org/cmdstanr/articles/cmdstanr.html

# ================================================================================
# Other needed R-packages
# ================================================================================
# CRAN packages
cran_packages = c(
  "bayesplot",
  "janitor",
  "readxl",
  "tidyverse"
)

installed = rownames(installed.packages())
cran_missing = setdiff(cran_packages, installed)

if (length(cran_missing) > 0) {
  install.packages(cran_missing)
} else {
  message("All CRAN packages already installed.")
}

