# ================================================================================
# 00_setup_packages.R
# Install required R packages and CmdStanR/CmdStan for the analysis workflow
#
# Purpose:
#   - Install required CRAN packages (if missing)
#   - Install CmdStanR (if missing)
#   - Check CmdStan toolchain setup
#   - Install CmdStan (if needed)
#
# Notes:
#   - Recommended to run interactively
#   - If installation fails due to namespace/session issues, restart R and rerun
#   - See CmdStanR installation guide if toolchain setup fails:
#     https://mc-stan.org/cmdstanr/articles/cmdstanr.html
# ================================================================================

# Install cmdstanR interface if needed
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

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

