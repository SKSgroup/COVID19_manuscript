# ================================================================================
# 00_setup_packages.R
# Install required R packages and CmdStanR/CmdStan for this analysis.
#
# Notes:
#   - Run interactively
#   - Restart R and rerun if installation fails due to session/namespace issues
#   - CmdStanR setup guide: https://mc-stan.org/cmdstanr/articles/cmdstanr.html
# ================================================================================

# ================================================================================
# 1) Install CmdStanR (R interface to CmdStan)
# ================================================================================
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))

# ================================================================================
# 2) Check CmdStan toolchain setup
# If issue with toolchain: check instructions here:
# https://mc-stan.org/cmdstanr/articles/cmdstanr.html
# ============================================================================================
library(cmdstanr)
check_cmdstan_toolchain()

# ================================================================================
# 3) Install CmdStan
# ================================================================================
install_cmdstan(cores = 2)

# ================================================================================
# 4) Install other required CRAN packages
# ================================================================================
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
